from __future__ import annotations

from dataclasses import replace
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from ads_benchmark.catalog import build_catalog
from ads_benchmark.framework.config import load_profiles
from ads_benchmark.framework.errors import DuplicateCaseError, ValidationError
from ads_benchmark.framework.filtering import CaseFilters
from ads_benchmark.framework.model import RepositoryState
from ads_benchmark.framework.planner import Planner, case_identity


BENCHMARKING_ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIRECTORY = BENCHMARKING_ROOT / "configs"


class PlannerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.catalog = build_catalog()
        self.profiles = load_profiles(CONFIG_DIRECTORY)
        self.planner = Planner(self.profiles, self.catalog)

    def test_registered_profiles_and_exact_temporal_matrix(self) -> None:
        self.assertEqual(
            self.profiles.names(),
            ("cluster-scaling", "local-scaling", "smoke", "temporal-full"),
        )
        plan = self.planner.plan("temporal-full")
        self.assertEqual(len(plan.cases), 792)
        self.assertEqual(len({case.case_id for case in plan.cases}), 792)
        self.assertEqual(
            {case.spec.problem for case in plan.cases},
            {"igrm_l2", "igrm_heat", "pure_diffusion_igrm"},
        )
        self.assertEqual({case.spec.scheme for case in plan.cases}, {"dg", "pr", "be"})
        self.assertEqual(
            {case.spec.time.steps for case in plan.cases},
            {4, 8, 16, 32, 64, 128, 256, 512},
        )
        self.assertEqual(len({(case.spec.test_degree, case.spec.trial_degree) for case in plan.cases}), 11)

    def test_every_registered_profile_expands_to_a_nonempty_valid_plan(self) -> None:
        expected_counts = {
            "smoke": 9,
            "temporal-full": 792,
            "local-scaling": 36,
            "cluster-scaling": 36,
        }
        for profile, expected_count in expected_counts.items():
            with self.subTest(profile=profile):
                self.assertEqual(len(self.planner.plan(profile).cases), expected_count)

    def test_identity_is_stable_under_axis_and_json_key_reordering(self) -> None:
        original = self.planner.plan("temporal-full")
        document = json.loads(
            (CONFIG_DIRECTORY / "temporal-full.json").read_text(encoding="utf-8")
        )
        reordered: dict[str, object] = {}
        for key in reversed(list(document)):
            value = document[key]
            if isinstance(value, list):
                value = list(reversed(value))
            reordered[key] = value
        with tempfile.TemporaryDirectory(prefix="ads-profile-order-") as temporary:
            path = Path(temporary) / "temporal-full.json"
            path.write_text(json.dumps(reordered), encoding="utf-8")
            alternate = Planner(load_profiles(path.parent), self.catalog).plan(
                "temporal-full"
            )
        self.assertEqual(
            [case.case_id for case in original.cases],
            [case.case_id for case in alternate.cases],
        )
        self.assertEqual(original.config_hash, alternate.config_hash)

    def test_identity_excludes_profile_and_provenance_but_tracks_semantics(self) -> None:
        case = self.planner.plan("smoke").cases[0]
        original_id, _ = case_identity(case.spec)
        changed_id, _ = case_identity(replace(case.spec, openmp_threads=2))
        self.assertNotEqual(original_id, changed_id)
        clean = self.planner.plan("smoke").manifest(
            run_id="one",
            repository=RepositoryState(commit="a" * 40, dirty=False),
            created_at="2026-01-01T00:00:00+00:00",
        )
        dirty = self.planner.plan("smoke").manifest(
            run_id="two",
            repository=RepositoryState(commit="b" * 40, dirty=True),
            created_at="2027-01-01T00:00:00+00:00",
        )
        self.assertEqual(
            [item["case_id"] for item in clean["cases"]],
            [item["case_id"] for item in dirty["cases"]],
        )

    def test_identifiers_do_not_depend_on_python_hash_seed(self) -> None:
        script = (
            "from pathlib import Path;"
            "from ads_benchmark.catalog import build_catalog;"
            "from ads_benchmark.framework.config import load_profiles;"
            "from ads_benchmark.framework.planner import Planner;"
            f"p=Planner(load_profiles(Path({str(CONFIG_DIRECTORY)!r})),build_catalog());"
            "print('\\n'.join(c.case_id for c in p.plan('smoke').cases))"
        )
        outputs = []
        for seed in ("1", "987654"):
            environment = os.environ.copy()
            environment["PYTHONHASHSEED"] = seed
            environment["PYTHONPATH"] = str(BENCHMARKING_ROOT)
            completed = subprocess.run(
                [sys.executable, "-c", script],
                check=True,
                capture_output=True,
                text=True,
                env=environment,
            )
            outputs.append(completed.stdout)
        self.assertEqual(outputs[0], outputs[1])

    def test_filters_are_conjunctive_and_zero_matches_fail(self) -> None:
        cases = (
            (CaseFilters(problems=frozenset({"igrm_heat"})), 264),
            (CaseFilters(schemes=frozenset({"be"})), 264),
            (
                CaseFilters(
                    degree_pairs=frozenset({((4, 4, 4), (3, 3, 3))})
                ),
                72,
            ),
            (CaseFilters(meshes=frozenset({(4, 4, 4)})), 792),
            (CaseFilters(mpi_grids=frozenset({(1, 1, 1)})), 792),
            (CaseFilters(mpi_ranks=frozenset({1})), 792),
            (CaseFilters(openmp_threads=frozenset({1})), 792),
            (
                CaseFilters(
                    problems=frozenset({"igrm_l2"}),
                    schemes=frozenset({"pr"}),
                    degree_pairs=frozenset({((5, 5, 5), (4, 4, 4))}),
                ),
                8,
            ),
        )
        for filters, expected in cases:
            with self.subTest(filters=filters):
                self.assertEqual(
                    len(self.planner.plan("temporal-full", filters).cases), expected
                )
        with self.assertRaisesRegex(ValidationError, "selected no cases"):
            self.planner.plan(
                "temporal-full", CaseFilters(meshes=frozenset({(99, 99, 99)}))
            )

    def test_duplicate_expansion_fails_before_filtering(self) -> None:
        document = json.loads(
            (CONFIG_DIRECTORY / "smoke.json").read_text(encoding="utf-8")
        )
        document["problems"].append(document["problems"][0])
        with tempfile.TemporaryDirectory(prefix="ads-profile-duplicate-") as temporary:
            path = Path(temporary) / "duplicate.json"
            path.write_text(json.dumps(document), encoding="utf-8")
            planner = Planner(load_profiles(path.parent), self.catalog)
            with self.assertRaisesRegex(DuplicateCaseError, "duplicate case"):
                planner.plan(
                    "smoke", CaseFilters(problems=frozenset({"igrm_heat"}))
                )


if __name__ == "__main__":
    unittest.main(verbosity=2)
