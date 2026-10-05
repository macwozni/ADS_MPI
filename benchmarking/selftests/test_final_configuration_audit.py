from __future__ import annotations

from itertools import product
from pathlib import Path
import tempfile
import unittest

from ads_benchmark.catalog import build_catalog
from ads_benchmark.framework.config import load_profiles
from ads_benchmark.framework.planner import Planner, case_identity
from ads_benchmark.framework.storage import ResultStore


BENCHMARKING_ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIRECTORY = BENCHMARKING_ROOT / "configs"

PROBLEMS = frozenset({"igrm_l2", "igrm_heat", "pure_diffusion_igrm"})
SCHEMES = frozenset({"dg", "pr", "be"})
PROBLEM_SCHEMES = frozenset(product(PROBLEMS, SCHEMES))
OPENMP_LEVELS = frozenset({1, 2, 4, 8})

TEMPORAL_TIME_STEPS = frozenset(
    {
        "0.025",
        "0.0125",
        "0.00625",
        "0.003125",
        "0.0015625",
        "0.00078125",
        "0.000390625",
        "0.0001953125",
    }
)
TEMPORAL_DEGREE_PAIRS = frozenset(
    [((degree + 1,) * 3, (degree,) * 3) for degree in range(3, 9)]
    + [((degree + 2,) * 3, (degree,) * 3) for degree in range(3, 8)]
)
SCALING_DEGREE_PAIRS = frozenset(
    [((degree + 1,) * 3, (degree,) * 3) for degree in range(1, 9)]
    + [((degree + 2,) * 3, (degree,) * 3) for degree in range(1, 8)]
)


class FinalConfigurationAuditTests(unittest.TestCase):
    """Audit the complete temporal, strong, and weak plans without execution."""

    @classmethod
    def setUpClass(cls) -> None:
        planner = Planner(load_profiles(CONFIG_DIRECTORY), build_catalog())
        cls.plans = {
            name: planner.plan(name)
            for name in (
                "temporal-full",
                "cluster-scaling",
                "cluster-weak-scaling",
            )
        }

    def test_full_profiles_cover_all_problems_schemes_and_families(self) -> None:
        expected_families = {
            "temporal-full": "temporal",
            "cluster-scaling": "strong",
            "cluster-weak-scaling": "weak",
        }
        expected_case_counts = {
            "temporal-full": 792,
            "cluster-scaling": 12_960,
            "cluster-weak-scaling": 14_985,
        }
        for profile_name, expected_family in expected_families.items():
            with self.subTest(profile=profile_name):
                cases = self.plans[profile_name].cases
                self.assertEqual(len(cases), expected_case_counts[profile_name])
                self.assertEqual(
                    {(case.spec.problem, case.spec.scheme) for case in cases},
                    PROBLEM_SCHEMES,
                )
                self.assertEqual(
                    {case.spec.family for case in cases}, {expected_family}
                )

        self.assertEqual(
            set(expected_families.values()), {"temporal", "strong", "weak"}
        )

    def test_temporal_and_scaling_levels_are_complete(self) -> None:
        temporal = self.plans["temporal-full"].cases
        self.assertEqual(
            {case.spec.time.time_step for case in temporal},
            TEMPORAL_TIME_STEPS,
        )
        self.assertEqual(
            {
                (case.spec.test_degree, case.spec.trial_degree)
                for case in temporal
            },
            TEMPORAL_DEGREE_PAIRS,
        )
        self.assertEqual(len(TEMPORAL_TIME_STEPS), 8)
        self.assertEqual(len(TEMPORAL_DEGREE_PAIRS), 11)

        for profile_name in ("cluster-scaling", "cluster-weak-scaling"):
            with self.subTest(profile=profile_name):
                cases = self.plans[profile_name].cases
                self.assertEqual(
                    {
                        (case.spec.test_degree, case.spec.trial_degree)
                        for case in cases
                    },
                    SCALING_DEGREE_PAIRS,
                )
                self.assertEqual(
                    {case.spec.openmp_threads for case in cases},
                    OPENMP_LEVELS,
                )
                grids = {case.spec.mpi.process_grid for case in cases}
                for axis, label in enumerate(("X", "Y", "Z")):
                    self.assertTrue(
                        any(grid[axis] > 1 for grid in grids),
                        f"{profile_name} must scale along {label}",
                    )
        self.assertEqual(len(SCALING_DEGREE_PAIRS), 15)

        strong_grids = {
            case.spec.mpi.process_grid
            for case in self.plans["cluster-scaling"].cases
        }
        self.assertTrue(
            {(2, 1, 1), (1, 2, 1), (1, 1, 2)} <= strong_grids,
            "cluster strong scaling must cover X, Y, and Z decompositions",
        )

    def test_case_ids_are_unique_and_bind_one_semantic_configuration(self) -> None:
        cases = tuple(
            case
            for plan in self.plans.values()
            for case in plan.cases
        )
        semantics_by_id: dict[str, str] = {}
        for case in cases:
            expected_id, canonical = case_identity(case.spec)
            self.assertEqual(case.case_id, expected_id)
            previous = semantics_by_id.setdefault(case.case_id, canonical)
            self.assertEqual(
                previous,
                canonical,
                f"case_id collision for {case.case_id}",
            )

        self.assertEqual(
            len(semantics_by_id),
            len(cases),
            "full profiles must not repeat a case_id",
        )

    def test_all_case_paths_remain_inside_the_owned_run_directory(self) -> None:
        with tempfile.TemporaryDirectory(prefix="ads-final-config-audit-") as temporary:
            repository = Path(temporary) / "repository"
            repository.mkdir()
            store = ResultStore(repository)
            expected_cases_root = (
                repository.resolve() / "benchmarks" / "audit-run" / "cases"
            )

            for plan in self.plans.values():
                for case in plan.cases:
                    case_path = store.case_directory_path(
                        "audit-run", case.case_id
                    )
                    self.assertEqual(case_path.parent, expected_cases_root)
                    self.assertEqual(case_path.name, case.case_id)


if __name__ == "__main__":
    unittest.main(verbosity=2)
