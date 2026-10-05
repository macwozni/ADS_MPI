from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import tempfile
import unittest

from ads_benchmark.catalog import build_catalog
from ads_benchmark.framework.config import load_profiles
from ads_benchmark.framework.errors import ConfigurationError, ValidationError
from ads_benchmark.framework.filtering import CaseFilters
from ads_benchmark.framework.model import RepositoryState, WeakScalingSpec
from ads_benchmark.framework.planner import Planner, frozen_plan_from_manifest
from ads_benchmark.framework.validation import validate_case


BENCHMARKING_ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIRECTORY = BENCHMARKING_ROOT / "configs"

DEGREE_PAIRS = frozenset(
    [
        ((degree + 1,) * 3, (degree,) * 3)
        for degree in range(1, 9)
    ]
    + [
        ((degree + 2,) * 3, (degree,) * 3)
        for degree in range(1, 8)
    ]
)
LOCAL_LAYOUTS = frozenset(
    {
        (1, (1, 1, 1)),
        (2, (2, 1, 1)),
        (4, (2, 2, 1)),
        (8, (2, 2, 2)),
    }
)
CLUSTER_LAYOUTS = LOCAL_LAYOUTS | frozenset(
    {
        (16, (4, 2, 2)),
        (27, (3, 3, 3)),
        (32, (4, 4, 2)),
        (64, (4, 4, 4)),
    }
)


def _field_key(case) -> tuple[object, ...]:
    spec = case.spec
    return (
        spec.problem,
        spec.scheme,
        spec.exact_case,
        spec.time,
        spec.mesh,
        spec.test_degree,
        spec.trial_degree,
    )


class WeakScalingProfileTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.catalog = build_catalog()
        cls.profiles = load_profiles(CONFIG_DIRECTORY)
        cls.planner = Planner(cls.profiles, cls.catalog)

    def assert_common_full_properties(self, profile_name: str) -> None:
        plan = self.planner.plan(profile_name)
        measurements = [
            case
            for case in plan.cases
            if case.spec.weak_scaling.role == "measurement"
        ]
        self.assertEqual(
            {case.spec.problem for case in measurements},
            {"igrm_l2", "igrm_heat", "pure_diffusion_igrm"},
        )
        self.assertEqual(
            {case.spec.scheme for case in measurements}, {"dg", "pr", "be"}
        )
        self.assertEqual(
            {
                (case.spec.test_degree, case.spec.trial_degree)
                for case in measurements
            },
            DEGREE_PAIRS,
        )
        self.assertEqual(
            {case.spec.weak_scaling.local_elements for case in measurements},
            {(8, 8, 8), (12, 12, 12), (16, 16, 16)},
        )
        self.assertEqual(
            {case.spec.weak_scaling.workload_basis for case in measurements},
            {"per-rank"},
        )
        self.assertEqual(
            {case.spec.openmp_threads for case in measurements}, {1, 2, 4, 8}
        )
        self.assertEqual(
            {case.spec.build_profile for case in plan.cases}, {"release"}
        )
        self.assertTrue(
            all(case.spec.sampling.write_samples for case in plan.cases)
        )
        self.assertTrue(
            all(
                case.spec.measurement.warmups == 2
                and case.spec.measurement.samples == 7
                and case.spec.measurement.minimum_sample_seconds == "0.05"
                and case.spec.openmp_dynamic is False
                and case.spec.openmp_proc_bind == "close"
                and case.spec.openmp_places == "cores"
                for case in plan.cases
            )
        )
        for case in measurements:
            local = case.spec.weak_scaling.local_elements
            grid = case.spec.mpi.process_grid
            self.assertEqual(
                case.spec.mesh,
                tuple(a * b for a, b in zip(local, grid, strict=True)),
            )

        modes = {
            (
                "mpi" if case.spec.mpi.ranks > 1 else "serial",
                "openmp" if case.spec.openmp_threads > 1 else "single-thread",
            )
            for case in measurements
        }
        self.assertIn(("mpi", "single-thread"), modes)
        self.assertIn(("serial", "openmp"), modes)
        self.assertIn(("mpi", "openmp"), modes)

        references: dict[tuple[object, ...], int] = {}
        for case in plan.cases:
            weak = case.spec.weak_scaling
            is_serial_reference = (
                case.spec.mpi.ranks == 1
                and case.spec.mpi.process_grid == (1, 1, 1)
                and case.spec.openmp_threads == 1
                and weak.role in {"measurement", "field-reference"}
            )
            if is_serial_reference:
                key = _field_key(case)
                references[key] = references.get(key, 0) + 1
            if weak.role == "field-reference":
                self.assertEqual(weak.local_elements, case.spec.mesh)
        self.assertTrue(references)
        self.assertEqual(set(references.values()), {1})
        self.assertTrue(all(_field_key(case) in references for case in measurements))

    def test_local_profile_is_complete(self) -> None:
        self.assert_common_full_properties("local-weak-scaling")
        plan = self.planner.plan("local-weak-scaling")
        self.assertEqual(len(plan.cases), 7560)
        measurements = [
            case for case in plan.cases if case.spec.weak_scaling.role == "measurement"
        ]
        self.assertEqual(len(measurements), 6480)
        self.assertEqual(
            {
                (case.spec.mpi.ranks, case.spec.mpi.process_grid)
                for case in measurements
            },
            LOCAL_LAYOUTS,
        )

    def test_cluster_profile_adds_balanced_and_uneven_layouts(self) -> None:
        self.assert_common_full_properties("cluster-weak-scaling")
        plan = self.planner.plan("cluster-weak-scaling")
        self.assertEqual(len(plan.cases), 14985)
        measurements = [
            case for case in plan.cases if case.spec.weak_scaling.role == "measurement"
        ]
        self.assertEqual(len(measurements), 12960)
        self.assertEqual(
            {
                (case.spec.mpi.ranks, case.spec.mpi.process_grid)
                for case in measurements
            },
            CLUSTER_LAYOUTS,
        )

    def test_smoke_has_two_measurements_and_one_required_helper(self) -> None:
        plan = self.planner.plan("weak-scaling-smoke")
        self.assertEqual(len(plan.cases), 3)
        self.assertEqual(
            sorted(case.spec.weak_scaling.role for case in plan.cases),
            ["field-reference", "measurement", "measurement"],
        )
        self.assertEqual(
            {
                (case.spec.mesh, case.spec.mpi.process_grid)
                for case in plan.cases
                if case.spec.weak_scaling.role == "measurement"
            },
            {((4, 4, 4), (1, 1, 1)), ((8, 4, 4), (2, 1, 1))},
        )

    def test_filters_select_measurements_then_add_field_helpers(self) -> None:
        plan = self.planner.plan(
            "weak-scaling-smoke",
            CaseFilters(mpi_grids=frozenset({(2, 1, 1)})),
        )
        self.assertEqual(len(plan.cases), 3)
        self.assertEqual(
            [case.spec.weak_scaling.role for case in plan.cases].count(
                "measurement"
            ),
            2,
        )
        self.assertEqual(
            [case.spec.weak_scaling.role for case in plan.cases].count(
                "field-reference"
            ),
            1,
        )
        with self.assertRaisesRegex(ValidationError, "selected no cases"):
            self.planner.plan(
                "weak-scaling-smoke",
                CaseFilters(problems=frozenset({"igrm_heat"})),
            )

    def test_manifest_round_trip_preserves_weak_metadata(self) -> None:
        plan = self.planner.plan("weak-scaling-smoke")
        manifest = plan.manifest(
            run_id="weak-round-trip",
            repository=RepositoryState(commit="a" * 40, dirty=False),
            created_at="2026-10-05T00:00:00+00:00",
        )
        frozen = frozen_plan_from_manifest(manifest, self.catalog)
        self.assertEqual(
            [case.to_dict() for case in frozen.cases],
            [case.to_dict() for case in plan.cases],
        )

    def test_weak_case_metadata_is_strictly_validated(self) -> None:
        case = next(
            case.spec
            for case in self.planner.plan("weak-scaling-smoke").cases
            if case.spec.weak_scaling.role == "measurement"
            and case.spec.mpi.ranks == 2
        )
        with self.assertRaisesRegex(ValidationError, "global mesh"):
            validate_case(replace(case, mesh=(4, 4, 4)), self.catalog)
        with self.assertRaisesRegex(ValidationError, "workload_basis"):
            validate_case(
                replace(
                    case,
                    weak_scaling=replace(
                        case.weak_scaling, workload_basis="per-core"
                    ),
                ),
                self.catalog,
            )
        with self.assertRaisesRegex(ValidationError, "field-reference"):
            validate_case(
                replace(
                    case,
                    weak_scaling=WeakScalingSpec(
                        local_elements=case.mesh,
                        workload_basis="per-rank",
                        role="field-reference",
                    ),
                ),
                self.catalog,
            )

    def test_profile_schema_rejects_implicit_or_per_core_work(self) -> None:
        document = json.loads(
            (CONFIG_DIRECTORY / "weak-scaling-smoke.json").read_text(
                encoding="utf-8"
            )
        )
        with tempfile.TemporaryDirectory(prefix="ads-weak-profile-") as temporary:
            directory = Path(temporary)
            missing = dict(document)
            del missing["weak_scaling"]
            (directory / "missing.json").write_text(
                json.dumps(missing), encoding="utf-8"
            )
            with self.assertRaisesRegex(ConfigurationError, "require weak_scaling"):
                load_profiles(directory)

        document["weak_scaling"]["workload_basis"] = "per-core"
        with tempfile.TemporaryDirectory(prefix="ads-weak-profile-") as temporary:
            path = Path(temporary) / "per-core.json"
            path.write_text(json.dumps(document), encoding="utf-8")
            with self.assertRaisesRegex(ConfigurationError, "must be per-rank"):
                load_profiles(path.parent)


if __name__ == "__main__":
    unittest.main(verbosity=2)
