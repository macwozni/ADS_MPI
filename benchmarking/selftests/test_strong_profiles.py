from __future__ import annotations

import json
from pathlib import Path
import unittest

from ads_benchmark.catalog import build_catalog
from ads_benchmark.framework.config import load_profiles
from ads_benchmark.framework.errors import ValidationError
from ads_benchmark.framework.filtering import CaseFilters
from ads_benchmark.framework.planner import Planner


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
        (2, (1, 2, 1)),
        (2, (1, 1, 2)),
    }
)
CLUSTER_LAYOUTS = LOCAL_LAYOUTS | frozenset(
    {
        (4, (2, 2, 1)),
        (4, (2, 1, 2)),
        (4, (1, 2, 2)),
        (8, (2, 2, 2)),
    }
)


class StrongScalingProfileTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.catalog = build_catalog()
        cls.profiles = load_profiles(CONFIG_DIRECTORY)
        cls.planner = Planner(cls.profiles, cls.catalog)

    def assert_common_strong_properties(self, profile_name: str) -> None:
        plan = self.planner.plan(profile_name)
        self.assertEqual(
            {case.spec.problem for case in plan.cases},
            {"igrm_l2", "igrm_heat", "pure_diffusion_igrm"},
        )
        self.assertEqual(
            {case.spec.scheme for case in plan.cases}, {"dg", "pr", "be"}
        )
        self.assertEqual(
            {
                (case.spec.test_degree, case.spec.trial_degree)
                for case in plan.cases
            },
            DEGREE_PAIRS,
        )
        self.assertEqual(
            {
                (
                    case.spec.time.final_time,
                    case.spec.time.time_step,
                    case.spec.time.steps,
                )
                for case in plan.cases
            },
            {("0.1", "0.025", 4)},
        )
        self.assertEqual(
            {case.spec.exact_case for case in plan.cases}, {"spatial-cosine"}
        )
        self.assertEqual(
            {case.spec.openmp_threads for case in plan.cases}, {1, 2, 4, 8}
        )
        self.assertEqual(
            {
                (
                    case.spec.sampling.points_per_axis,
                    case.spec.sampling.write_samples,
                )
                for case in plan.cases
            },
            {(17, True)},
        )
        self.assertEqual(
            {case.spec.build_profile for case in plan.cases}, {"release"}
        )
        self.assertEqual({case.spec.launcher for case in plan.cases}, {"mpi"})

        normalized = plan.cases[0].spec.to_dict()
        self.assertEqual(
            {
                key: normalized["openmp"][key]
                for key in ("dynamic", "proc_bind", "places")
            },
            {"dynamic": False, "proc_bind": "close", "places": "cores"},
        )
        self.assertEqual(
            normalized["measurement"],
            {
                "warmups": 2,
                "samples": 7,
                "minimum_sample_seconds": "0.05",
                "timeout_seconds": (
                    "3600" if profile_name == "local-scaling" else "14400"
                ),
            },
        )

    def test_local_profile_is_the_complete_2160_case_local_matrix(self) -> None:
        self.assert_common_strong_properties("local-scaling")
        plan = self.planner.plan("local-scaling")
        self.assertEqual(len(plan.cases), 2160)
        self.assertEqual({case.spec.mesh for case in plan.cases}, {(16, 16, 16)})
        self.assertEqual(
            {
                (case.spec.mpi.ranks, case.spec.mpi.process_grid)
                for case in plan.cases
            },
            LOCAL_LAYOUTS,
        )

    def test_cluster_profile_is_the_complete_12960_case_full_matrix(self) -> None:
        self.assert_common_strong_properties("cluster-scaling")
        plan = self.planner.plan("cluster-scaling")
        self.assertEqual(len(plan.cases), 12960)
        self.assertEqual(
            {case.spec.mesh for case in plan.cases},
            {(16, 16, 16), (32, 32, 32), (64, 64, 64)},
        )
        self.assertEqual(
            {
                (case.spec.mpi.ranks, case.spec.mpi.process_grid)
                for case in plan.cases
            },
            CLUSTER_LAYOUTS,
        )

    def test_profile_sources_declare_measurement_and_openmp_policy(self) -> None:
        for profile_name, timeout in (
            ("local-scaling", "3600"),
            ("cluster-scaling", "14400"),
        ):
            with self.subTest(profile=profile_name):
                document = json.loads(
                    (CONFIG_DIRECTORY / f"{profile_name}.json").read_text(
                        encoding="utf-8"
                    )
                )
                self.assertEqual(
                    document["openmp"],
                    {
                        "dynamic": False,
                        "proc_bind": "close",
                        "places": "cores",
                    },
                )
                self.assertEqual(
                    document["execution"],
                    {
                        "warmups": 2,
                        "samples": 7,
                        "minimum_sample_seconds": "0.05",
                        "timeout_seconds": timeout,
                    },
                )

    def test_local_verification_slice_selects_exactly_eight_cases(self) -> None:
        selected_degrees = frozenset(
            {
                ((2, 2, 2), (1, 1, 1)),
                ((3, 3, 3), (2, 2, 2)),
            }
        )
        filters = CaseFilters(
            problems=frozenset({"igrm_l2"}),
            schemes=frozenset({"dg", "pr"}),
            degree_pairs=selected_degrees,
            meshes=frozenset({(16, 16, 16)}),
            mpi_grids=frozenset({(1, 1, 1), (2, 1, 1)}),
            openmp_threads=frozenset({1}),
        )
        plan = self.planner.plan("cluster-scaling", filters)

        self.assertEqual(len(plan.cases), 8)
        self.assertEqual(
            {
                (
                    case.spec.scheme,
                    (case.spec.test_degree, case.spec.trial_degree),
                    case.spec.mpi.process_grid,
                    case.spec.openmp_threads,
                )
                for case in plan.cases
            },
            {
                (scheme, degrees, grid, 1)
                for scheme in ("dg", "pr")
                for degrees in selected_degrees
                for grid in ((1, 1, 1), (2, 1, 1))
            },
        )

    def test_smoke_profile_is_the_eight_case_real_verification_slice(self) -> None:
        plan = self.planner.plan("strong-scaling-smoke")
        self.assertEqual(len(plan.cases), 8)
        self.assertEqual({case.spec.problem for case in plan.cases}, {"igrm_l2"})
        self.assertEqual({case.spec.scheme for case in plan.cases}, {"dg", "pr"})
        self.assertEqual({case.spec.mesh for case in plan.cases}, {(4, 4, 4)})
        self.assertEqual(
            {
                (case.spec.test_degree, case.spec.trial_degree)
                for case in plan.cases
            },
            {
                ((2, 2, 2), (1, 1, 1)),
                ((3, 3, 3), (2, 2, 2)),
            },
        )
        self.assertEqual(
            {
                (case.spec.mpi.ranks, case.spec.mpi.process_grid)
                for case in plan.cases
            },
            {(1, (1, 1, 1)), (2, (2, 1, 1))},
        )
        self.assertEqual({case.spec.openmp_threads for case in plan.cases}, {1})
        self.assertTrue(all(case.spec.sampling.write_samples for case in plan.cases))
        self.assertEqual(
            {
                (
                    case.spec.measurement.warmups,
                    case.spec.measurement.samples,
                    case.spec.measurement.minimum_sample_seconds,
                )
                for case in plan.cases
            },
            {(2, 7, "0.001")},
        )

    def test_filters_cannot_remove_the_required_full_field_reference(self) -> None:
        filters = CaseFilters(mpi_grids=frozenset({(2, 1, 1)}))
        with self.assertRaisesRegex(
            ValidationError, "exactly one MPI=1, OMP=1"
        ):
            self.planner.plan("strong-scaling-smoke", filters)



if __name__ == "__main__":
    unittest.main()
