from __future__ import annotations

from dataclasses import replace
from math import prod
from pathlib import Path
import unittest

from ads_benchmark.catalog import build_catalog
from ads_benchmark.framework.config import load_profiles
from ads_benchmark.framework.errors import ValidationError
from ads_benchmark.framework.model import MpiSpec
from ads_benchmark.framework.planner import Planner
from ads_benchmark.framework.validation import validate_case


BENCHMARKING_ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIRECTORY = BENCHMARKING_ROOT / "configs"
PROBLEMS = {"igrm_l2", "igrm_heat", "pure_diffusion_igrm"}
SCHEMES = {"dg", "pr", "be"}
LAYOUTS = {
    (1, (1, 1, 1)),
    (2, (2, 1, 1)),
    (2, (1, 2, 1)),
    (2, (1, 1, 2)),
    (4, (2, 2, 1)),
    (4, (2, 1, 2)),
    (4, (1, 2, 2)),
    (8, (2, 2, 2)),
    (6, (3, 2, 1)),
}


class ValidationProfileTests(unittest.TestCase):
    def setUp(self) -> None:
        self.catalog = build_catalog()
        profiles = load_profiles(CONFIG_DIRECTORY)
        self.planner = Planner(profiles, self.catalog)

    def assert_common_validation_configuration(self, profile_name: str) -> None:
        plan = self.planner.plan(profile_name)
        self.assertEqual({case.spec.family for case in plan.cases}, {"validation"})
        self.assertEqual(
            {case.spec.exact_case for case in plan.cases}, {"spatial-cosine"}
        )
        self.assertEqual({case.spec.time.final_time for case in plan.cases}, {"0.01"})
        self.assertEqual({case.spec.time.steps for case in plan.cases}, {128, 256})
        self.assertEqual({case.spec.mesh for case in plan.cases}, {(2, 2, 2)})
        self.assertEqual(
            {
                (case.spec.test_degree, case.spec.trial_degree)
                for case in plan.cases
            },
            {((3, 3, 3), (2, 2, 2))},
        )
        self.assertEqual(
            {
                (case.spec.sampling.points_per_axis, case.spec.sampling.write_samples)
                for case in plan.cases
            },
            {(17, True)},
        )
        self.assertEqual({case.spec.build_profile for case in plan.cases}, {"release"})
        self.assertEqual({case.spec.launcher for case in plan.cases}, {"mpi"})
        self.assertEqual(
            {
                (
                    case.spec.measurement.warmups,
                    case.spec.measurement.samples,
                    case.spec.measurement.timeout_seconds,
                )
                for case in plan.cases
            },
            {(0, 1, "1800")},
        )
        self.assertEqual(len({case.case_id for case in plan.cases}), len(plan.cases))

    def assert_parallel_matrix(self, profile_name: str) -> None:
        plan = self.planner.plan(profile_name)
        expected_parallel_points = {
            (ranks, grid, threads)
            for ranks, grid in LAYOUTS
            for threads in (1, 4)
        }
        groups: dict[tuple[object, ...], set[tuple[int, tuple[int, int, int], int]]] = {}

        for case in plan.cases:
            spec = case.spec
            validate_case(spec, self.catalog)
            self.assertEqual(spec.mpi.ranks, prod(spec.mpi.process_grid))
            for processes, elements, trial_degree in zip(
                spec.mpi.process_grid,
                spec.mesh,
                spec.trial_degree,
                strict=True,
            ):
                self.assertLessEqual(processes, elements + trial_degree)

            physical_key = (
                spec.problem,
                spec.scheme,
                spec.exact_case,
                spec.time,
                spec.mesh,
                spec.test_degree,
                spec.trial_degree,
                spec.sampling,
                spec.measurement,
                spec.build_profile,
                spec.launcher,
            )
            groups.setdefault(physical_key, set()).add(
                (spec.mpi.ranks, spec.mpi.process_grid, spec.openmp_threads)
            )

        self.assertTrue(groups)
        for parallel_points in groups.values():
            self.assertEqual(parallel_points, expected_parallel_points)
            self.assertIn((1, (1, 1, 1), 1), parallel_points)

    def test_smoke_is_the_complete_parallel_matrix_for_one_case(self) -> None:
        self.assert_common_validation_configuration("validation-smoke")
        self.assert_parallel_matrix("validation-smoke")
        plan = self.planner.plan("validation-smoke")
        self.assertEqual(len(plan.cases), 36)
        self.assertEqual({case.spec.problem for case in plan.cases}, {"igrm_l2"})
        self.assertEqual({case.spec.scheme for case in plan.cases}, {"dg"})

    def test_full_covers_all_problem_scheme_pairs(self) -> None:
        self.assert_common_validation_configuration("validation-full")
        self.assert_parallel_matrix("validation-full")
        plan = self.planner.plan("validation-full")
        self.assertEqual(len(plan.cases), 324)
        self.assertEqual({case.spec.problem for case in plan.cases}, PROBLEMS)
        self.assertEqual({case.spec.scheme for case in plan.cases}, SCHEMES)
        self.assertEqual(
            {(case.spec.problem, case.spec.scheme) for case in plan.cases},
            {(problem, scheme) for problem in PROBLEMS for scheme in SCHEMES},
        )

    def test_irregular_layout_and_openmp_only_cases_are_explicit(self) -> None:
        plan = self.planner.plan("validation-smoke")
        parallel_points = {
            (
                case.spec.mpi.ranks,
                case.spec.mpi.process_grid,
                case.spec.openmp_threads,
            )
            for case in plan.cases
        }
        self.assertIn((6, (3, 2, 1), 1), parallel_points)
        self.assertIn((6, (3, 2, 1), 4), parallel_points)
        self.assertIn((1, (1, 1, 1), 4), parallel_points)

    def test_existing_validation_rejects_bad_np_and_distribution(self) -> None:
        reference = next(
            case.spec
            for case in self.planner.plan("validation-smoke").cases
            if case.spec.mpi == MpiSpec(1, (1, 1, 1))
            and case.spec.openmp_threads == 1
        )
        with self.assertRaisesRegex(ValidationError, "must equal"):
            validate_case(
                replace(reference, mpi=MpiSpec(2, (1, 1, 1))),
                self.catalog,
            )
        with self.assertRaisesRegex(ValidationError, "trial-space DOFs"):
            validate_case(
                replace(reference, mpi=MpiSpec(5, (5, 1, 1))),
                self.catalog,
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)
