from __future__ import annotations

from itertools import permutations
from pathlib import Path
import unittest

from ads_benchmark.catalog import build_catalog
from ads_benchmark.framework.config import load_profiles
from ads_benchmark.framework.planner import Planner
from benchmark_paths import CONFIG_DIRECTORY


PROBLEMS = {"igrm_l2", "igrm_heat", "pure_diffusion_igrm"}
SCHEMES = {"dg", "pr", "be"}


class SpatialProfileTests(unittest.TestCase):
    def setUp(self) -> None:
        profiles = load_profiles(CONFIG_DIRECTORY)
        self.planner = Planner(profiles, build_catalog())

    def assert_dt_pair_per_spatial_point(self, profile_name: str) -> None:
        plan = self.planner.plan(profile_name)
        groups: dict[tuple[object, ...], set[int]] = {}
        for case in plan.cases:
            spec = case.spec
            self.assertEqual(spec.time.final_time, "0.01")
            key = (
                spec.problem,
                spec.scheme,
                spec.exact_case,
                spec.mesh,
                spec.test_degree,
                spec.trial_degree,
                spec.mpi,
                spec.openmp_threads,
                spec.sampling,
                spec.measurement,
                spec.build_profile,
                spec.launcher,
            )
            groups.setdefault(key, set()).add(spec.time.steps)
        self.assertTrue(groups)
        self.assertEqual(
            {frozenset(steps) for steps in groups.values()},
            {frozenset({128, 256})},
        )

    def assert_common_spatial_scope(self, profile_name: str) -> None:
        plan = self.planner.plan(profile_name)
        self.assertEqual({case.spec.problem for case in plan.cases}, PROBLEMS)
        self.assertEqual({case.spec.scheme for case in plan.cases}, SCHEMES)
        self.assertEqual({case.spec.exact_case for case in plan.cases}, {"spatial-cosine"})
        self.assertEqual({case.spec.family for case in plan.cases}, {profile_name[0]})
        self.assertEqual({case.spec.build_profile for case in plan.cases}, {"release"})
        self.assertEqual({case.spec.mpi.ranks for case in plan.cases}, {1})
        self.assertEqual({case.spec.openmp_threads for case in plan.cases}, {1})
        self.assertEqual(
            {case.spec.sampling.write_samples for case in plan.cases}, {True}
        )
        expected_sample_points = 65 if profile_name == "h-convergence-full" else 33
        self.assertEqual(
            {case.spec.sampling.points_per_axis for case in plan.cases},
            {expected_sample_points},
        )
        self.assertEqual(len({case.case_id for case in plan.cases}), len(plan.cases))
        self.assert_dt_pair_per_spatial_point(profile_name)

    def test_h_profiles_refine_only_the_mesh_for_each_series(self) -> None:
        expected = {
            "h-convergence-smoke": ({2, 4, 8}, 54),
            "h-convergence-full": ({2, 4, 8, 16, 32}, 90),
        }
        for profile_name, (levels, count) in expected.items():
            with self.subTest(profile=profile_name):
                self.assert_common_spatial_scope(profile_name)
                plan = self.planner.plan(profile_name)
                self.assertEqual(len(plan.cases), count)
                self.assertEqual(
                    {case.spec.mesh for case in plan.cases},
                    {(level, level, level) for level in levels},
                )
                self.assertEqual(
                    {
                        (case.spec.test_degree, case.spec.trial_degree)
                        for case in plan.cases
                    },
                    {((3, 3, 3), (2, 2, 2))},
                )

    def test_full_isotropic_p_profile_has_every_admissible_pair(self) -> None:
        self.assert_common_spatial_scope("p-convergence-full")
        plan = self.planner.plan("p-convergence-full")
        self.assertEqual(len(plan.cases), 270)
        actual = {
            (case.spec.test_degree[0], case.spec.trial_degree[0])
            for case in plan.cases
        }
        expected = {(trial + 1, trial) for trial in range(1, 9)} | {
            (trial + 2, trial) for trial in range(1, 8)
        }
        self.assertEqual(actual, expected)
        for case in plan.cases:
            self.assertEqual(case.spec.mesh, (2, 2, 2))
            self.assertEqual(len(set(case.spec.test_degree)), 1)
            self.assertEqual(len(set(case.spec.trial_degree)), 1)
            self.assertLessEqual(max(case.spec.test_degree), 9)

    def test_isotropic_p_smoke_has_three_levels_for_both_enrichments(self) -> None:
        self.assert_common_spatial_scope("p-convergence-smoke")
        plan = self.planner.plan("p-convergence-smoke")
        self.assertEqual(len(plan.cases), 108)
        self.assertEqual({case.spec.mesh for case in plan.cases}, {(2, 2, 2)})
        actual = {
            (case.spec.test_degree[0], case.spec.trial_degree[0])
            for case in plan.cases
        }
        self.assertEqual(
            actual,
            {(2, 1), (3, 2), (4, 3), (3, 1), (4, 2), (5, 3)},
        )

    def test_anisotropic_full_preserves_all_vectors_and_rotations(self) -> None:
        self.assert_common_spatial_scope("p-anisotropic-full")
        plan = self.planner.plan("p-anisotropic-full")
        self.assertEqual(len(plan.cases), 216)
        self.assertEqual({case.spec.mesh for case in plan.cases}, {(2, 2, 2)})
        pairs = {
            (case.spec.test_degree, case.spec.trial_degree) for case in plan.cases
        }
        rotations = set(permutations((3, 4, 5)))
        self.assertEqual({trial for _, trial in pairs}, rotations)
        self.assertEqual(len(pairs), 12)
        for test, trial in pairs:
            differences = tuple(
                test_value - trial_value
                for test_value, trial_value in zip(test, trial, strict=True)
            )
            self.assertIn(differences, {(1, 1, 1), (2, 2, 2)})
            self.assertLessEqual(max(test), 9)

    def test_anisotropic_smoke_is_a_vector_preserving_rotation_subset(self) -> None:
        self.assert_common_spatial_scope("p-anisotropic-smoke")
        smoke = self.planner.plan("p-anisotropic-smoke")
        full = self.planner.plan("p-anisotropic-full")
        self.assertEqual(len(smoke.cases), 108)
        self.assertEqual({case.spec.mesh for case in smoke.cases}, {(2, 2, 2)})
        smoke_pairs = {
            (case.spec.test_degree, case.spec.trial_degree) for case in smoke.cases
        }
        full_pairs = {
            (case.spec.test_degree, case.spec.trial_degree) for case in full.cases
        }
        self.assertEqual(len(smoke_pairs), 6)
        self.assertLess(smoke_pairs, full_pairs)
        self.assertEqual(
            {trial for _, trial in smoke_pairs},
            {(3, 4, 5), (4, 5, 3), (5, 3, 4)},
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
