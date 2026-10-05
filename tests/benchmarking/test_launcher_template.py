from __future__ import annotations

from dataclasses import replace
import os
from pathlib import Path
import unittest
from unittest.mock import patch

from ads_benchmark.catalog import build_catalog
from ads_benchmark.components.planning import (
    MpiLauncher,
    TemplateLauncher,
    default_mpi_launcher,
)
from ads_benchmark.framework.config import load_profiles
from ads_benchmark.framework.planner import Planner
from benchmark_paths import BENCHMARKING_ROOT


class LauncherTemplateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        planner = Planner(
            load_profiles(BENCHMARKING_ROOT / "configs"), build_catalog()
        )
        base = next(
            case.spec
            for case in planner.plan("strong-scaling-smoke").cases
            if case.spec.mpi.ranks == 2
        )
        cls.hybrid_case = replace(base, openmp_threads=4)

    def test_scheduler_template_expands_to_shell_free_final_argv(self) -> None:
        launcher = default_mpi_launcher(
            command_template=(
                "srun --ntasks={ranks} --cpus-per-task={threads} "
                "--distribution={procx}x{procy}x{procz} {payload}"
            ),
            available_mpi_slots=8,
            available_cpu_slots=32,
        )
        self.assertIsInstance(launcher, TemplateLauncher)
        launcher.validate_resources(self.hybrid_case)
        self.assertEqual(
            tuple(launcher.command(("solver", "arg with space"), self.hybrid_case)),
            (
                "srun",
                "--ntasks=2",
                "--cpus-per-task=4",
                "--distribution=2x1x1",
                "solver",
                "arg with space",
            ),
        )
        self.assertEqual(launcher.name, "mpi")

    def test_environment_template_and_slurm_allocation_are_honored(self) -> None:
        environment = {
            "BENCHMARK_LAUNCHER_TEMPLATE": (
                "srun -n {ranks} -c {threads} {payload}"
            ),
            "SLURM_NTASKS": "2",
            "SLURM_CPUS_PER_TASK": "4",
        }
        with patch.dict(os.environ, environment, clear=True):
            launcher = default_mpi_launcher()
        self.assertIsInstance(launcher, TemplateLauncher)
        launcher.validate_resources(self.hybrid_case)
        with self.assertRaisesRegex(ValueError, "only 8 CPU slots"):
            launcher.validate_resources(
                replace(self.hybrid_case, openmp_threads=8)
            )

    def test_explicit_capacities_cannot_exceed_detected_slurm_allocation(self) -> None:
        environment = {
            "SLURM_NTASKS": "2",
            "SLURM_CPUS_PER_TASK": "4",
        }
        with patch.dict(os.environ, environment, clear=True):
            launcher = default_mpi_launcher(
                command_template="srun -n {ranks} -c {threads} {payload}",
                available_mpi_slots=8,
                available_cpu_slots=64,
            )
        self.assertIsInstance(launcher, TemplateLauncher)
        self.assertEqual(launcher.available_mpi_slots, 2)
        self.assertEqual(launcher.available_cpu_slots, 8)
        with self.assertRaisesRegex(ValueError, "only 2 MPI slots"):
            launcher.validate_resources(
                replace(
                    self.hybrid_case,
                    mpi=replace(
                        self.hybrid_case.mpi,
                        ranks=4,
                        process_grid=(2, 2, 1),
                    ),
                )
            )

    def test_resource_shortages_are_reported_before_execution(self) -> None:
        launcher = TemplateLauncher(
            arguments=("srun", "-n", "{ranks}", "-c", "{threads}", "{payload}"),
            available_mpi_slots=1,
            available_cpu_slots=4,
            allocated_threads_per_rank=2,
        )
        with self.assertRaisesRegex(ValueError, "only 1 MPI slots"):
            launcher.validate_resources(self.hybrid_case)

        launcher = replace(launcher, available_mpi_slots=2)
        with self.assertRaisesRegex(ValueError, "only 4 CPU slots"):
            launcher.validate_resources(self.hybrid_case)

        launcher = replace(launcher, available_cpu_slots=8)
        with self.assertRaisesRegex(ValueError, "only 2"):
            launcher.validate_resources(self.hybrid_case)

    def test_templates_reject_ambiguous_rank_and_binding_contracts(self) -> None:
        invalid = (
            (
                TemplateLauncher(("srun", "{ranks}", "{payload}", "tail")),
                "final standalone",
            ),
            (
                TemplateLauncher(("srun", "{threads}", "{payload}")),
                "expose.*ranks",
            ),
            (
                TemplateLauncher(("srun", "{ranks}", "{mystery}", "{payload}")),
                "unknown.*mystery",
            ),
        )
        for launcher, message in invalid:
            with self.subTest(arguments=launcher.arguments):
                with self.assertRaisesRegex(ValueError, message):
                    launcher.validate_case(self.hybrid_case)

        no_thread_binding = TemplateLauncher(
            ("srun", "--ntasks={ranks}", "{payload}")
        )
        with self.assertRaisesRegex(ValueError, "requires.*threads"):
            no_thread_binding.validate_case(self.hybrid_case)

    def test_empty_public_make_value_keeps_local_mpiexec(self) -> None:
        with patch.dict(
            os.environ,
            {"BENCHMARK_LAUNCHER_TEMPLATE": "   ", "MPIEXEC": "mpiexec"},
            clear=True,
        ):
            launcher = default_mpi_launcher()
        self.assertIsInstance(launcher, MpiLauncher)

    def test_invalid_slurm_capacity_is_not_silently_ignored(self) -> None:
        with patch.dict(
            os.environ,
            {
                "BENCHMARK_LAUNCHER_TEMPLATE": "srun -n {ranks} {payload}",
                "SLURM_NTASKS": "many",
            },
            clear=True,
        ):
            with self.assertRaisesRegex(ValueError, "SLURM_NTASKS"):
                default_mpi_launcher()


if __name__ == "__main__":
    unittest.main(verbosity=2)
