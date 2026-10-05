from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from ads_benchmark import cli
from ads_benchmark.framework.errors import ValidationError
from ads_benchmark.framework.model import RepositoryState


def _mpi_profile() -> dict[str, object]:
    return {
        "schema_version": 1,
        "name": "resource-capacity",
        "description": "one hybrid case for declared resource checks",
        "family": "temporal",
        "exact_cases": ["temporal-polynomial"],
        "problems": ["igrm_l2"],
        "schemes": ["dg"],
        "time_discretizations": [{"final_time": "0.1", "steps": 4}],
        "meshes": [[2, 2, 2]],
        "degree_pairs": [{"test": [4, 4, 4], "trial": [3, 3, 3]}],
        "process_layouts": [{"ranks": 4, "grid": [2, 2, 1]}],
        "thread_counts": [2],
        "sampling": {"points_per_axis": 5, "write_samples": False},
        "execution": {"warmups": 0, "samples": 1, "timeout_seconds": "2"},
        "build_profiles": ["debug"],
        "launcher": "mpi",
    }


class ResourceCapacityTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory(prefix="ads-resource-capacity-")
        self.addCleanup(temporary.cleanup)
        self.repository = Path(temporary.name) / "repository"
        self.repository.mkdir()
        self.configuration = self.repository / "profiles"
        self.configuration.mkdir()
        (self.configuration / "resource-capacity.json").write_text(
            json.dumps(_mpi_profile()), encoding="utf-8"
        )

    def _plan(self, *resource_arguments: str) -> tuple[int, str, str]:
        stdout = io.StringIO()
        stderr = io.StringIO()
        arguments = [
            "plan",
            "--repository-root",
            str(self.repository),
            "--config-dir",
            str(self.configuration),
            "--profile",
            "resource-capacity",
            *resource_arguments,
        ]
        with (
            patch(
                "ads_benchmark.cli.inspect_repository",
                return_value=RepositoryState(commit="0" * 40, dirty=False),
            ),
            redirect_stdout(stdout),
            redirect_stderr(stderr),
        ):
            return_code = cli.main(arguments)
        return return_code, stdout.getvalue(), stderr.getvalue()

    def test_omitted_capacity_keeps_structural_planning_behavior(self) -> None:
        return_code, stdout, stderr = self._plan()
        self.assertEqual(return_code, 0, stderr)
        self.assertIn("cases:        1", stdout)
        self.assertEqual(stderr, "")

    def test_exact_declared_capacities_accept_hybrid_case(self) -> None:
        return_code, _, stderr = self._plan(
            "--available-mpi-slots",
            "4",
            "--available-cpu-slots",
            "8",
        )
        self.assertEqual(return_code, 0, stderr)

    def test_mpi_capacity_rejects_rank_count_with_explicit_error(self) -> None:
        return_code, _, stderr = self._plan("--available-mpi-slots", "3")
        self.assertEqual(return_code, 2)
        self.assertIn("case requests 4 MPI ranks", stderr)
        self.assertIn("only 3 MPI slots are available", stderr)

    def test_cpu_capacity_rejects_rank_thread_product_with_explicit_error(self) -> None:
        return_code, _, stderr = self._plan("--available-cpu-slots", "7")
        self.assertEqual(return_code, 2)
        self.assertIn("4 MPI ranks * 2 OpenMP threads = 8 CPU slots", stderr)
        self.assertIn("only 7 CPU slots are available", stderr)

    def test_capacity_is_checked_after_case_filters(self) -> None:
        profile = _mpi_profile()
        profile["process_layouts"] = [
            {"ranks": 1, "grid": [1, 1, 1]},
            {"ranks": 4, "grid": [2, 2, 1]},
        ]
        (self.configuration / "resource-capacity.json").write_text(
            json.dumps(profile), encoding="utf-8"
        )
        return_code, stdout, stderr = self._plan(
            "--mpi-grid",
            "1x1x1",
            "--available-mpi-slots",
            "1",
            "--available-cpu-slots",
            "2",
        )
        self.assertEqual(return_code, 0, stderr)
        self.assertIn("cases:        1", stdout)

    def test_capacity_arguments_must_be_positive(self) -> None:
        stderr = io.StringIO()
        with redirect_stderr(stderr), self.assertRaises(SystemExit) as raised:
            cli.main(["plan", "--available-mpi-slots", "0"])
        self.assertEqual(raised.exception.code, 2)
        self.assertIn("--available-mpi-slots", stderr.getvalue())
        self.assertIn("must be positive", stderr.getvalue())

    def test_real_validation_requires_both_capacity_declarations(self) -> None:
        plan = SimpleNamespace(
            cases=(
                SimpleNamespace(spec=SimpleNamespace(family="validation")),
            )
        )
        options = SimpleNamespace(
            available_mpi_slots=6,
            available_cpu_slots=None,
        )
        with self.assertRaisesRegex(
            ValidationError,
            "validation execution requires.*--available-cpu-slots",
        ):
            cli._require_validation_execution_resources(options, plan)

        options.available_cpu_slots = 24
        cli._require_validation_execution_resources(options, plan)

    def test_real_strong_scaling_requires_both_capacity_declarations(self) -> None:
        plan = SimpleNamespace(
            cases=(SimpleNamespace(spec=SimpleNamespace(family="strong")),)
        )
        options = SimpleNamespace(
            available_mpi_slots=None,
            available_cpu_slots=8,
        )
        with self.assertRaisesRegex(
            ValidationError,
            "strong-scaling execution requires.*--available-mpi-slots",
        ):
            cli._require_validation_execution_resources(options, plan)

        options.available_mpi_slots = 2
        cli._require_validation_execution_resources(options, plan)

    def test_real_weak_scaling_requires_both_capacity_declarations(self) -> None:
        plan = SimpleNamespace(
            cases=(SimpleNamespace(spec=SimpleNamespace(family="weak")),)
        )
        options = SimpleNamespace(
            available_mpi_slots=2,
            available_cpu_slots=None,
        )
        with self.assertRaisesRegex(
            ValidationError,
            "weak-scaling execution requires.*--available-cpu-slots",
        ):
            cli._require_validation_execution_resources(options, plan)

        options.available_cpu_slots = 2
        cli._require_validation_execution_resources(options, plan)


if __name__ == "__main__":
    unittest.main(verbosity=2)
