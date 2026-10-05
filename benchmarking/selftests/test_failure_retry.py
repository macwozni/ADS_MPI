from __future__ import annotations

import errno
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from ads_benchmark.catalog import build_catalog
from ads_benchmark.framework.config import load_profiles
from ads_benchmark.framework.errors import ExecutionError
from ads_benchmark.framework.executor import (
    FAILURE_KINDS,
    RETRYABLE_FAILURE_KINDS,
    Executor,
)
from ads_benchmark.framework.model import ExecutionContext, RepositoryState
from ads_benchmark.framework.planner import Planner
from ads_benchmark.framework.storage import ResultStore
from fake_adapter import FakeAdapter


def fake_profile(
    *,
    timeout: str = "2",
    launcher: str = "direct",
    minimum_sample_seconds: str | None = None,
) -> dict[str, object]:
    execution: dict[str, object] = {
        "warmups": 0,
        "samples": 1,
        "timeout_seconds": timeout,
    }
    if minimum_sample_seconds is not None:
        execution["minimum_sample_seconds"] = minimum_sample_seconds
    return {
        "schema_version": 1,
        "name": "failure-profile",
        "description": "failure classification and retry contract",
        "family": "temporal",
        "exact_cases": ["temporal-polynomial"],
        "problems": ["fake"],
        "schemes": ["dg"],
        "time_discretizations": [{"final_time": "0.1", "steps": 4}],
        "meshes": [[2, 2, 2]],
        "degree_pairs": [{"test": [4, 4, 4], "trial": [3, 3, 3]}],
        "process_layouts": [{"ranks": 1, "grid": [1, 1, 1]}],
        "thread_counts": [1],
        "sampling": {"points_per_axis": 5, "write_samples": False},
        "execution": execution,
        "build_profiles": ["debug"],
        "launcher": launcher,
    }


class PassthroughMpiLauncher:
    """Exercise MPI failure policy without requiring an MPI installation."""

    name = "mpi"

    def validate_case(self, case) -> None:
        return None

    def validate_resources(self, case) -> None:
        return None

    def command(self, payload, case):
        return tuple(payload)

    def environment(self, case):
        return {}


class MissingMeasurementAdapter:
    name = "fake"
    execution_ready = True

    def __init__(self) -> None:
        self.delegate = FakeAdapter(mode="success")

    def validate_case(self, case) -> None:
        self.delegate.validate_case(case)

    def build_payload_command(self, case, context: ExecutionContext):
        return self.delegate.build_payload_command(case, context)

    def parse_result(self, stdout: str, stderr: str):
        result = dict(self.delegate.parse_result(stdout, stderr))
        result.pop("physical_step_wall_seconds")
        return result

    def validate_result(self, case, result) -> None:
        self.delegate.validate_result(case, result)


class FailureClassificationAndRetryTests(unittest.TestCase):
    def _prepared_executor(
        self,
        *,
        mode: str = "success",
        timeout: str = "2",
        launcher=None,
        adapter=None,
        minimum_sample_seconds: str | None = None,
    ):
        temporary = tempfile.TemporaryDirectory(prefix="ads-failure-retry-")
        self.addCleanup(temporary.cleanup)
        repository = Path(temporary.name) / "repository"
        repository.mkdir()
        configuration = repository / "profiles"
        configuration.mkdir()

        launcher_key = "direct"
        catalog = build_catalog()
        catalog.adapters.register("fake", adapter or FakeAdapter(mode=mode))
        if launcher is not None:
            launcher_key = "test-mpi"
            catalog.launchers.register(launcher_key, launcher)
        (configuration / "failure.json").write_text(
            json.dumps(
                fake_profile(
                    timeout=timeout,
                    launcher=launcher_key,
                    minimum_sample_seconds=minimum_sample_seconds,
                )
            ),
            encoding="utf-8",
        )
        plan = Planner(load_profiles(configuration), catalog).plan(
            "failure-profile"
        )
        store = ResultStore(repository)
        store.create_run(
            "failure-run",
            plan.manifest(
                run_id="failure-run",
                repository=RepositoryState(commit="0" * 40, dirty=False),
                created_at="2026-01-01T00:00:00+00:00",
            ),
        )
        case = plan.cases[0]
        case_directory = (
            repository / "benchmarks" / "failure-run" / "cases" / case.case_id
        )
        return Executor(catalog, store, repository), case, case_directory

    @staticmethod
    def _persisted_status(case_directory: Path) -> dict[str, object]:
        return json.loads(
            (case_directory / "status.json").read_text(encoding="utf-8")
        )

    def test_failure_and_retry_kind_sets_are_exact(self) -> None:
        self.assertEqual(
            FAILURE_KINDS,
            {"numerical", "mpi", "timeout", "resource", "configuration"},
        )
        self.assertEqual(
            RETRYABLE_FAILURE_KINDS,
            {"mpi", "timeout", "resource"},
        )

    def test_configuration_and_numerical_failures_are_not_retried(self) -> None:
        for mode, expected_kind in (
            ("nul-argv", "configuration"),
            ("bad", "numerical"),
            ("reject-result", "numerical"),
        ):
            with self.subTest(mode=mode):
                executor, case, case_directory = self._prepared_executor(mode=mode)
                summary = executor.execute_frozen(
                    "failure-run",
                    (case,),
                    preflight=False,
                    max_retries=3,
                )
                status = self._persisted_status(case_directory)
                self.assertEqual(status["state"], "failed")
                self.assertEqual(status["failure_kind"], expected_kind)
                self.assertEqual(len(status["attempts"]), 1)
                self.assertEqual((summary.passed, summary.failed), (0, 1))
                self.assertFalse((case_directory / "result.json").exists())

    def test_measurement_failure_is_numerical_and_not_retried(self) -> None:
        executor, case, case_directory = self._prepared_executor(
            adapter=MissingMeasurementAdapter(),
            minimum_sample_seconds="0.1",
        )
        summary = executor.execute_frozen(
            "failure-run", (case,), preflight=False, max_retries=2
        )
        status = self._persisted_status(case_directory)
        self.assertEqual(status["state"], "failed")
        self.assertEqual(status["failure_kind"], "numerical")
        self.assertIn("result measurement failed", status["error"])
        self.assertEqual(len(status["attempts"]), 1)
        self.assertEqual((summary.passed, summary.failed), (0, 1))

    def test_timeout_and_mpi_exit_retry_but_never_become_passed(self) -> None:
        cases = (
            ("sleep", "0.05", None, "timeout", "timeout"),
            ("nonzero", "2", PassthroughMpiLauncher(), "failed", "mpi"),
        )
        for mode, timeout, launcher, expected_state, expected_kind in cases:
            with self.subTest(mode=mode):
                executor, case, case_directory = self._prepared_executor(
                    mode=mode, timeout=timeout, launcher=launcher
                )
                summary = executor.execute_frozen(
                    "failure-run",
                    (case,),
                    preflight=False,
                    max_retries=1,
                )
                status = self._persisted_status(case_directory)
                self.assertEqual(status["state"], expected_state)
                self.assertEqual(status["failure_kind"], expected_kind)
                self.assertEqual(len(status["attempts"]), 2)
                self.assertEqual(
                    {attempt["failure_kind"] for attempt in status["attempts"]},
                    {expected_kind},
                )
                self.assertEqual((summary.passed, summary.failed), (0, 1))
                self.assertFalse((case_directory / "result.json").exists())

    def test_resource_errnos_are_classified_as_resource_failures(self) -> None:
        for error_number in (
            errno.ENOMEM,
            errno.ENOSPC,
            errno.EMFILE,
            errno.ENFILE,
            errno.EAGAIN,
        ):
            with self.subTest(errno=error_number):
                executor, case, case_directory = self._prepared_executor()
                failure = OSError(error_number, os.strerror(error_number))
                with patch(
                    "ads_benchmark.framework.executor.subprocess.Popen",
                    side_effect=failure,
                ):
                    summary = executor.execute_frozen(
                        "failure-run", (case,), preflight=False
                    )
                status = self._persisted_status(case_directory)
                self.assertEqual(status["state"], "failed")
                self.assertEqual(status["failure_kind"], "resource")
                self.assertEqual((summary.passed, summary.failed), (0, 1))
                self.assertFalse((case_directory / "result.json").exists())

    def test_transient_resource_failure_retries_and_then_passes(self) -> None:
        executor, case, case_directory = self._prepared_executor()
        real_popen = subprocess.Popen
        invocation_count = 0

        def transient_popen(*args, **kwargs):
            nonlocal invocation_count
            invocation_count += 1
            if invocation_count == 1:
                raise OSError(errno.EAGAIN, os.strerror(errno.EAGAIN))
            return real_popen(*args, **kwargs)

        with patch(
            "ads_benchmark.framework.executor.subprocess.Popen",
            side_effect=transient_popen,
        ):
            summary = executor.execute_frozen(
                "failure-run", (case,), preflight=False, max_retries=1
            )

        status = self._persisted_status(case_directory)
        self.assertEqual(status["state"], "passed")
        self.assertNotIn("failure_kind", status)
        self.assertEqual(
            [attempt["state"] for attempt in status["attempts"]],
            ["failed", "passed"],
        )
        self.assertEqual(status["attempts"][0]["failure_kind"], "resource")
        self.assertIsNone(status["attempts"][1]["failure_kind"])
        self.assertEqual((summary.passed, summary.failed), (1, 0))
        self.assertTrue((case_directory / "result.json").is_file())
        self.assertIsNotNone(executor.verified_result(case_directory, case))

    def test_invalid_retry_limit_is_configuration_error(self) -> None:
        for invalid in (-1, True, 1.5):
            with self.subTest(max_retries=invalid):
                executor, case, _ = self._prepared_executor()
                with self.assertRaisesRegex(
                    ExecutionError, "max_retries must be a nonnegative integer"
                ):
                    executor.execute_frozen(
                        "failure-run", (case,), max_retries=invalid
                    )


if __name__ == "__main__":
    unittest.main(verbosity=2)
