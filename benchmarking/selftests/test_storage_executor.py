from __future__ import annotations

import json
from pathlib import Path
import tempfile
import time
import unittest

from ads_benchmark.catalog import build_catalog
from ads_benchmark.framework.config import load_profiles
from ads_benchmark.framework.errors import (
    ExecutionError,
    StorageError,
    ValidationError,
)
from ads_benchmark.framework.executor import Executor
from ads_benchmark.framework.model import RepositoryState
from ads_benchmark.framework.planner import Planner
from ads_benchmark.framework.storage import ResultStore
from fake_adapter import FakeAdapter


def owned_manifest(run_id: str) -> dict[str, object]:
    return {
        "schema_version": 1,
        "kind": "ads-benchmark-plan",
        "run_id": run_id,
    }


def fake_profile(timeout: str = "2") -> dict[str, object]:
    return {
        "schema_version": 1,
        "name": "fake-profile",
        "description": "one executable fake case",
        "family": "temporal",
        "exact_cases": ["temporal-polynomial"],
        "problems": ["fake"],
        "schemes": ["dg"],
        "time_discretizations": [{"final_time": "0.1", "steps": 4}],
        "meshes": [[2, 2, 2]],
        "degree_pairs": [{"test": [4, 4, 4], "trial": [3, 3, 3]}],
        "process_layouts": [{"ranks": 1, "grid": [1, 1, 1]}],
        "thread_counts": [2],
        "execution": {"warmups": 0, "samples": 1, "timeout_seconds": timeout},
        "build_profiles": ["debug"],
        "launcher": "direct",
    }


class BrokenEnvironmentLauncher:
    name = "broken-environment"

    def validate_case(self, case) -> None:
        return None

    def command(self, payload, case):
        return payload

    def environment(self, case):
        raise ValueError("intentional launcher environment failure")


class StorageTests(unittest.TestCase):
    def test_results_root_is_fixed_and_symlinks_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="ads-store-path-") as temporary:
            repository = Path(temporary) / "repository"
            repository.mkdir()
            invalid_roots = (
                Path(temporary),
                repository,
                repository / "benchmarks-escape",
                repository / "benchmarks" / "nested",
            )
            for root in invalid_roots:
                with self.subTest(root=root):
                    with self.assertRaisesRegex(StorageError, "exactly"):
                        ResultStore(repository, root)

            outside = Path(temporary) / "outside"
            outside.mkdir()
            (repository / "benchmarks").symlink_to(outside, target_is_directory=True)
            with self.assertRaisesRegex(StorageError, "exactly"):
                ResultStore(repository)

    def test_ids_overwrite_and_legacy_data_are_protected(self) -> None:
        invalid_ids = (
            "",
            ".",
            "..",
            "../escape",
            "/absolute",
            "with/slash",
            "with\\backslash",
            "control\nline",
        )
        with tempfile.TemporaryDirectory(prefix="ads-store-id-") as temporary:
            repository = Path(temporary) / "repository"
            repository.mkdir()
            store = ResultStore(repository)
            for run_id in invalid_ids:
                with self.subTest(run_id=run_id):
                    with self.assertRaisesRegex(StorageError, "unsafe run_id"):
                        store.create_run(run_id, {"kind": "test"})

            legacy = repository / "benchmarks" / "igrm_strong_scaling"
            legacy.mkdir(parents=True)
            sentinel = legacy / "results.csv"
            sentinel.write_text("user-data\n", encoding="utf-8")
            with self.assertRaisesRegex(StorageError, "already exists"):
                store.create_run(
                    "igrm_strong_scaling", owned_manifest("igrm_strong_scaling")
                )
            with self.assertRaisesRegex(StorageError, "ownership manifest"):
                store.create_case_directory("igrm_strong_scaling", "fake-case")
            self.assertEqual(sentinel.read_text(encoding="utf-8"), "user-data\n")

            manifest = owned_manifest("new-run")
            run_path = store.create_run("new-run", manifest)
            self.assertEqual(
                json.loads((run_path / "manifest.json").read_text(encoding="utf-8")),
                manifest,
            )
            self.assertFalse((run_path / "cases").exists())
            with self.assertRaisesRegex(StorageError, "already exists"):
                store.create_run("new-run", manifest)

    def test_case_writes_require_owned_exact_layout_and_reject_swapped_symlinks(self) -> None:
        with tempfile.TemporaryDirectory(prefix="ads-store-case-") as temporary:
            repository = Path(temporary) / "repository"
            repository.mkdir()
            store = ResultStore(repository)
            store.create_run("owned", owned_manifest("owned"))
            case_directory = store.create_case_directory("owned", "case-one")

            with self.assertRaisesRegex(StorageError, "invalid result layout"):
                store.write_status(
                    repository / "benchmarks" / "owned", {"state": "bad"}
                )

            outside = Path(temporary) / "outside"
            outside.mkdir()
            original = case_directory.with_name("case-original")
            case_directory.rename(original)
            case_directory.symlink_to(outside, target_is_directory=True)
            with self.assertRaisesRegex(StorageError, "missing or unsafe"):
                store.write_status(case_directory, {"state": "bad"})
            self.assertFalse((outside / "status.json").exists())


class ExecutorContractTests(unittest.TestCase):
    def _prepare_mode(
        self,
        mode: str,
        timeout: str = "2",
        launcher: BrokenEnvironmentLauncher | None = None,
    ):
        temporary = tempfile.TemporaryDirectory(prefix=f"ads-executor-{mode}-")
        self.addCleanup(temporary.cleanup)
        repository = Path(temporary.name) / "repository"
        repository.mkdir()
        configuration = repository / "profiles"
        configuration.mkdir()
        catalog = build_catalog()
        catalog.adapters.register("fake", FakeAdapter(mode=mode))
        profile = fake_profile(timeout)
        if launcher is not None:
            catalog.launchers.register(launcher.name, launcher)
            profile["launcher"] = launcher.name
        (configuration / "fake.json").write_text(
            json.dumps(profile), encoding="utf-8"
        )
        plan = Planner(load_profiles(configuration), catalog).plan("fake-profile")
        case = plan.cases[0]
        store = ResultStore(repository)
        manifest = plan.manifest(
            run_id="contract",
            repository=RepositoryState(commit="0" * 40, dirty=False),
            created_at="2026-01-01T00:00:00+00:00",
        )
        store.create_run("contract", manifest)
        case_directory = repository / "benchmarks" / "contract" / "cases" / case.case_id
        return Executor(catalog, store, repository), case, case_directory

    def _run_mode(self, mode: str, timeout: str = "2"):
        executor, case, case_directory = self._prepare_mode(mode, timeout)
        status = executor.execute("contract", case)
        return status, case_directory

    def test_registered_fake_adapter_runs_and_parses_without_engine_changes(self) -> None:
        status, case_directory = self._run_mode("success")
        self.assertEqual(status["state"], "passed")
        persisted_status = json.loads(
            (case_directory / "status.json").read_text(encoding="utf-8")
        )
        result = json.loads(
            (case_directory / "result.json").read_text(encoding="utf-8")
        )
        self.assertEqual(persisted_status["state"], "passed")
        self.assertEqual(
            result["domain_result"],
            {
                "checksum": "fake-ok",
                "steps": 4,
                "working_directory": str(case_directory),
            },
        )
        self.assertIn("fake stdout", (case_directory / "stdout.log").read_text())
        self.assertIn("fake stderr", (case_directory / "stderr.log").read_text())
        self.assertEqual(result["configuration"]["openmp"]["threads"], 2)

    def test_nonzero_parser_failure_and_timeout_have_no_fake_result(self) -> None:
        cases = (
            ("nonzero", "2", "failed", "status 7"),
            ("bad", "2", "failed", "parser failed"),
            ("sleep", "0.05", "timeout", "exceeded timeout"),
            ("child", "0.05", "timeout", "exceeded timeout"),
            ("nan", "2", "failed", "serialization failed"),
        )
        for mode, timeout, expected_state, error_text in cases:
            with self.subTest(mode=mode):
                started = time.monotonic()
                status, case_directory = self._run_mode(mode, timeout)
                self.assertLess(time.monotonic() - started, 3)
                self.assertEqual(status["state"], expected_state)
                self.assertIn(error_text, status["error"])
                self.assertFalse((case_directory / "result.json").exists())
                self.assertTrue((case_directory / "stdout.log").is_file())
                self.assertTrue((case_directory / "stderr.log").is_file())

    def test_extension_command_and_environment_failures_are_recorded(self) -> None:
        cases = (
            ("nul-argv", None, "NUL-free"),
            (
                "success",
                BrokenEnvironmentLauncher(),
                "intentional launcher environment failure",
            ),
        )
        for mode, launcher, message in cases:
            with self.subTest(mode=mode):
                executor, case, case_directory = self._prepare_mode(
                    mode, launcher=launcher
                )
                with self.assertRaisesRegex(ExecutionError, message):
                    executor.execute("contract", case)
                status = json.loads(
                    (case_directory / "status.json").read_text(encoding="utf-8")
                )
                self.assertEqual(status["state"], "failed")
                self.assertEqual((case_directory / "stdout.log").read_text(), "")
                self.assertEqual((case_directory / "stderr.log").read_text(), "")

    def test_adapter_validation_errors_are_normalized(self) -> None:
        with self.assertRaisesRegex(ValidationError, "intentional fake validation"):
            self._prepare_mode("reject")

    def test_real_stage_one_adapter_is_explicitly_planning_only(self) -> None:
        benchmarking_root = Path(__file__).resolve().parents[1]
        catalog = build_catalog()
        plan = Planner(
            load_profiles(benchmarking_root / "configs"), catalog
        ).plan("smoke")
        with tempfile.TemporaryDirectory(prefix="ads-planning-only-") as temporary:
            repository = Path(temporary) / "repository"
            repository.mkdir()
            store = ResultStore(repository)
            store.create_run("planning-only", owned_manifest("planning-only"))
            with self.assertRaisesRegex(ExecutionError, "planning-only"):
                Executor(catalog, store, repository).execute(
                    "planning-only", plan.cases[0]
                )
            self.assertFalse(
                (repository / "benchmarks" / "planning-only" / "cases").exists()
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)
