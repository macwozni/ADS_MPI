from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from copy import deepcopy

from ads_benchmark.catalog import build_catalog
from ads_benchmark.framework.config import load_profiles
from ads_benchmark.framework.errors import (
    ExecutionError,
    StorageError,
    ValidationError,
)
from ads_benchmark.framework.executor import Executor
from ads_benchmark.framework.model import RepositoryState
from ads_benchmark.framework.planner import (
    Planner,
    frozen_plan_from_manifest,
    validate_resume_request,
)
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
        "sampling": {"points_per_axis": 5, "write_samples": False},
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


class PrefixLauncher:
    name = "prefix-launcher"

    def __init__(self, prefix: tuple[str, ...]) -> None:
        self.prefix = prefix

    def validate_case(self, case) -> None:
        return None

    def command(self, payload, case):
        return (*self.prefix, *payload)

    def environment(self, case):
        return {}


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

    def test_execution_lock_refuses_concurrent_writers(self) -> None:
        with tempfile.TemporaryDirectory(prefix="ads-store-lock-") as temporary:
            repository = Path(temporary) / "repository"
            repository.mkdir()
            store = ResultStore(repository)
            store.create_run("owned", owned_manifest("owned"))
            with store.execution_lock("owned"):
                with self.assertRaisesRegex(StorageError, "already being executed"):
                    with store.execution_lock("owned"):
                        self.fail("a second execution lock was acquired")

    def test_open_case_directory_anchors_child_cwd_across_symlink_swap(self) -> None:
        with tempfile.TemporaryDirectory(prefix="ads-store-cwd-") as temporary:
            repository = Path(temporary) / "repository"
            repository.mkdir()
            store = ResultStore(repository)
            store.create_run("owned", owned_manifest("owned"))
            case_directory = store.create_case_directory("owned", "case-one")
            outside = Path(temporary) / "outside"
            outside.mkdir()
            original = case_directory.with_name("case-original")
            with store.open_case_directory(case_directory) as case_fd:
                case_directory.rename(original)
                case_directory.symlink_to(outside, target_is_directory=True)
                subprocess.run(
                    [
                        sys.executable,
                        "-c",
                        "from pathlib import Path; Path('marker').write_text('safe')",
                    ],
                    cwd=f"/proc/self/fd/{case_fd}",
                    pass_fds=(case_fd,),
                    check=True,
                )
            self.assertEqual((original / "marker").read_text(), "safe")
            self.assertFalse((outside / "marker").exists())

    def test_generated_artifact_reader_is_allowlisted_bounded_and_nofollow(self) -> None:
        with tempfile.TemporaryDirectory(prefix="ads-store-artifact-") as temporary:
            repository = Path(temporary) / "repository"
            repository.mkdir()
            store = ResultStore(repository)
            store.create_run("owned", owned_manifest("owned"))
            case_directory = store.create_case_directory("owned", "case-one")
            artifact = case_directory / "field_samples.csv"
            artifact.write_text("x,y,z,numerical,exact,error\n", encoding="utf-8")

            self.assertEqual(
                store.read_generated_artifact(case_directory, "field_samples.csv"),
                "x,y,z,numerical,exact,error\n",
            )
            with self.assertRaisesRegex(StorageError, "unsupported generated"):
                store.read_generated_artifact(case_directory, "stdout.log")

            outside = Path(temporary) / "outside.csv"
            outside.write_text("sentinel\n", encoding="utf-8")
            artifact.unlink()
            artifact.symlink_to(outside)
            with self.assertRaisesRegex(StorageError, "cannot read"):
                store.read_generated_artifact(case_directory, "field_samples.csv")

            artifact.unlink()
            with artifact.open("wb") as stream:
                stream.truncate(64 * 1024 * 1024 + 1)
            with self.assertRaisesRegex(StorageError, "regular file below 64 MiB"):
                store.read_generated_artifact(case_directory, "field_samples.csv")

            artifact.unlink()
            artifact.mkdir()
            with self.assertRaisesRegex(StorageError, "regular file below 64 MiB"):
                store.read_generated_artifact(case_directory, "field_samples.csv")

            artifact.rmdir()
            artifact.write_bytes(b"\xff")
            with self.assertRaisesRegex(StorageError, "cannot read"):
                store.read_generated_artifact(case_directory, "field_samples.csv")


class ExecutorContractTests(unittest.TestCase):
    def _prepare_mode(
        self,
        mode: str,
        timeout: str = "2",
        launcher: BrokenEnvironmentLauncher | None = None,
        *,
        warmups: int = 0,
        samples: int = 1,
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
        profile["execution"] = {
            "warmups": warmups,
            "samples": samples,
            "timeout_seconds": timeout,
        }
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
            ("detached", "0.05", "timeout", "exceeded timeout"),
            ("nan", "2", "failed", "serialization failed"),
            (
                "reject-result",
                "2",
                "failed",
                "result validation failed: intentional fake result rejection",
            ),
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

    def test_real_stage_two_adapters_are_execution_ready(self) -> None:
        catalog = build_catalog()
        for name in ("igrm_l2", "igrm_heat", "pure_diffusion_igrm"):
            with self.subTest(name=name):
                adapter = catalog.adapters.get(name)
                self.assertEqual(adapter.name, name)
                self.assertTrue(adapter.execution_ready)

    def test_direct_execute_rejects_unimplemented_repetitions(self) -> None:
        executor, case, case_directory = self._prepare_mode(
            "success", warmups=2, samples=7
        )
        with self.assertRaisesRegex(
            ExecutionError, "supports only warmups=0 and samples=1"
        ):
            executor.execute("contract", case)
        status = json.loads(
            (case_directory / "status.json").read_text(encoding="utf-8")
        )
        self.assertEqual(status["state"], "failed")
        self.assertIn("requests warmups=2 and samples=7", status["error"])


class FrozenManifestAndResumeTests(unittest.TestCase):
    def _prepared_run(
        self,
        *,
        schemes: list[str] | None = None,
        timeout: str = "2",
        write_samples: bool = False,
        launcher: PrefixLauncher | None = None,
    ):
        temporary = tempfile.TemporaryDirectory(prefix="ads-resume-")
        self.addCleanup(temporary.cleanup)
        repository = Path(temporary.name) / "repository"
        repository.mkdir()
        configuration = repository / "profiles"
        configuration.mkdir()
        profile = fake_profile(timeout)
        profile["schemes"] = schemes or ["dg"]
        profile["sampling"]["write_samples"] = write_samples
        if launcher is not None:
            profile["launcher"] = launcher.name
        (configuration / "fake.json").write_text(
            json.dumps(profile), encoding="utf-8"
        )
        catalog = build_catalog()
        catalog.adapters.register("fake", FakeAdapter(mode="success"))
        if launcher is not None:
            catalog.launchers.register(launcher.name, launcher)
        plan = Planner(load_profiles(configuration), catalog).plan("fake-profile")
        repository_state = RepositoryState(
            commit="0" * 40,
            dirty=False,
            worktree_fingerprint="sha256:" + "a" * 64,
        )
        manifest = plan.manifest(
            run_id="resume-contract",
            repository=repository_state,
            created_at="2026-01-01T00:00:00+00:00",
        )
        store = ResultStore(repository)
        store.create_run("resume-contract", manifest)
        frozen = frozen_plan_from_manifest(store.read_manifest("resume-contract"), catalog)
        return repository, catalog, plan, repository_state, store, manifest, frozen

    def test_frozen_manifest_is_the_validated_execution_source(self) -> None:
        _, catalog, plan, repository_state, _, manifest, frozen = self._prepared_run()
        self.assertEqual(frozen.cases, plan.cases)
        validate_resume_request(frozen, plan, repository_state)

        legacy = deepcopy(manifest)
        legacy["filters"].pop("steps")
        legacy["repository"].pop("worktree_fingerprint")
        legacy_frozen = frozen_plan_from_manifest(legacy, catalog)
        validate_resume_request(
            legacy_frozen,
            plan,
            RepositoryState(commit=repository_state.commit, dirty=False),
        )

        bad_schema = deepcopy(manifest)
        bad_schema["schema_version"] = True
        with self.assertRaisesRegex(ValidationError, "schema_version"):
            frozen_plan_from_manifest(bad_schema, catalog)

        bad_case = deepcopy(manifest)
        bad_case["cases"][0]["configuration"]["openmp"]["threads"] = 3
        with self.assertRaisesRegex(ValidationError, "case_id"):
            frozen_plan_from_manifest(bad_case, catalog)

    def test_resume_rejects_changed_sha_dirty_state_and_configuration(self) -> None:
        _, _, plan, repository_state, _, _, frozen = self._prepared_run()
        with self.assertRaisesRegex(ValidationError, "SHA mismatch"):
            validate_resume_request(
                frozen,
                plan,
                RepositoryState(commit="1" * 40, dirty=False),
            )
        with self.assertRaisesRegex(ValidationError, "dirty-state mismatch"):
            validate_resume_request(
                frozen,
                plan,
                RepositoryState(commit=repository_state.commit, dirty=True),
            )
        with self.assertRaisesRegex(ValidationError, "fingerprint mismatch"):
            validate_resume_request(
                frozen,
                plan,
                RepositoryState(
                    commit=repository_state.commit,
                    dirty=False,
                    worktree_fingerprint="sha256:" + "b" * 64,
                ),
            )
        dirty_frozen = type(frozen)(
            run_id=frozen.run_id,
            profile_name=frozen.profile_name,
            profile_description=frozen.profile_description,
            filters=frozen.filters,
            config_hash=frozen.config_hash,
            repository=RepositoryState(
                commit=repository_state.commit,
                dirty=True,
                worktree_fingerprint="sha256:" + "a" * 64,
            ),
            cases=frozen.cases,
        )
        with self.assertRaisesRegex(ValidationError, "fingerprint mismatch"):
            validate_resume_request(
                dirty_frozen,
                plan,
                RepositoryState(
                    commit=repository_state.commit,
                    dirty=True,
                    worktree_fingerprint="sha256:" + "b" * 64,
                ),
            )
        legacy_dirty = type(frozen)(
            run_id=frozen.run_id,
            profile_name=frozen.profile_name,
            profile_description=frozen.profile_description,
            filters=frozen.filters,
            config_hash=frozen.config_hash,
            repository=RepositoryState(
                commit=repository_state.commit,
                dirty=True,
            ),
            cases=frozen.cases,
        )
        with self.assertRaisesRegex(ValidationError, "legacy dirty manifest"):
            validate_resume_request(
                legacy_dirty,
                plan,
                RepositoryState(
                    commit=repository_state.commit,
                    dirty=True,
                    worktree_fingerprint="sha256:" + "b" * 64,
                ),
            )
        changed_plan = type(plan)(
            profile=plan.profile,
            cases=plan.cases,
            filters=plan.filters,
            config_hash="f" * 64,
        )
        with self.assertRaisesRegex(ValidationError, "configuration mismatch"):
            validate_resume_request(frozen, changed_plan, repository_state)

    def test_resume_skips_only_complete_verified_results_and_retries_the_rest(self) -> None:
        repository, catalog, _, _, store, _, frozen = self._prepared_run(
            schemes=["dg", "pr", "be"]
        )
        executor = Executor(catalog, store, repository)
        first = executor.execute_frozen("resume-contract", frozen.cases)
        self.assertEqual((first.passed, first.skipped, first.failed), (3, 0, 0))

        directories = [
            repository / "benchmarks" / "resume-contract" / "cases" / case.case_id
            for case in frozen.cases
        ]
        verified = executor.verified_result(directories[0], frozen.cases[0])
        self.assertIsNotNone(verified)
        assert verified is not None
        self.assertEqual(verified["domain_result"]["checksum"], "fake-ok")
        # Missing result and an explicitly failed status must both be retried.
        (directories[1] / "result.json").unlink()
        failed_status = json.loads(
            (directories[2] / "status.json").read_text(encoding="utf-8")
        )
        failed_status["state"] = "failed"
        (directories[2] / "status.json").write_text(
            json.dumps(failed_status), encoding="utf-8"
        )

        resumed = executor.execute_frozen(
            "resume-contract", frozen.cases, resume=True
        )
        self.assertEqual((resumed.passed, resumed.skipped, resumed.failed), (2, 1, 0))
        for directory in directories:
            status = json.loads(
                (directory / "status.json").read_text(encoding="utf-8")
            )
            self.assertEqual(status["state"], "passed")
            self.assertTrue((directory / "result.json").is_file())

        # A torn/corrupt status is not trusted as a completed case either.
        (directories[0] / "status.json").write_text("{", encoding="utf-8")
        repaired = executor.execute_frozen(
            "resume-contract", frozen.cases, resume=True
        )
        self.assertEqual((repaired.passed, repaired.skipped, repaired.failed), (1, 2, 0))

        # A valid-looking persisted result is still untrusted when either log
        # reparses differently or the saved domain record disagrees with it.
        stdout = (directories[0] / "stdout.log").read_text(encoding="utf-8")
        (directories[0] / "stdout.log").write_text(
            stdout.replace("fake-ok", "forged-checksum"), encoding="utf-8"
        )
        result = json.loads(
            (directories[1] / "result.json").read_text(encoding="utf-8")
        )
        result["domain_result"]["checksum"] = "forged-result"
        (directories[1] / "result.json").write_text(
            json.dumps(result), encoding="utf-8"
        )
        with (directories[2] / "stderr.log").open("a", encoding="utf-8") as stream:
            stream.write("CORRUPTED_FAKE_STDERR\n")
        reparsed = executor.execute_frozen(
            "resume-contract", frozen.cases, resume=True
        )
        self.assertEqual((reparsed.passed, reparsed.skipped, reparsed.failed), (3, 0, 0))

        (directories[0] / "status.json").write_text(
            '{"oversized_integer":' + "9" * 5000 + "}",
            encoding="utf-8",
        )
        huge_integer = executor.execute_frozen(
            "resume-contract", frozen.cases, resume=True
        )
        self.assertEqual(
            (huge_integer.passed, huge_integer.skipped, huge_integer.failed),
            (1, 2, 0),
        )

        (directories[1] / "result.json").write_text(
            '{"nested":' + "[" * 2000 + "0" + "]" * 2000 + "}",
            encoding="utf-8",
        )
        recursive = executor.execute_frozen(
            "resume-contract", frozen.cases, resume=True
        )
        self.assertEqual(
            (recursive.passed, recursive.skipped, recursive.failed),
            (1, 2, 0),
        )

    def test_resume_requires_regular_requested_sample_file(self) -> None:
        repository, catalog, _, _, store, _, frozen = self._prepared_run(
            write_samples=True
        )
        executor = Executor(catalog, store, repository)
        first = executor.execute_frozen("resume-contract", frozen.cases)
        self.assertEqual((first.passed, first.failed), (1, 0))
        case_directory = (
            repository
            / "benchmarks"
            / "resume-contract"
            / "cases"
            / frozen.cases[0].case_id
        )
        samples = case_directory / "field_samples.csv"
        self.assertTrue(samples.is_file())

        outside = repository / "outside-samples.csv"
        outside.write_text("sentinel\n", encoding="utf-8")
        samples.unlink()
        samples.symlink_to(outside)
        resumed = executor.execute_frozen(
            "resume-contract", frozen.cases, resume=True
        )
        self.assertEqual((resumed.passed, resumed.skipped, resumed.failed), (1, 0, 0))
        self.assertFalse(samples.is_symlink())
        self.assertTrue(samples.is_file())
        self.assertEqual(outside.read_text(encoding="utf-8"), "sentinel\n")

    def test_tampered_command_payload_is_rejected_and_retried(self) -> None:
        repository, catalog, _, _, store, _, frozen = self._prepared_run()
        executor = Executor(catalog, store, repository)
        first = executor.execute_frozen("resume-contract", frozen.cases)
        self.assertEqual((first.passed, first.failed), (1, 0))
        case = frozen.cases[0]
        case_directory = (
            repository
            / "benchmarks"
            / "resume-contract"
            / "cases"
            / case.case_id
        )
        status = json.loads(
            (case_directory / "status.json").read_text(encoding="utf-8")
        )
        status["command"][-1] = "forged-payload-argument"
        (case_directory / "status.json").write_text(
            json.dumps(status), encoding="utf-8"
        )

        self.assertIsNone(executor.verified_result(case_directory, case))
        resumed = executor.execute_frozen(
            "resume-contract", frozen.cases, resume=True
        )
        self.assertEqual((resumed.passed, resumed.skipped, resumed.failed), (1, 0, 0))
        self.assertIsNotNone(executor.verified_result(case_directory, case))

    def test_resume_refuses_changed_launcher_prefix_before_retry(self) -> None:
        initial_launcher = PrefixLauncher(())
        repository, _, _, _, store, _, frozen = self._prepared_run(
            schemes=["dg", "pr"], launcher=initial_launcher
        )
        initial_catalog = build_catalog()
        initial_catalog.adapters.register("fake", FakeAdapter(mode="success"))
        initial_catalog.launchers.register(initial_launcher.name, initial_launcher)
        initial_executor = Executor(initial_catalog, store, repository)
        first = initial_executor.execute_frozen(
            "resume-contract", frozen.cases
        )
        self.assertEqual((first.passed, first.failed), (2, 0))

        directories = [
            repository / "benchmarks" / "resume-contract" / "cases" / case.case_id
            for case in frozen.cases
        ]
        # The first case is entirely absent; the second is reusable but was
        # executed with a different launcher prefix.  Refusal must not recreate
        # the missing directory or planned status as a side effect.
        for artifact in directories[0].iterdir():
            artifact.unlink()
        directories[0].rmdir()
        second_status = (directories[1] / "status.json").read_bytes()
        second_stdout = (directories[1] / "stdout.log").read_bytes()

        changed_launcher = PrefixLauncher(("/usr/bin/env",))
        changed_catalog = build_catalog()
        changed_catalog.adapters.register("fake", FakeAdapter(mode="success"))
        changed_catalog.launchers.register(changed_launcher.name, changed_launcher)
        changed_executor = Executor(changed_catalog, store, repository)
        # Offline verification deliberately remains payload-suffix-only.
        self.assertIsNotNone(
            changed_executor.verified_result(directories[1], frozen.cases[1])
        )
        observed: list[str] = []
        with self.assertRaisesRegex(
            ExecutionError, "resume launcher command mismatch"
        ):
            changed_executor.execute_frozen(
                "resume-contract",
                frozen.cases,
                resume=True,
                observer=lambda _position, _total, outcome: observed.append(
                    outcome.case.case_id
                ),
            )

        self.assertEqual(observed, [])
        self.assertFalse(directories[0].exists())
        self.assertEqual((directories[1] / "status.json").read_bytes(), second_status)
        self.assertEqual((directories[1] / "stdout.log").read_bytes(), second_stdout)

    def test_timeout_is_retried_instead_of_being_treated_as_complete(self) -> None:
        repository, catalog, _, _, store, _, frozen = self._prepared_run(
            timeout="0.05"
        )
        timeout_catalog = build_catalog()
        timeout_catalog.adapters.register("fake", FakeAdapter(mode="sleep"))
        timeout_executor = Executor(timeout_catalog, store, repository)
        case = frozen.cases[0]
        first = timeout_executor.execute_frozen("resume-contract", [case])
        self.assertEqual(first.outcomes[0].status["state"], "timeout")

        success_executor = Executor(catalog, store, repository)
        resumed = success_executor.execute_frozen(
            "resume-contract", [case], resume=True
        )
        self.assertEqual((resumed.passed, resumed.skipped, resumed.failed), (1, 0, 0))


if __name__ == "__main__":
    unittest.main(verbosity=2)
