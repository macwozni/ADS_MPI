from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
from fractions import Fraction
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

from ads_benchmark import cli
from ads_benchmark.catalog import build_catalog
from ads_benchmark.components.manufactured import RESULT_PREFIX
from ads_benchmark.components.planning import MpiLauncher
from ads_benchmark.framework.config import load_profiles
from ads_benchmark.framework.executor import Executor
from ads_benchmark.framework.filtering import CaseFilters
from ads_benchmark.framework.model import ExecutionContext, RepositoryState
from ads_benchmark.framework.planner import Planner
from ads_benchmark.framework.storage import ResultStore
from fake_adapter import FakeAdapter


CONFIG_DIRECTORY = Path(__file__).resolve().parents[1] / "configs"


def fake_profile() -> dict[str, object]:
    return {
        "schema_version": 1,
        "name": "cli-profile",
        "description": "one executable CLI contract case",
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
        "execution": {"warmups": 0, "samples": 1, "timeout_seconds": "2"},
        "build_profiles": ["debug"],
        "launcher": "direct",
    }


class RunCommandTests(unittest.TestCase):
    def test_run_requires_an_explicit_run_id(self) -> None:
        stderr = io.StringIO()
        with redirect_stderr(stderr), self.assertRaises(SystemExit) as raised:
            cli.main(["run"])
        self.assertEqual(raised.exception.code, 2)
        self.assertIn("--run-id", stderr.getvalue())

    def test_steps_filter_rejects_nonpositive_values(self) -> None:
        temporary = tempfile.TemporaryDirectory(prefix="ads-cli-steps-")
        self.addCleanup(temporary.cleanup)
        repository = Path(temporary.name) / "repository"
        repository.mkdir()
        configuration = repository / "profiles"
        configuration.mkdir()
        (configuration / "cli-profile.json").write_text(
            json.dumps(fake_profile()), encoding="utf-8"
        )
        stderr = io.StringIO()
        with redirect_stderr(stderr):
            return_code = cli.main(
                [
                    "plan",
                    "--repository-root",
                    str(repository),
                    "--config-dir",
                    str(configuration),
                    "--profile",
                    "cli-profile",
                    "--steps",
                    "0",
                ]
            )
        self.assertEqual(return_code, 2)
        self.assertIn("time steps filters must be positive", stderr.getvalue())

    def _run_fake(
        self,
        mode: str,
        run_id: str,
        *,
        warmups: int = 0,
        samples: int = 1,
    ) -> tuple[int, Path, str, str]:
        temporary = tempfile.TemporaryDirectory(prefix=f"ads-cli-{mode}-")
        self.addCleanup(temporary.cleanup)
        repository = Path(temporary.name) / "repository"
        repository.mkdir()
        configuration = repository / "profiles"
        configuration.mkdir()
        profile = fake_profile()
        profile["execution"] = {
            "warmups": warmups,
            "samples": samples,
            "timeout_seconds": "2",
        }
        (configuration / "cli-profile.json").write_text(
            json.dumps(profile), encoding="utf-8"
        )

        catalog = build_catalog()
        catalog.adapters.register("fake", FakeAdapter(mode=mode))
        stdout = io.StringIO()
        stderr = io.StringIO()
        arguments = [
            "run",
            "--repository-root",
            str(repository),
            "--config-dir",
            str(configuration),
            "--profile",
            "cli-profile",
            "--run-id",
            run_id,
            "--steps",
            "4",
        ]
        with (
            patch("ads_benchmark.cli.build_catalog", return_value=catalog),
            patch(
                "ads_benchmark.cli.inspect_repository",
                return_value=RepositoryState(commit="0" * 40, dirty=False),
            ),
            redirect_stdout(stdout),
            redirect_stderr(stderr),
        ):
            return_code = cli.main(arguments)
        return return_code, repository, stdout.getvalue(), stderr.getvalue()

    def test_run_creates_manifest_and_executes_every_selected_case(self) -> None:
        return_code, repository, stdout, stderr = self._run_fake(
            "success", "cli-success"
        )
        self.assertEqual(return_code, 0, stderr)
        self.assertEqual(stderr, "")
        self.assertIn("summary:      passed=1 failed=0 total=1", stdout)

        run_directory = repository / "benchmarks" / "cli-success"
        manifest = json.loads(
            (run_directory / "manifest.json").read_text(encoding="utf-8")
        )
        self.assertEqual(manifest["run_id"], "cli-success")
        self.assertEqual(manifest["case_count"], 1)
        self.assertEqual(manifest["filters"]["steps"], [4])
        case_id = manifest["cases"][0]["case_id"]
        case_directory = run_directory / "cases" / case_id
        status = json.loads(
            (case_directory / "status.json").read_text(encoding="utf-8")
        )
        self.assertEqual(status["state"], "passed")
        self.assertTrue((case_directory / "result.json").is_file())

    def test_run_returns_failure_when_a_case_fails(self) -> None:
        return_code, repository, stdout, stderr = self._run_fake(
            "nonzero", "cli-failure"
        )
        self.assertEqual(return_code, 1, stderr)
        self.assertIn("summary:      passed=0 failed=1 total=1", stdout)
        self.assertIn("process exited with status 7", stderr)

        run_directory = repository / "benchmarks" / "cli-failure"
        manifest = json.loads(
            (run_directory / "manifest.json").read_text(encoding="utf-8")
        )
        case_id = manifest["cases"][0]["case_id"]
        case_directory = run_directory / "cases" / case_id
        status = json.loads(
            (case_directory / "status.json").read_text(encoding="utf-8")
        )
        self.assertEqual(status["state"], "failed")
        self.assertFalse((case_directory / "result.json").exists())

    def test_execution_preflight_refuses_missing_payload_before_run_creation(self) -> None:
        return_code, repository, _, stderr = self._run_fake(
            "missing-executable", "missing-payload"
        )
        self.assertEqual(return_code, 2)
        self.assertIn("execution preflight", stderr)
        self.assertFalse((repository / "benchmarks" / "missing-payload").exists())

    def test_execution_preflight_refuses_unimplemented_repetitions(self) -> None:
        return_code, repository, _, stderr = self._run_fake(
            "success",
            "unsupported-repetitions",
            warmups=2,
            samples=7,
        )
        self.assertEqual(return_code, 2)
        self.assertIn("supports only warmups=0 and samples=1", stderr)
        self.assertIn("requests warmups=2 and samples=7", stderr)
        self.assertFalse(
            (repository / "benchmarks" / "unsupported-repetitions").exists()
        )

    def test_execution_preflight_refuses_missing_launcher_before_run_creation(self) -> None:
        temporary = tempfile.TemporaryDirectory(prefix="ads-cli-launcher-")
        self.addCleanup(temporary.cleanup)
        repository = Path(temporary.name) / "repository"
        repository.mkdir()
        configuration = repository / "profiles"
        configuration.mkdir()
        profile = fake_profile()
        profile["launcher"] = "missing-mpi"
        (configuration / "cli-profile.json").write_text(
            json.dumps(profile), encoding="utf-8"
        )
        catalog = build_catalog()
        catalog.adapters.register("fake", FakeAdapter(mode="success"))
        catalog.launchers.register(
            "missing-mpi",
            MpiLauncher(
                executable="/definitely/missing/ads-mpiexec",
                rank_flag="-n",
                name="missing-mpi",
            ),
        )
        stderr = io.StringIO()
        with (
            patch("ads_benchmark.cli.build_catalog", return_value=catalog),
            patch(
                "ads_benchmark.cli.inspect_repository",
                return_value=RepositoryState(commit="0" * 40, dirty=False),
            ),
            redirect_stdout(io.StringIO()),
            redirect_stderr(stderr),
        ):
            return_code = cli.main(
                [
                    "run",
                    "--repository-root",
                    str(repository),
                    "--config-dir",
                    str(configuration),
                    "--profile",
                    "cli-profile",
                    "--run-id",
                    "missing-launcher",
                ]
            )
        self.assertEqual(return_code, 2)
        self.assertIn("launcher missing-mpi", stderr.getvalue())
        self.assertFalse((repository / "benchmarks" / "missing-launcher").exists())

    def test_dry_run_plan_does_not_require_execution_preflight(self) -> None:
        temporary = tempfile.TemporaryDirectory(prefix="ads-cli-dry-plan-")
        self.addCleanup(temporary.cleanup)
        repository = Path(temporary.name) / "repository"
        repository.mkdir()
        configuration = repository / "profiles"
        configuration.mkdir()
        (configuration / "cli-profile.json").write_text(
            json.dumps(fake_profile()), encoding="utf-8"
        )
        catalog = build_catalog()
        catalog.adapters.register("fake", FakeAdapter(mode="missing-executable"))
        with (
            patch("ads_benchmark.cli.build_catalog", return_value=catalog),
            patch(
                "ads_benchmark.cli.inspect_repository",
                return_value=RepositoryState(commit="0" * 40, dirty=False),
            ),
            patch.object(
                Executor,
                "preflight",
                side_effect=AssertionError("dry-run called execution preflight"),
            ),
            redirect_stdout(io.StringIO()),
            redirect_stderr(io.StringIO()),
        ):
            self.assertEqual(
                cli.main(
                    [
                        "plan",
                        "--repository-root",
                        str(repository),
                        "--config-dir",
                        str(configuration),
                        "--profile",
                        "cli-profile",
                        "--dry-run",
                    ]
                ),
                0,
            )

    def test_resume_command_uses_frozen_manifest_and_skips_verified_result(self) -> None:
        temporary = tempfile.TemporaryDirectory(prefix="ads-cli-resume-")
        self.addCleanup(temporary.cleanup)
        repository = Path(temporary.name) / "repository"
        repository.mkdir()
        configuration = repository / "profiles"
        configuration.mkdir()
        (configuration / "cli-profile.json").write_text(
            json.dumps(fake_profile()), encoding="utf-8"
        )
        catalog = build_catalog()
        catalog.adapters.register("fake", FakeAdapter(mode="success"))
        common = [
            "--repository-root",
            str(repository),
            "--config-dir",
            str(configuration),
            "--profile",
            "cli-profile",
            "--run-id",
            "cli-resume",
        ]

        with (
            patch("ads_benchmark.cli.build_catalog", return_value=catalog),
            patch(
                "ads_benchmark.cli.inspect_repository",
                return_value=RepositoryState(commit="0" * 40, dirty=False),
            ),
            redirect_stdout(io.StringIO()),
            redirect_stderr(io.StringIO()),
        ):
            self.assertEqual(cli.main(["run", *common]), 0)

        stdout = io.StringIO()
        stderr = io.StringIO()
        with (
            patch("ads_benchmark.cli.build_catalog", return_value=catalog),
            patch(
                "ads_benchmark.cli.inspect_repository",
                return_value=RepositoryState(commit="0" * 40, dirty=False),
            ),
            redirect_stdout(stdout),
            redirect_stderr(stderr),
        ):
            self.assertEqual(cli.main(["resume", *common]), 0)
        self.assertEqual(stderr.getvalue(), "")
        self.assertIn("passed (verified, skipped)", stdout.getvalue())
        self.assertIn("skipped=1 retried=0", stdout.getvalue())


class AnalyzeCommandTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory(prefix="ads-cli-analyze-")
        self.addCleanup(temporary.cleanup)
        self.repository = Path(temporary.name) / "repository"
        self.repository.mkdir()

    def _create_run(self, run_id: str) -> tuple[ResultStore, tuple[object, ...]]:
        catalog = build_catalog()
        filters = CaseFilters(
            problems=frozenset({"igrm_l2"}),
            schemes=frozenset({"dg"}),
            degree_pairs=frozenset(
                {((4, 4, 4), (3, 3, 3))}
            ),
            steps=frozenset({4, 8, 16, 32}),
        )
        plan = Planner(load_profiles(CONFIG_DIRECTORY), catalog).plan(
            "temporal-full", filters
        )
        store = ResultStore(self.repository)
        store.create_run(
            run_id,
            plan.manifest(
                run_id=run_id,
                repository=RepositoryState(commit="0" * 40, dirty=False),
            ),
        )
        for case in plan.cases:
            case_directory = store.create_case_directory(run_id, case.case_id)
            context = ExecutionContext(
                repository_root=self.repository,
                case_directory=case_directory,
            )
            payload = catalog.adapters.get(case.spec.problem).build_payload_command(
                case.spec, context
            )
            command = catalog.launchers.get(case.spec.launcher).command(
                payload, case.spec
            )
            dt = float(Fraction(case.spec.time.time_step))
            error = dt * dt
            domain_result = {
                "schema_version": 1,
                "kind": "ads-manufactured-transient-result",
                "problem": case.spec.problem,
                "scheme": case.spec.scheme,
                "exact_case": case.spec.exact_case,
                "requested_final_time": 0.1,
                "actual_final_time": 0.1,
                "time_step": dt,
                "steps": case.spec.time.steps,
                "initial_l2_error": 1.0e-15,
                "initial_linf_error": 2.0e-15,
                "l2_error": error,
                "linf_error": 2.0 * error,
                "solution_l2_norm": 0.2,
                "field_checksum": 1.0,
                "sample_points_per_axis": case.spec.sampling.points_per_axis,
                "field_samples_written": case.spec.sampling.write_samples,
                "physical_step_wall_seconds": 0.5,
                "solver_status": 0,
            }
            store.write_status(
                case_directory,
                {
                    "schema_version": 1,
                    "case_id": case.case_id,
                    "state": "passed",
                    "started_at": "2026-01-01T00:00:00+00:00",
                    "finished_at": "2026-01-01T00:00:01+00:00",
                    "duration_seconds": 1.0,
                    "return_code": 0,
                    "command": list(command),
                },
            )
            store.write_log(
                case_directory,
                "stdout.log",
                RESULT_PREFIX + json.dumps(domain_result, sort_keys=True) + "\n",
            )
            store.write_log(case_directory, "stderr.log", "")
            store.write_result(
                case_directory,
                {
                    "schema_version": 1,
                    "kind": "ads-benchmark-case-result",
                    "case_id": case.case_id,
                    "status": "passed",
                    "configuration": case.spec.to_dict(),
                    "timing": {"wall_seconds": 1.0},
                    "domain_result": domain_result,
                },
            )
        return store, plan.cases

    def test_analyze_dispatches_registered_family_analyzer(self) -> None:
        self._create_run("filtered-four-levels")
        stdout = io.StringIO()
        stderr = io.StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            return_code = cli.main(
                [
                    "analyze",
                    "--repository-root",
                    str(self.repository),
                    "--run-id",
                    "filtered-four-levels",
                ]
            )
        self.assertEqual(return_code, 0, stderr.getvalue())
        self.assertEqual(stderr.getvalue(), "")
        self.assertIn("analyzer:     temporal-convergence", stdout.getvalue())
        output = self.repository / "benchmarks" / "filtered-four-levels" / "analysis"
        report = json.loads((output / "analysis.json").read_text(encoding="utf-8"))
        self.assertEqual(report["status"], "passed")
        self.assertEqual(report["summary"]["levels_per_series"], 4)
        self.assertTrue((output / "analysis.csv").is_file())

    def test_analyze_rejects_result_that_disagrees_with_tagged_stdout(self) -> None:
        self._create_run("forged-result")
        run = self.repository / "benchmarks" / "forged-result"
        result_path = next((run / "cases").glob("*/result.json"))
        result = json.loads(result_path.read_text(encoding="utf-8"))
        result["domain_result"]["l2_error"] *= 10.0
        result_path.write_text(json.dumps(result), encoding="utf-8")
        stderr = io.StringIO()
        with redirect_stdout(io.StringIO()), redirect_stderr(stderr):
            return_code = cli.main(
                [
                    "analyze",
                    "--repository-root",
                    str(self.repository),
                    "--run-id",
                    "forged-result",
                ]
            )
        self.assertEqual(return_code, 2)
        self.assertIn("not a fully verified completed result", stderr.getvalue())
        self.assertFalse((run / "analysis").exists())

    def test_analyze_refuses_symlinked_output_directory(self) -> None:
        self._create_run("unsafe-analysis")
        outside = self.repository / "outside"
        outside.mkdir()
        run = self.repository / "benchmarks" / "unsafe-analysis"
        (run / "analysis").symlink_to(outside, target_is_directory=True)
        stderr = io.StringIO()
        with redirect_stdout(io.StringIO()), redirect_stderr(stderr):
            return_code = cli.main(
                [
                    "analyze",
                    "--repository-root",
                    str(self.repository),
                    "--run-id",
                    "unsafe-analysis",
                ]
            )
        self.assertEqual(return_code, 2)
        self.assertIn("unsafe analysis directory", stderr.getvalue())
        self.assertEqual(list(outside.iterdir()), [])

    def test_analyze_removes_stale_plot_without_following_symlink(self) -> None:
        self._create_run("stale-plot")
        analysis = self.repository / "benchmarks" / "stale-plot" / "analysis"
        analysis.mkdir()
        outside = self.repository / "outside.png"
        outside.write_bytes(b"sentinel")
        plot = analysis / "convergence.png"
        plot.symlink_to(outside)

        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            self.assertEqual(
                cli.main(
                    [
                        "analyze",
                        "--repository-root",
                        str(self.repository),
                        "--run-id",
                        "stale-plot",
                    ]
                ),
                0,
            )
        self.assertFalse(plot.exists())
        self.assertEqual(outside.read_bytes(), b"sentinel")

        plot.symlink_to(outside)
        stdout = io.StringIO()
        with (
            patch.dict(sys.modules, {"matplotlib": None}),
            redirect_stdout(stdout),
            redirect_stderr(io.StringIO()),
        ):
            self.assertEqual(
                cli.main(
                    [
                        "analyze",
                        "--repository-root",
                        str(self.repository),
                        "--run-id",
                        "stale-plot",
                        "--plot",
                    ]
                ),
                0,
            )
        self.assertFalse(plot.exists())
        self.assertEqual(outside.read_bytes(), b"sentinel")
        self.assertIn("plot omitted", stdout.getvalue())


if __name__ == "__main__":
    unittest.main(verbosity=2)
