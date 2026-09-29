from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from ads_benchmark import cli
from ads_benchmark.catalog import build_catalog
from ads_benchmark.framework.model import RepositoryState
from fake_adapter import FakeAdapter


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

    def _run_fake(self, mode: str, run_id: str) -> tuple[int, Path, str, str]:
        temporary = tempfile.TemporaryDirectory(prefix=f"ads-cli-{mode}-")
        self.addCleanup(temporary.cleanup)
        repository = Path(temporary.name) / "repository"
        repository.mkdir()
        configuration = repository / "profiles"
        configuration.mkdir()
        (configuration / "cli-profile.json").write_text(
            json.dumps(fake_profile()), encoding="utf-8"
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


if __name__ == "__main__":
    unittest.main(verbosity=2)
