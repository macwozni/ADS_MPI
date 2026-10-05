from __future__ import annotations

import csv
from contextlib import redirect_stdout
from datetime import datetime, timezone
from fractions import Fraction
import io
import json
import math
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from ads_benchmark import cli
from ads_benchmark.framework.errors import StorageError
from ads_benchmark.framework.provenance import collect_execution_provenance


def _git(repository: Path, *arguments: str) -> str:
    completed = subprocess.run(
        ("git", "-C", str(repository), *arguments),
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return completed.stdout.strip()


def _profile() -> dict[str, object]:
    """Two cases are enough to exercise both AB and BA schedule entries."""

    return {
        "schema_version": 1,
        "name": "ab-cli-test",
        "description": "Hermetic two-case A/B orchestration test",
        "family": "strong",
        "exact_cases": ["spatial-cosine"],
        "problems": ["igrm_l2"],
        "schemes": ["dg", "pr"],
        "time_discretizations": [{"final_time": "0.1", "steps": 4}],
        "meshes": [[2, 2, 2]],
        "degree_pairs": [
            {"test": [2, 2, 2], "trial": [1, 1, 1]},
        ],
        "process_layouts": [{"ranks": 1, "grid": [1, 1, 1]}],
        "thread_counts": [1],
        "openmp": {
            "dynamic": False,
            "proc_bind": "close",
            "places": "cores",
        },
        "sampling": {"points_per_axis": 3, "write_samples": True},
        "execution": {
            "warmups": 1,
            "samples": 5,
            "minimum_sample_seconds": "0.001",
            "timeout_seconds": "60",
        },
        "build_profiles": ["release"],
        "launcher": "mpi",
    }


def _checksum(values: list[float]) -> float:
    ordinary = math.fsum(values)
    weighted = math.fsum(
        (index + 1) * value for index, value in enumerate(values)
    )
    return ordinary + weighted / (len(values) + 1)


def _result(case: object) -> dict[str, object]:
    spec = case.spec
    final_time = float(Fraction(spec.time.final_time))
    coordinates: list[tuple[float, float, float]] = []
    numerical: list[float] = []
    for iz in range(3):
        z = iz / 2
        for iy in range(3):
            y = iy / 2
            for ix in range(3):
                x = ix / 2
                coordinates.append((x, y, z))
                numerical.append(
                    math.exp(-final_time)
                    * math.cos(math.pi * x)
                    * math.cos(math.pi * y)
                    * math.cos(math.pi * z)
                )

    stream = io.StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(("x", "y", "z", "numerical", "exact", "error"))
    for coordinate, value in zip(coordinates, numerical, strict=True):
        writer.writerow(
            tuple(format(component, ".17g") for component in coordinate)
            + (format(value, ".17g"), format(value, ".17g"), "0")
        )

    samples = [1.0] * spec.measurement.samples
    return {
        "schema_version": 1,
        "kind": "ads-benchmark-case-result",
        "case_id": case.case_id,
        "status": "passed",
        "configuration": spec.to_dict(),
        "timing": {
            "wall_seconds": math.fsum(samples),
            "metric": "physical_step_wall_seconds",
            "measured_samples": samples,
            "reliable": True,
        },
        "domain_result": {
            "problem": spec.problem,
            "scheme": spec.scheme,
            "exact_case": spec.exact_case,
            "solver_status": 0,
            "steps": spec.time.steps,
            "time_step": float(Fraction(spec.time.time_step)),
            "requested_final_time": final_time,
            "actual_final_time": final_time,
            "sample_points_per_axis": 3,
            "field_samples_written": True,
            "l2_error": 0.0,
            "linf_error": 0.0,
            "solution_l2_norm": math.sqrt(
                math.fsum(value * value for value in numerical)
            ),
            "field_checksum": _checksum(numerical),
            "physical_step_wall_seconds": samples[-1],
        },
        "analysis_artifacts": {"field_samples_csv": stream.getvalue()},
    }


class _FakeExecutor:
    calls: list[tuple[str, str, str]] = []

    def __init__(self, catalog: object, store: object, repository_root: Path) -> None:
        self.catalog = catalog
        self.store = store
        self.repository_root = repository_root
        self.side = repository_root.name

    def preflight(self, cases: object) -> None:
        tuple(cases)

    def execute_new_case_locked(
        self,
        run_id: str,
        case: object,
        *,
        max_retries: int = 0,
    ) -> dict[str, object]:
        self.calls.append((self.side, case.case_id, run_id))
        return {"state": "passed"}


class CompareRefsCliTests(unittest.TestCase):
    def test_cli_orchestrates_pair_schedule_report_and_cleanup(self) -> None:
        with tempfile.TemporaryDirectory(prefix="ads-ab-cli-test-") as temporary:
            root = Path(temporary)
            repository = root / "repository"
            configs = repository / "benchmarking" / "configs"
            configs.mkdir(parents=True)
            (configs / "ab-cli-test.json").write_text(
                json.dumps(_profile(), indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            _git(repository, "init", "-q")
            _git(repository, "config", "user.name", "Benchmark Test")
            _git(repository, "config", "user.email", "benchmark@example.invalid")
            _git(repository, "add", "benchmarking/configs/ab-cli-test.json")
            _git(repository, "commit", "-q", "-m", "benchmark profile")

            _FakeExecutor.calls = []

            def execution_record(
                repository_root: Path,
                catalog: object,
                manifest: dict[str, object],
                cases: object,
            ) -> dict[str, object]:
                # Missing fake binaries remain explicit observations.  A fixed,
                # empty environment avoids probing a host compiler/launcher.
                return collect_execution_provenance(
                    repository_root,
                    manifest,
                    environment={},
                    recorded_at=datetime(2026, 1, 2, tzinfo=timezone.utc),
                )

            def loaded_results(
                executor: _FakeExecutor,
                run_id: str,
                cases: object,
            ) -> tuple[dict[str, object], ...]:
                with self.assertRaisesRegex(StorageError, "already being"):
                    with executor.store.execution_lock(run_id):
                        pass
                return tuple(_result(case) for case in cases)

            with (
                patch.object(cli, "Executor", _FakeExecutor),
                patch.object(cli, "_build_ref_cases") as build,
                patch.object(
                    cli, "_execution_record", side_effect=execution_record
                ) as provenance,
                patch.object(
                    cli, "load_run_results", side_effect=loaded_results
                ),
            ):
                output = io.StringIO()
                with redirect_stdout(output):
                    return_code = cli.main(
                        [
                            "compare-refs",
                            "--repository-root",
                            str(repository),
                            "--profile",
                            "ab-cli-test",
                            "--run-id",
                            "hermetic-ab",
                            "--baseline-ref",
                            "HEAD",
                            "--candidate-ref",
                            "HEAD",
                            "--workspace-parent",
                            str(root),
                            "--available-mpi-slots",
                            "1",
                            "--available-cpu-slots",
                            "1",
                        ]
                    )

            self.assertEqual(return_code, 0)
            self.assertIn("status:       no-regression-detected", output.getvalue())
            self.assertEqual(build.call_count, 2)
            self.assertEqual(provenance.call_count, 2)

            baseline = repository / "benchmarks" / "hermetic-ab-baseline"
            candidate = repository / "benchmarks" / "hermetic-ab-candidate"
            self.assertTrue((baseline / "manifest.json").is_file())
            self.assertTrue((baseline / "execution.json").is_file())
            self.assertTrue((candidate / "manifest.json").is_file())
            self.assertTrue((candidate / "execution.json").is_file())

            baseline_schedule = json.loads(
                (baseline / "analysis" / "ab-schedule.json").read_text(
                    encoding="utf-8"
                )
            )
            candidate_schedule = json.loads(
                (candidate / "analysis" / "ab-schedule.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(baseline_schedule, candidate_schedule)
            entries = baseline_schedule["entries"]
            self.assertEqual(
                [entry["order"] for entry in entries],
                [["baseline", "candidate"], ["candidate", "baseline"]],
            )
            expected_calls = [
                (
                    side,
                    entry["case_id"],
                    f"hermetic-ab-{side}",
                )
                for entry in entries
                for side in entry["order"]
            ]
            self.assertEqual(_FakeExecutor.calls, expected_calls)

            report = json.loads(
                (candidate / "analysis" / "comparison.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(report["status"], "no-regression-detected")
            self.assertEqual(report["summary"]["case_count"], 2)
            self.assertEqual(report["summary"]["compared_timing_count"], 2)
            self.assertTrue(
                (candidate / "analysis" / "comparison.csv").is_file()
            )

            self.assertEqual(
                list(root.glob("ads-ab-hermetic-ab-*")),
                [],
                "successful orchestration must remove its owned workspace",
            )
            worktrees = [
                line
                for line in _git(repository, "worktree", "list", "--porcelain").splitlines()
                if line.startswith("worktree ")
            ]
            self.assertEqual(worktrees, [f"worktree {repository}"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
