from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
import io
import os
from pathlib import Path
import shlex
import tempfile
import unittest
from unittest.mock import patch

from ads_benchmark import cli
from ads_benchmark.framework.model import RepositoryState
from ads_benchmark.framework.sharding import (
    merge_shard_manifests,
    shard_subset_manifest,
    validate_shard_manifests,
)
from ads_benchmark.framework.storage import ResultStore
from benchmark_paths import CONFIG_DIRECTORY


class StageSevenShardCliTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory(prefix="ads-shard-cli-")
        self.addCleanup(temporary.cleanup)
        self.repository = Path(temporary.name) / "repository"
        self.repository.mkdir()
        self.repository_state = RepositoryState(
            commit="a" * 40,
            dirty=False,
            worktree_fingerprint=None,
        )

    def _selection_arguments(self) -> list[str]:
        return [
            "--repository-root",
            str(self.repository),
            "--config-dir",
            str(CONFIG_DIRECTORY),
            "--profile",
            "weak-scaling-smoke",
        ]

    def _main(self, arguments: list[str]) -> tuple[int, str, str]:
        stdout = io.StringIO()
        stderr = io.StringIO()
        with (
            patch.dict(os.environ, {"BENCHMARK_LAUNCHER_TEMPLATE": ""}),
            patch(
                "ads_benchmark.cli.inspect_repository",
                return_value=self.repository_state,
            ),
            redirect_stdout(stdout),
            redirect_stderr(stderr),
        ):
            return_code = cli.main(arguments)
        return return_code, stdout.getvalue(), stderr.getvalue()

    def _create_parent(self, run_id: str = "weak-parent") -> ResultStore:
        return_code, stdout, stderr = self._main(
            [
                "shard",
                *self._selection_arguments(),
                "--run-id",
                run_id,
                "--strategy",
                "index",
                "--shard-count",
                "2",
            ]
        )
        self.assertEqual(return_code, 0, stderr)
        self.assertEqual(stderr, "")
        self.assertIn("shards:       2", stdout)
        return ResultStore(self.repository)

    def test_plan_show_commands_expands_neutral_template_to_final_argv(self) -> None:
        template = (
            "srun --ntasks={ranks} --cpus-per-task={threads} "
            "--distribution=block:block --comment 'weak scaling' "
            "--grid={procx}x{procy}x{procz} {payload}"
        )
        return_code, stdout, stderr = self._main(
            [
                "plan",
                *self._selection_arguments(),
                "--launcher-template",
                template,
                "--show-commands",
            ]
        )

        self.assertEqual(return_code, 0, stderr)
        self.assertEqual(stderr, "")
        self.assertIn("mode:         dry-run (no files written)", stdout)
        lines = stdout.splitlines()
        command_start = lines.index("commands:") + 1
        rendered = [
            shlex.split(line.strip().split(": ", 1)[1])
            for line in lines[command_start:]
        ]
        self.assertEqual(len(rendered), 3)
        expected_payload = str(
            self.repository
            / "benchmarking"
            / "build"
            / "release"
            / "EXEC"
            / "igrm_l2_manufactured"
        )
        for argv in rendered:
            self.assertEqual(argv[0], "srun")
            self.assertIn("--cpus-per-task=1", argv)
            self.assertIn("weak scaling", argv)
            self.assertIn(expected_payload, argv)
            self.assertFalse(any("{" in argument for argument in argv))
        self.assertEqual(
            sum("--ntasks=2" in argv and "--grid=2x1x1" in argv for argv in rendered),
            1,
        )
        self.assertEqual(sum("--ntasks=1" in argv for argv in rendered), 2)
        self.assertFalse((self.repository / "benchmarks").exists())

    def test_shard_creates_parent_and_complete_index_shard_set(self) -> None:
        store = self._create_parent()
        parent = store.read_manifest("weak-parent")
        shards = store.read_all_shard_manifests("weak-parent")

        validated = validate_shard_manifests(shards)
        self.assertEqual(validated.shard_count, 2)
        self.assertEqual(validated.parent_manifest, parent)
        self.assertEqual(merge_shard_manifests(shards), parent)
        self.assertEqual(parent["case_count"], 3)
        self.assertEqual(sum(shard["case_count"] for shard in shards), 3)
        self.assertEqual(
            sorted(
                path.name
                for path in (
                    store.results_root / "weak-parent" / "shards"
                ).iterdir()
            ),
            ["shard-000000.json", "shard-000001.json"],
        )

    def test_run_shard_persists_and_executes_only_the_exact_subset(self) -> None:
        store = self._create_parent()
        shard = store.read_shard_manifest("weak-parent", 0)
        expected = shard_subset_manifest(shard, run_id="weak-child-0")

        with (
            patch.object(cli.Executor, "preflight", autospec=True) as preflight,
            patch("ads_benchmark.cli._execute_frozen_run", return_value=0) as execute,
        ):
            return_code, _, stderr = self._main(
                [
                    "run-shard",
                    "--repository-root",
                    str(self.repository),
                    "--parent-run-id",
                    "weak-parent",
                    "--shard-index",
                    "0",
                    "--run-id",
                    "weak-child-0",
                    "--available-mpi-slots",
                    "2",
                    "--available-cpu-slots",
                    "2",
                ]
            )

        self.assertEqual(return_code, 0, stderr)
        self.assertEqual(stderr, "")
        self.assertEqual(store.read_manifest("weak-child-0"), expected)
        self.assertLess(
            expected["case_count"],
            store.read_manifest("weak-parent")["case_count"],
        )
        preflight.assert_called_once()
        execute.assert_called_once()
        frozen = execute.call_args.args[1]
        self.assertEqual(frozen.run_id, "weak-child-0")
        self.assertEqual(
            [case.case_id for case in frozen.cases],
            [case["case_id"] for case in shard["cases"]],
        )
        self.assertEqual(execute.call_args.kwargs["parent_run_id"], "weak-parent")
        self.assertEqual(execute.call_args.kwargs["shard_index"], 0)

        wrong = shard_subset_manifest(
            store.read_shard_manifest("weak-parent", 1),
            run_id="wrong-child",
        )
        store.create_run("wrong-child", wrong)
        return_code, _, stderr = self._main(
            [
                "run-shard",
                "--repository-root",
                str(self.repository),
                "--parent-run-id",
                "weak-parent",
                "--shard-index",
                "0",
                "--run-id",
                "wrong-child",
                "--resume",
                "--available-mpi-slots",
                "2",
                "--available-cpu-slots",
                "2",
            ]
        )
        self.assertEqual(return_code, 2)
        self.assertIn("does not exactly match the generated subset", stderr)

    def test_merge_rejects_missing_run_then_accepts_all_exact_children(self) -> None:
        store = self._create_parent()
        shards = store.read_all_shard_manifests("weak-parent")
        child_ids: list[str] = []
        for index, shard in enumerate(shards):
            child_id = f"weak-child-{index}"
            child_ids.append(child_id)
            store.create_run(
                child_id,
                shard_subset_manifest(shard, run_id=child_id),
            )

        return_code, _, stderr = self._main(
            [
                "merge-shards",
                "--repository-root",
                str(self.repository),
                "--parent-run-id",
                "weak-parent",
                "--shard-run",
                child_ids[0],
            ]
        )
        self.assertEqual(return_code, 2)
        self.assertIn("missing completed shard runs for indices: 1", stderr)

        def fake_load(_executor, run_id, cases):
            return tuple(
                {"case_id": case.case_id, "child_run": run_id}
                for case in cases
            )

        report = object()
        written = object()
        arguments = [
            "merge-shards",
            "--repository-root",
            str(self.repository),
            "--parent-run-id",
            "weak-parent",
        ]
        for child_id in child_ids:
            arguments.extend(("--shard-run", child_id))
        with (
            patch(
                "ads_benchmark.cli.load_run_results", side_effect=fake_load
            ) as load,
            patch(
                "ads_benchmark.cli._analyze_results", return_value=report
            ) as analyze,
            patch(
                "ads_benchmark.cli.write_run_report", return_value=written
            ) as write,
            patch("ads_benchmark.cli._print_analysis", return_value=0) as display,
        ):
            return_code, stdout, stderr = self._main(arguments)

        self.assertEqual(return_code, 0, stderr)
        self.assertEqual(stderr, "")
        self.assertIn("shards:       2", stdout)
        self.assertIn("child runs:   2", stdout)
        self.assertEqual(load.call_count, 2)
        analyze.assert_called_once()
        analyzed_frozen = analyze.call_args.args[1]
        analyzed_results = analyze.call_args.args[2]
        self.assertEqual(len(analyzed_frozen.cases), 3)
        self.assertEqual(len(analyzed_results), 3)
        self.assertEqual(
            {result["child_run"] for result in analyzed_results},
            set(child_ids),
        )
        self.assertEqual(analyze.call_args.kwargs["source_run"], "weak-parent")
        write.assert_called_once()
        self.assertEqual(write.call_args.args[2], "weak-parent")
        display.assert_called_once_with("weak-parent", report, written)


if __name__ == "__main__":
    unittest.main()
