from __future__ import annotations

import json
from pathlib import Path
import subprocess
import tempfile
import unittest

from ads_benchmark.framework.errors import ExecutionError, ValidationError
from ads_benchmark.orchestration import (
    OwnedWorktreeWorkspace,
    alternating_schedule,
    resolve_commit,
)


def git(repository: Path, *arguments: str) -> str:
    completed = subprocess.run(
        ("git", "-C", str(repository), *arguments),
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return completed.stdout.strip()


class OrchestrationTests(unittest.TestCase):
    def _repository(self) -> tuple[tempfile.TemporaryDirectory[str], Path, str, str]:
        temporary = tempfile.TemporaryDirectory(prefix="ads-ab-orchestrator-")
        self.addCleanup(temporary.cleanup)
        repository = Path(temporary.name) / "repository"
        repository.mkdir()
        git(repository, "init", "-q")
        git(repository, "config", "user.name", "Benchmark Test")
        git(repository, "config", "user.email", "benchmark@example.invalid")
        source = repository / "source.txt"
        source.write_text("first\n", encoding="utf-8")
        git(repository, "add", "source.txt")
        git(repository, "commit", "-q", "-m", "first")
        first = git(repository, "rev-parse", "HEAD")
        source.write_text("second\n", encoding="utf-8")
        git(repository, "commit", "-q", "-am", "second")
        second = git(repository, "rev-parse", "HEAD")
        return temporary, repository, first, second

    def test_schedule_alternates_first_side_and_rejects_duplicates(self) -> None:
        schedule = alternating_schedule(("case-a", "case-b", "case-c"))
        self.assertEqual(
            [entry.order for entry in schedule],
            [
                ("baseline", "candidate"),
                ("candidate", "baseline"),
                ("baseline", "candidate"),
            ],
        )
        self.assertEqual(schedule[0].to_dict()["order"], ["baseline", "candidate"])
        with self.assertRaisesRegex(ValidationError, "duplicate"):
            alternating_schedule(("same", "same"))
        with self.assertRaisesRegex(ValidationError, "at least one"):
            alternating_schedule(())

    def test_refs_resolve_to_full_commits_and_reject_option_injection(self) -> None:
        _, repository, first, second = self._repository()
        self.assertEqual(resolve_commit(repository, "HEAD"), second)
        self.assertEqual(resolve_commit(repository, "HEAD~1"), first)
        with self.assertRaisesRegex(ValidationError, "invalid Git reference"):
            resolve_commit(repository, "--help")

    def test_workspace_parent_inside_repository_is_refused(self) -> None:
        _, repository, _, _ = self._repository()
        owner = OwnedWorktreeWorkspace(
            repository,
            run_id="inside-repository",
            baseline_ref="HEAD",
            candidate_ref="HEAD",
            parent=repository,
        )
        with self.assertRaisesRegex(ExecutionError, "outside"):
            owner.prepare()
        worktrees = [
            line
            for line in git(
                repository, "worktree", "list", "--porcelain"
            ).splitlines()
            if line.startswith("worktree ")
        ]
        self.assertEqual(worktrees, [f"worktree {repository}"])

    def test_owned_worktrees_preserve_dirty_main_tree_and_cleanup_selectively(self) -> None:
        temporary, repository, first, second = self._repository()
        source = repository / "source.txt"
        source.write_text("dirty main tree\n", encoding="utf-8")
        untracked = repository / "untracked.txt"
        untracked.write_text("keep me\n", encoding="utf-8")
        before_status = git(repository, "status", "--porcelain=v1")
        before_source = source.read_bytes()
        before_untracked = untracked.read_bytes()

        owner = OwnedWorktreeWorkspace(
            repository,
            run_id="same-machine",
            baseline_ref=first,
            candidate_ref=second,
            parent=Path(temporary.name),
        )
        pair = owner.prepare()
        self.assertEqual(git(pair.baseline, "rev-parse", "HEAD"), first)
        self.assertEqual(git(pair.candidate, "rev-parse", "HEAD"), second)
        marker = json.loads(pair.marker.read_text(encoding="utf-8"))
        self.assertEqual(marker["baseline_commit"], first)
        self.assertEqual(marker["candidate_commit"], second)
        self.assertNotEqual(pair.baseline, pair.candidate)
        owner.cleanup()

        self.assertFalse(pair.workspace.exists())
        self.assertEqual(git(repository, "status", "--porcelain=v1"), before_status)
        self.assertEqual(source.read_bytes(), before_source)
        self.assertEqual(untracked.read_bytes(), before_untracked)

    def test_same_ref_uses_two_distinct_detached_worktrees(self) -> None:
        temporary, repository, _, second = self._repository()
        owner = OwnedWorktreeWorkspace(
            repository,
            run_id="identical-builds",
            baseline_ref=second,
            candidate_ref=second,
            parent=Path(temporary.name),
        )
        pair = owner.prepare()
        self.assertNotEqual(pair.baseline, pair.candidate)
        self.assertEqual(pair.baseline_commit, second)
        self.assertEqual(pair.candidate_commit, second)
        self.assertEqual(git(pair.baseline, "rev-parse", "HEAD"), second)
        self.assertEqual(git(pair.candidate, "rev-parse", "HEAD"), second)
        owner.cleanup()

    def test_cleanup_refuses_tampered_marker_or_modified_worktree(self) -> None:
        temporary, repository, first, second = self._repository()
        owner = OwnedWorktreeWorkspace(
            repository,
            run_id="tamper",
            baseline_ref=first,
            candidate_ref=second,
            parent=Path(temporary.name),
        )
        pair = owner.prepare()
        original_marker = pair.marker.read_text(encoding="utf-8")
        pair.marker.write_text("{}\n", encoding="utf-8")
        with self.assertRaisesRegex(ExecutionError, "invalid marker"):
            owner.cleanup()
        pair.marker.write_text(original_marker, encoding="utf-8")
        (pair.baseline / "foreign.txt").write_text("do not delete\n", encoding="utf-8")
        with self.assertRaisesRegex(ExecutionError, "modified"):
            owner.cleanup()
        self.assertTrue(pair.workspace.is_dir())
        self.assertTrue((pair.baseline / "foreign.txt").is_file())

        # The temporary test root owns this disposable repository.  Remove the
        # test modification through Git's normal worktree command only after
        # proving the production cleanup refused it.
        (pair.baseline / "foreign.txt").unlink()
        owner.cleanup()


if __name__ == "__main__":
    unittest.main(verbosity=2)
