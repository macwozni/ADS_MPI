from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile
import unittest

from ads_benchmark.framework.errors import ProvenanceError
from ads_benchmark.framework.provenance import inspect_repository


def git(repository: Path, *arguments: str) -> None:
    subprocess.run(
        ["git", "-C", str(repository), *arguments],
        check=True,
        capture_output=True,
        text=True,
    )


def initialize_repository(parent: Path) -> Path:
    repository = parent / "repository"
    repository.mkdir()
    git(repository, "init", "-q")
    git(repository, "config", "user.name", "Benchmark Test")
    git(repository, "config", "user.email", "benchmark@example.invalid")
    (repository / ".gitignore").write_text("/benchmarks/\n", encoding="utf-8")
    (repository / "tracked.txt").write_text("original\n", encoding="utf-8")
    git(repository, "add", ".gitignore", "tracked.txt")
    git(repository, "commit", "-q", "-m", "initial")
    return repository


class ProvenanceTests(unittest.TestCase):
    def test_clean_and_ignored_benchmark_outputs_are_not_dirty(self) -> None:
        with tempfile.TemporaryDirectory(prefix="ads-provenance-clean-") as temporary:
            repository = initialize_repository(Path(temporary))
            clean = inspect_repository(repository)
            self.assertFalse(clean.dirty)
            self.assertEqual(len(clean.commit), 40)
            output = repository / "benchmarks" / "run-one"
            output.mkdir(parents=True)
            (output / "manifest.json").write_text("{}\n", encoding="utf-8")
            self.assertFalse(inspect_repository(repository).dirty)

    def test_every_relevant_worktree_change_sets_dirty(self) -> None:
        actions = {
            "modified": lambda repository: (repository / "tracked.txt").write_text(
                "modified\n", encoding="utf-8"
            ),
            "staged": lambda repository: (
                (repository / "tracked.txt").write_text("staged\n", encoding="utf-8"),
                git(repository, "add", "tracked.txt"),
            ),
            "deleted": lambda repository: (repository / "tracked.txt").unlink(),
            "untracked": lambda repository: (repository / "new.txt").write_text(
                "new\n", encoding="utf-8"
            ),
        }
        for name, action in actions.items():
            with self.subTest(change=name):
                with tempfile.TemporaryDirectory(
                    prefix=f"ads-provenance-{name}-"
                ) as temporary:
                    repository = initialize_repository(Path(temporary))
                    action(repository)
                    self.assertTrue(inspect_repository(repository).dirty)

    def test_non_repository_fails_explicitly(self) -> None:
        with tempfile.TemporaryDirectory(prefix="ads-provenance-missing-") as temporary:
            with self.assertRaisesRegex(ProvenanceError, "cannot inspect"):
                inspect_repository(Path(temporary))


if __name__ == "__main__":
    unittest.main(verbosity=2)
