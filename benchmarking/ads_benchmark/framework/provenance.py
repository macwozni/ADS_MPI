"""Git provenance captured for every written benchmark plan."""

from __future__ import annotations

from pathlib import Path
import subprocess

from .errors import ProvenanceError
from .model import RepositoryState


def _git(repository_root: Path, *arguments: str) -> str:
    command = ["git", "-C", str(repository_root), *arguments]
    try:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ProvenanceError(f"cannot execute git provenance command: {error}") from error
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip() or "git failed"
        raise ProvenanceError(f"cannot inspect repository provenance: {detail}")
    return completed.stdout


def inspect_repository(repository_root: Path) -> RepositoryState:
    root = repository_root.resolve()
    commit = _git(root, "rev-parse", "--verify", "HEAD").strip()
    if len(commit) != 40 or any(character not in "0123456789abcdef" for character in commit):
        raise ProvenanceError(f"git returned an invalid commit id: {commit!r}")
    status = _git(root, "status", "--porcelain=v1", "--untracked-files=normal")
    return RepositoryState(commit=commit, dirty=bool(status.strip()))
