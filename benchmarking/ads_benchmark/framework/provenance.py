"""Git provenance captured for every written benchmark plan."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
import stat
import subprocess
from typing import Any

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


def _git_bytes(repository_root: Path, *arguments: str) -> bytes:
    command = ["git", "-C", os.fsencode(repository_root), *arguments]
    try:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ProvenanceError(f"cannot execute git provenance command: {error}") from error
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or b"git failed").decode(
            "utf-8", errors="replace"
        ).strip()
        raise ProvenanceError(f"cannot inspect repository provenance: {detail}")
    return completed.stdout


def _record(digest: Any, label: bytes, value: bytes) -> None:
    # Length framing makes concatenated paths and contents unambiguous.
    digest.update(len(label).to_bytes(8, "big"))
    digest.update(label)
    digest.update(len(value).to_bytes(8, "big"))
    digest.update(value)


def _snapshot_once(root: Path) -> tuple[bool, str]:
    status = _git_bytes(
        root,
        "status",
        "--porcelain=v1",
        "-z",
        "--untracked-files=all",
    )
    listed = _git_bytes(
        root,
        "ls-files",
        "-z",
        "--cached",
        "--others",
        "--exclude-standard",
    )
    paths = sorted({entry for entry in listed.split(b"\0") if entry})
    digest = hashlib.sha256()
    _record(digest, b"format", b"ads-worktree-v1")
    _record(digest, b"git-status", status)

    nofollow = getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_CLOEXEC", 0)
    for raw_path in paths:
        path = root / os.fsdecode(raw_path)
        _record(digest, b"path", raw_path)
        try:
            metadata = path.lstat()
        except FileNotFoundError:
            _record(digest, b"kind", b"missing")
            continue
        except OSError as error:
            raise ProvenanceError(f"cannot inspect worktree path {path}: {error}") from error

        mode = metadata.st_mode
        _record(digest, b"mode", str(mode & 0o177777).encode("ascii"))
        if stat.S_ISLNK(mode):
            try:
                target = os.readlink(path)
            except OSError as error:
                raise ProvenanceError(
                    f"cannot read worktree symlink {path}: {error}"
                ) from error
            _record(digest, b"symlink", os.fsencode(target))
            continue
        if stat.S_ISREG(mode):
            descriptor: int | None = None
            try:
                descriptor = os.open(path, os.O_RDONLY | nofollow)
                opened = os.fstat(descriptor)
                if (
                    not stat.S_ISREG(opened.st_mode)
                    or (metadata.st_dev, metadata.st_ino)
                    != (opened.st_dev, opened.st_ino)
                ):
                    raise ProvenanceError(
                        f"worktree path changed type while fingerprinting: {path}"
                    )
                _record(digest, b"size", str(opened.st_size).encode("ascii"))
                while True:
                    chunk = os.read(descriptor, 1024 * 1024)
                    if not chunk:
                        break
                    digest.update(chunk)
                finished = os.fstat(descriptor)
                if (
                    opened.st_dev,
                    opened.st_ino,
                    opened.st_size,
                    opened.st_mtime_ns,
                ) != (
                    finished.st_dev,
                    finished.st_ino,
                    finished.st_size,
                    finished.st_mtime_ns,
                ):
                    raise ProvenanceError(
                        f"worktree path changed while fingerprinting: {path}"
                    )
            except ProvenanceError:
                raise
            except OSError as error:
                raise ProvenanceError(
                    f"cannot fingerprint worktree file {path}: {error}"
                ) from error
            finally:
                if descriptor is not None:
                    os.close(descriptor)
            continue
        if stat.S_ISDIR(mode):
            # Gitlinks appear as directories in the superproject.  Include the
            # nested commit and nonignored worktree fingerprint when available.
            try:
                nested_commit = _git(path, "rev-parse", "--verify", "HEAD").strip()
                nested_dirty, nested_fingerprint = _stable_snapshot(path)
            except ProvenanceError as error:
                raise ProvenanceError(
                    f"cannot fingerprint worktree directory {path}: {error}"
                ) from error
            _record(digest, b"gitlink-commit", nested_commit.encode("ascii"))
            _record(
                digest,
                b"gitlink-worktree",
                f"{int(nested_dirty)}:{nested_fingerprint}".encode("ascii"),
            )
            continue
        raise ProvenanceError(f"unsupported special worktree path: {path}")
    return bool(status), f"sha256:{digest.hexdigest()}"


def _stable_snapshot(root: Path) -> tuple[bool, str]:
    previous: tuple[bool, str] | None = None
    for _ in range(3):
        current = _snapshot_once(root)
        if current == previous:
            return current
        previous = current
    raise ProvenanceError("worktree changed repeatedly while computing fingerprint")


def inspect_repository(repository_root: Path) -> RepositoryState:
    root = repository_root.resolve()
    commit = _git(root, "rev-parse", "--verify", "HEAD").strip()
    if len(commit) != 40 or any(character not in "0123456789abcdef" for character in commit):
        raise ProvenanceError(f"git returned an invalid commit id: {commit!r}")
    dirty, fingerprint = _stable_snapshot(root)
    return RepositoryState(
        commit=commit,
        dirty=dirty,
        worktree_fingerprint=fingerprint,
    )
