"""Safe detached-worktree ownership and deterministic A/B scheduling."""

from __future__ import annotations

from dataclasses import dataclass
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import tempfile
import uuid

from .framework.errors import ExecutionError, ValidationError
from .framework.storage import validate_identifier


WORKSPACE_SCHEMA_VERSION = 1
WORKSPACE_KIND = "ads-benchmark-ab-workspace"
_COMMIT = re.compile(r"^[0-9a-f]{40}(?:[0-9a-f]{24})?$")


def _git(repository: Path, *arguments: str) -> str:
    """Run one argument-vector Git query with bounded diagnostic output."""

    try:
        completed = subprocess.run(
            ("git", "-C", str(repository), *arguments),
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=60,
        )
    except (OSError, subprocess.SubprocessError) as error:
        raise ExecutionError(f"cannot run git {' '.join(arguments)}: {error}") from error
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip()
        if len(detail) > 2000:
            detail = detail[:2000] + "..."
        raise ExecutionError(
            f"git {' '.join(arguments)} failed with status "
            f"{completed.returncode}: {detail}"
        )
    return completed.stdout.strip()


def resolve_commit(repository: Path, reference: str) -> str:
    """Resolve a user ref once and return only its immutable full commit SHA."""

    if (
        not isinstance(reference, str)
        or not reference
        or reference.startswith("-")
        or "\0" in reference
        or "\n" in reference
        or "\r" in reference
    ):
        raise ValidationError(f"invalid Git reference: {reference!r}")
    commit = _git(
        repository.resolve(),
        "rev-parse",
        "--verify",
        "--end-of-options",
        f"{reference}^{{commit}}",
    ).lower()
    if not _COMMIT.fullmatch(commit):
        raise ExecutionError(
            f"Git returned an invalid full commit ID for {reference!r}: {commit!r}"
        )
    return commit


@dataclass(frozen=True)
class ScheduleEntry:
    case_id: str
    order: tuple[str, str]

    def to_dict(self) -> dict[str, object]:
        return {"case_id": self.case_id, "order": list(self.order)}


def alternating_schedule(case_ids: tuple[str, ...]) -> tuple[ScheduleEntry, ...]:
    """Alternate the first side for each matched case: AB, BA, AB, ..."""

    entries: list[ScheduleEntry] = []
    observed: set[str] = set()
    for index, case_id in enumerate(case_ids):
        validate_identifier(case_id, "case_id")
        if case_id in observed:
            raise ValidationError(f"duplicate case ID in A/B schedule: {case_id}")
        observed.add(case_id)
        order = ("baseline", "candidate") if index % 2 == 0 else (
            "candidate",
            "baseline",
        )
        entries.append(ScheduleEntry(case_id=case_id, order=order))
    if not entries:
        raise ValidationError("A/B schedule requires at least one case")
    return tuple(entries)


@dataclass(frozen=True)
class WorktreePair:
    workspace: Path
    baseline: Path
    candidate: Path
    baseline_commit: str
    candidate_commit: str
    marker: Path


class OwnedWorktreeWorkspace:
    """Create and remove only marker-owned detached A/B worktrees.

    Failure intentionally retains the workspace for diagnosis.  Cleanup never
    uses ``--force``, ``git clean``, pruning, or a recursive filesystem delete.
    """

    def __init__(
        self,
        repository_root: Path,
        *,
        run_id: str,
        baseline_ref: str,
        candidate_ref: str,
        parent: Path | None = None,
    ) -> None:
        self.repository_root = repository_root.resolve()
        self.run_id = validate_identifier(run_id, "run_id")
        self.baseline_ref = baseline_ref
        self.candidate_ref = candidate_ref
        self.parent = (parent or Path(tempfile.gettempdir())).resolve()
        self._pair: WorktreePair | None = None

    @property
    def pair(self) -> WorktreePair:
        if self._pair is None:
            raise ExecutionError("A/B worktree workspace has not been prepared")
        return self._pair

    def prepare(self) -> WorktreePair:
        if self._pair is not None:
            raise ExecutionError("A/B worktree workspace is already prepared")
        if self.parent.is_symlink() or not self.parent.is_dir():
            raise ExecutionError(f"workspace parent is missing or unsafe: {self.parent}")
        try:
            self.parent.relative_to(self.repository_root)
        except ValueError:
            pass
        else:
            raise ExecutionError(
                "workspace parent must be outside the benchmark repository: "
                f"{self.parent}"
            )

        baseline_commit = resolve_commit(self.repository_root, self.baseline_ref)
        candidate_commit = resolve_commit(self.repository_root, self.candidate_ref)
        common_git_dir = Path(
            _git(
                self.repository_root,
                "rev-parse",
                "--path-format=absolute",
                "--git-common-dir",
            )
        ).resolve()
        workspace = Path(
            tempfile.mkdtemp(prefix=f"ads-ab-{self.run_id}-", dir=self.parent)
        )
        marker = workspace / ".ads-benchmark-ab-workspace.json"
        nonce = uuid.uuid4().hex
        document = {
            "schema_version": WORKSPACE_SCHEMA_VERSION,
            "kind": WORKSPACE_KIND,
            "workspace": str(workspace.resolve()),
            "common_git_dir": str(common_git_dir),
            "run_id": self.run_id,
            "baseline_commit": baseline_commit,
            "candidate_commit": candidate_commit,
            "nonce": nonce,
        }
        descriptor: int | None = None
        try:
            descriptor = os.open(
                marker,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
                0o600,
            )
            with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
                descriptor = None
                json.dump(document, stream, indent=2, sort_keys=True, allow_nan=False)
                stream.write("\n")
                stream.flush()
                os.fsync(stream.fileno())
            workspace_fd = os.open(
                workspace, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
            )
            try:
                os.fsync(workspace_fd)
            finally:
                os.close(workspace_fd)
        except Exception as error:
            if descriptor is not None:
                os.close(descriptor)
            raise ExecutionError(
                f"cannot initialize owned A/B workspace {workspace}: {error}; "
                "workspace retained"
            ) from error

        baseline = workspace / "baseline"
        candidate = workspace / "candidate"
        try:
            _git(
                self.repository_root,
                "worktree",
                "add",
                "--detach",
                str(baseline),
                baseline_commit,
            )
            _git(
                self.repository_root,
                "worktree",
                "add",
                "--detach",
                str(candidate),
                candidate_commit,
            )
        except Exception as error:
            raise ExecutionError(
                f"cannot prepare A/B worktrees: {error}; diagnostic workspace "
                f"retained at {workspace}"
            ) from error

        self._pair = WorktreePair(
            workspace=workspace,
            baseline=baseline,
            candidate=candidate,
            baseline_commit=baseline_commit,
            candidate_commit=candidate_commit,
            marker=marker,
        )
        return self._pair

    def _validated_marker(self) -> dict[str, object]:
        pair = self.pair
        try:
            metadata = pair.marker.lstat()
            if not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1:
                raise ValueError("marker is not a singly linked regular file")
            document = json.loads(pair.marker.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as error:
            raise ExecutionError(
                f"refusing cleanup of unverified A/B workspace {pair.workspace}: {error}"
            ) from error
        expected = {
            "schema_version": WORKSPACE_SCHEMA_VERSION,
            "kind": WORKSPACE_KIND,
            "workspace": str(pair.workspace.resolve()),
            "common_git_dir": str(
                Path(
                    _git(
                        self.repository_root,
                        "rev-parse",
                        "--path-format=absolute",
                        "--git-common-dir",
                    )
                ).resolve()
            ),
            "run_id": self.run_id,
            "baseline_commit": pair.baseline_commit,
            "candidate_commit": pair.candidate_commit,
        }
        if not isinstance(document, dict) or set(document) != set(expected) | {"nonce"}:
            raise ExecutionError(
                f"refusing cleanup of A/B workspace with invalid marker: {pair.workspace}"
            )
        if any(document.get(key) != value for key, value in expected.items()):
            raise ExecutionError(
                f"refusing cleanup of A/B workspace with mismatched marker: {pair.workspace}"
            )
        nonce = document.get("nonce")
        if not isinstance(nonce, str) or not re.fullmatch(r"[0-9a-f]{32}", nonce):
            raise ExecutionError(
                f"refusing cleanup of A/B workspace with invalid nonce: {pair.workspace}"
            )
        return document

    def cleanup(self) -> None:
        """Remove clean owned worktrees and then the now-empty owned directory."""

        pair = self.pair
        self._validated_marker()
        for path in (pair.baseline, pair.candidate):
            if path.is_symlink() or not path.is_dir():
                raise ExecutionError(
                    f"refusing cleanup because owned worktree is missing or unsafe: {path}"
                )
            if _git(path, "status", "--porcelain=v1", "--untracked-files=normal"):
                raise ExecutionError(
                    f"refusing cleanup of modified A/B worktree; retained at {path}"
                )

        # All preflight checks happen before the first mutation.  Git removes
        # only the two paths it registered; no global prune or force is used.
        for path in (pair.baseline, pair.candidate):
            _git(self.repository_root, "worktree", "remove", str(path))

        remaining = {entry.name for entry in pair.workspace.iterdir()}
        if remaining != {pair.marker.name}:
            raise ExecutionError(
                f"refusing to remove nonempty A/B workspace {pair.workspace}; "
                f"unexpected entries: {sorted(remaining)}"
            )
        pair.marker.unlink()
        pair.workspace.rmdir()
        self._pair = None

    def __enter__(self) -> WorktreePair:
        return self.prepare()

    def __exit__(self, exception_type, exception, traceback) -> bool:
        if exception_type is None:
            self.cleanup()
        # Any failed orchestration intentionally leaves both worktrees and the
        # ownership marker intact for inspection and an explicit later retry.
        return False
