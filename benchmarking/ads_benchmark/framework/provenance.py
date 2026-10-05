"""Bounded repository, build, launcher, and host provenance collection.

The immutable plan manifest intentionally remains at schema version 1.  An
executed run gets a separate ``execution.json`` record whose hashes bind it to
the exact manifest while keeping machine observations out of semantic case
identity and ``config_hash``.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from copy import deepcopy
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import socket
import stat
import subprocess
import tempfile
from typing import Any, Callable

from .errors import ProvenanceError
from .model import RepositoryState


EXECUTION_SCHEMA_VERSION = 1
EXECUTION_KIND = "ads-benchmark-execution"

_SHA256 = re.compile(r"^sha256:[0-9a-f]{64}$")
_MUMPS_VERSION = re.compile(
    r'^\s*#\s*define\s+MUMPS_VERSION\s+"([^"\r\n]+)"', re.MULTILINE
)
_OMP_NAME = re.compile(r"^(?:OMP|GOMP|KMP)_[A-Z0-9_]+$")
_MAX_HASH_BYTES = 512 * 1024 * 1024
_MAX_STAMP_BYTES = 1024 * 1024
_MAX_SYSTEM_TEXT_BYTES = 4 * 1024 * 1024
_MAX_PROBE_OUTPUT_BYTES = 64 * 1024
_PROBE_TIMEOUT_SECONDS = 5.0
_MAX_BUILD_UNITS = 128
_MAX_TRACKED_LIBRARIES = 256
_MAX_ENVIRONMENT_VARIABLES = 256
_MAX_LAUNCHERS = 64
_MAX_LAUNCHER_ARGUMENTS = 1024
_MAX_TEXT_LENGTH = 1024 * 1024

_IMPORTANT_LIBRARY_KEYS = (
    "mumps",
    "blas",
    "lapack",
    "scalapack",
    "parmetis",
    "metis",
    "gklib",
    "other",
)
_LIBRARY_DIRECTORY_ENVIRONMENT = {
    "MUMPS_DIR": "mumps",
    "BLAS_DIR": "blas",
    "LAPACK_DIR": "lapack",
    "SCALAPACK_DIR": "scalapack",
    "PARMETIS_DIR": "parmetis",
    "METIS_DIR": "metis",
    "GKLIB_DIR": "gklib",
}
_IMPORTANT_OMP_ENVIRONMENT = frozenset(
    {
        "OMP_NUM_THREADS",
        "OMP_DYNAMIC",
        "OMP_PROC_BIND",
        "OMP_PLACES",
        "OMP_WAIT_POLICY",
        "OMP_STACKSIZE",
        "OMP_THREAD_LIMIT",
        "OMP_MAX_ACTIVE_LEVELS",
        "OMP_DISPLAY_ENV",
        "GOMP_CPU_AFFINITY",
        "GOMP_STACKSIZE",
        "KMP_AFFINITY",
        "KMP_BLOCKTIME",
        "KMP_SETTINGS",
        "KMP_HW_SUBSET",
    }
)

_EXECUTION_KEYS = {
    "schema_version",
    "kind",
    "run_id",
    "manifest_hash",
    "config_hash",
    "recorded_at",
    "timezone",
    "system",
    "topology",
    "environment",
    "build",
    "launchers",
    "build_identity_hash",
    "execution_compatibility_hash",
    "record_hash",
}
_TIMEZONE_KEYS = {"name", "utc_offset_seconds"}
_SYSTEM_KEYS = {
    "hostname",
    "operating_system",
    "kernel",
    "architecture",
    "cpu_model",
    "physical_cores",
    "logical_cores",
    "memory_bytes",
}
_TOPOLOGY_KEYS = {"rank_grids", "bindings"}
_RANK_GRID_KEYS = {"ranks", "process_grid"}
_BINDING_KEYS = {"threads", "dynamic", "proc_bind", "places"}
_BUILD_KEYS = {
    "profiles",
    "units",
    "compiler_versions",
    "libraries",
    "configured_libraries",
    "mumps_version",
}
_BUILD_UNIT_KEYS = {
    "profile",
    "problem",
    "executable",
    "core_stamp",
    "adapter_stamp",
    "compiler_command",
    "compile_flags",
    "link_flags",
}
_FILE_KEYS = {"path", "size_bytes", "sha256", "reason"}
_COMPILER_VERSION_KEYS = {"command", "resolved_executable", "version"}
_LAUNCHER_KEYS = {
    "name",
    "kind",
    "executable",
    "resolved_executable",
    "arguments",
    "version",
    "available_mpi_slots",
    "available_cpu_slots",
    "allocated_threads_per_rank",
}
_OBSERVATION_KEYS = {"value", "reason"}


def _canonical_json(value: object) -> str:
    try:
        return json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
    except (TypeError, ValueError, RecursionError) as error:
        raise ProvenanceError(f"provenance data is not canonical JSON: {error}") from error


def _digest(value: object) -> str:
    return "sha256:" + hashlib.sha256(
        _canonical_json(value).encode("utf-8")
    ).hexdigest()


def manifest_digest(manifest: Mapping[str, object]) -> str:
    """Return the stable digest used to bind execution metadata to a plan."""

    if not isinstance(manifest, Mapping):
        raise ProvenanceError("manifest must be an object")
    return _digest(dict(manifest))


def _available(value: object) -> dict[str, object]:
    return {"value": value, "reason": None}


def _unavailable(reason: str) -> dict[str, object]:
    if not reason:
        reason = "unavailable"
    return {"value": None, "reason": reason[:4096]}


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
    if len(commit) not in {40, 64} or any(
        character not in "0123456789abcdef" for character in commit
    ):
        raise ProvenanceError(f"git returned an invalid commit id: {commit!r}")
    dirty, fingerprint = _stable_snapshot(root)
    return RepositoryState(
        commit=commit,
        dirty=dirty,
        worktree_fingerprint=fingerprint,
    )


def _read_bounded_text(path: Path, maximum: int) -> tuple[str | None, str | None]:
    descriptor: int | None = None
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_CLOEXEC", 0)
    try:
        descriptor = os.open(path, flags)
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            return None, "not a regular file"
        if metadata.st_size > maximum:
            return None, f"file exceeds {maximum} byte provenance limit"
        chunks: list[bytes] = []
        remaining = maximum + 1
        while remaining > 0:
            chunk = os.read(descriptor, min(1024 * 1024, remaining))
            if not chunk:
                break
            chunks.append(chunk)
            remaining -= len(chunk)
        content = b"".join(chunks)
        if len(content) > maximum:
            return None, f"file exceeds {maximum} byte provenance limit"
        return content.decode("utf-8"), None
    except FileNotFoundError:
        return None, "file does not exist"
    except (OSError, UnicodeError) as error:
        return None, f"cannot read file: {error}"
    finally:
        if descriptor is not None:
            os.close(descriptor)


def _file_identity(path: Path) -> dict[str, object]:
    absolute = path if path.is_absolute() else path.absolute()
    descriptor: int | None = None
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_CLOEXEC", 0)
    try:
        descriptor = os.open(absolute, flags)
        opened = os.fstat(descriptor)
        if not stat.S_ISREG(opened.st_mode):
            return {
                "path": str(absolute),
                "size_bytes": None,
                "sha256": None,
                "reason": "not a regular file",
            }
        if opened.st_size > _MAX_HASH_BYTES:
            return {
                "path": str(absolute),
                "size_bytes": opened.st_size,
                "sha256": None,
                "reason": f"file exceeds {_MAX_HASH_BYTES} byte hash limit",
            }
        digest = hashlib.sha256()
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
            return {
                "path": str(absolute),
                "size_bytes": None,
                "sha256": None,
                "reason": "file changed while hashing",
            }
        return {
            "path": str(absolute),
            "size_bytes": opened.st_size,
            "sha256": f"sha256:{digest.hexdigest()}",
            "reason": None,
        }
    except FileNotFoundError:
        return {
            "path": str(absolute),
            "size_bytes": None,
            "sha256": None,
            "reason": "file does not exist",
        }
    except OSError as error:
        return {
            "path": str(absolute),
            "size_bytes": None,
            "sha256": None,
            "reason": f"cannot hash file: {error}",
        }
    finally:
        if descriptor is not None:
            os.close(descriptor)


def _parse_stamp(path: Path) -> tuple[dict[str, str] | None, str | None]:
    content, error = _read_bounded_text(path, _MAX_STAMP_BYTES)
    if content is None:
        return None, error
    values: dict[str, str] = {}
    for line_number, raw_line in enumerate(content.splitlines(), start=1):
        if not raw_line:
            continue
        if "=" not in raw_line:
            return None, f"stamp line {line_number} has no '=' separator"
        key, value = raw_line.split("=", 1)
        if not key or key in values:
            return None, f"stamp line {line_number} has invalid key {key!r}"
        values[key] = value
    return values, None


def _resolved_executable(command: str, environment: Mapping[str, str]) -> dict[str, object]:
    try:
        argv = shlex.split(command)
    except ValueError as error:
        return _unavailable(f"cannot parse command: {error}")
    if not argv:
        return _unavailable("command is empty")
    token = argv[0]
    if os.sep in token or (os.altsep and os.altsep in token):
        candidate = Path(token)
        if not candidate.is_absolute():
            candidate = Path.cwd() / candidate
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return _available(str(candidate.absolute()))
        return _unavailable(f"executable does not exist or is not executable: {candidate}")
    located = shutil.which(token, path=environment.get("PATH"))
    if located is None:
        return _unavailable(f"executable is not available on PATH: {token}")
    return _available(located)


def _probe_version(resolved: Mapping[str, object]) -> dict[str, object]:
    executable = resolved.get("value")
    if not isinstance(executable, str):
        reason = resolved.get("reason")
        return _unavailable(
            f"version probe skipped because executable is unavailable: {reason}"
        )
    try:
        with tempfile.TemporaryFile() as output:
            completed = subprocess.run(
                [executable, "--version"],
                check=False,
                stdin=subprocess.DEVNULL,
                stdout=output,
                stderr=subprocess.STDOUT,
                timeout=_PROBE_TIMEOUT_SECONDS,
            )
            output.seek(0, os.SEEK_END)
            size = output.tell()
            if size > _MAX_PROBE_OUTPUT_BYTES:
                return _unavailable(
                    f"version output exceeds {_MAX_PROBE_OUTPUT_BYTES} byte limit"
                )
            output.seek(0)
            text = output.read(_MAX_PROBE_OUTPUT_BYTES + 1).decode(
                "utf-8", errors="replace"
            ).strip()
    except subprocess.TimeoutExpired:
        return _unavailable(
            f"version probe exceeded {_PROBE_TIMEOUT_SECONDS:g} second timeout"
        )
    except OSError as error:
        return _unavailable(f"version probe failed: {error}")
    if completed.returncode != 0:
        return _unavailable(
            f"version probe exited with status {completed.returncode}"
        )
    if not text:
        return _unavailable("version probe produced no output")
    return _available(text)


def _library_group(path: str) -> str:
    name = Path(path).name.lower()
    if "mumps" in name or name == "libpord.a":
        return "mumps"
    if "scalapack" in name:
        return "scalapack"
    if "lapack" in name:
        return "lapack"
    if "blas" in name:
        return "blas"
    if "parmetis" in name:
        return "parmetis"
    if "metis" in name:
        return "metis"
    if "gklib" in name:
        return "gklib"
    return "other"


def _split_link_inputs(value: str | None) -> tuple[str, ...]:
    if not value:
        return ()
    try:
        tokens = shlex.split(value)
    except ValueError:
        return ()
    return tuple(
        token
        for token in tokens
        if token and not token.startswith("-")
    )


def _manifest_matrix(
    manifest: Mapping[str, object],
) -> tuple[str, str, list[tuple[str, str]], list[dict[str, object]], list[dict[str, object]]]:
    run_id = manifest.get("run_id")
    config_hash = manifest.get("config_hash")
    raw_cases = manifest.get("cases")
    if not isinstance(run_id, str) or not run_id:
        raise ProvenanceError("manifest run_id must be a nonempty string")
    if not isinstance(config_hash, str) or _SHA256.fullmatch(config_hash) is None:
        raise ProvenanceError("manifest config_hash must be a SHA-256 digest")
    if not isinstance(raw_cases, list) or not raw_cases:
        raise ProvenanceError("manifest cases must be a nonempty array")
    units: set[tuple[str, str]] = set()
    rank_grids: set[tuple[int, int, int, int]] = set()
    bindings: set[tuple[int, bool, str, str]] = set()
    for index, entry in enumerate(raw_cases):
        if not isinstance(entry, dict) or not isinstance(entry.get("configuration"), dict):
            raise ProvenanceError(f"manifest cases[{index}] has no configuration")
        configuration = entry["configuration"]
        problem = configuration.get("problem")
        build = configuration.get("build")
        mpi = configuration.get("mpi")
        openmp = configuration.get("openmp")
        if not isinstance(problem, str) or not problem:
            raise ProvenanceError(f"manifest cases[{index}].problem is invalid")
        if not isinstance(build, dict) or not isinstance(build.get("profile"), str):
            raise ProvenanceError(f"manifest cases[{index}].build is invalid")
        profile = build["profile"]
        assert isinstance(profile, str)
        units.add((profile, problem))
        if not isinstance(mpi, dict):
            raise ProvenanceError(f"manifest cases[{index}].mpi is invalid")
        ranks = mpi.get("ranks")
        grid = mpi.get("process_grid")
        if (
            type(ranks) is not int
            or ranks <= 0
            or not isinstance(grid, list)
            or len(grid) != 3
            or any(type(item) is not int or item <= 0 for item in grid)
        ):
            raise ProvenanceError(f"manifest cases[{index}].mpi is invalid")
        rank_grids.add((ranks, grid[0], grid[1], grid[2]))
        if not isinstance(openmp, dict):
            raise ProvenanceError(f"manifest cases[{index}].openmp is invalid")
        threads = openmp.get("threads")
        if type(threads) is not int or threads <= 0:
            raise ProvenanceError(f"manifest cases[{index}].openmp.threads is invalid")
        dynamic = openmp.get("dynamic", False)
        proc_bind = openmp.get("proc_bind", "close")
        places = openmp.get("places", "cores")
        if (
            type(dynamic) is not bool
            or not isinstance(proc_bind, str)
            or not isinstance(places, str)
        ):
            raise ProvenanceError(f"manifest cases[{index}].openmp binding is invalid")
        bindings.add((threads, dynamic, proc_bind, places))
    if len(units) > _MAX_BUILD_UNITS:
        raise ProvenanceError(
            f"manifest requires more than {_MAX_BUILD_UNITS} build units"
        )
    binding_documents = [
        {
            "threads": threads,
            "dynamic": dynamic,
            "proc_bind": proc_bind,
            "places": places,
        }
        for threads, dynamic, proc_bind, places in bindings
    ]
    binding_documents.sort(key=_canonical_json)
    return (
        run_id,
        config_hash,
        sorted(units),
        [
            {"ranks": ranks, "process_grid": [x, y, z]}
            for ranks, x, y, z in sorted(rank_grids)
        ],
        binding_documents,
    )


def _collect_build(
    root: Path,
    units: Sequence[tuple[str, str]],
    environment: Mapping[str, str],
) -> dict[str, object]:
    unit_documents: list[dict[str, object]] = []
    compiler_commands: set[str] = set()
    link_inputs: set[str] = set()
    configured: dict[str, set[str]] = {
        name: set() for name in _IMPORTANT_LIBRARY_KEYS
    }
    mumps_header_candidates: set[Path] = set()

    for variable, group in _LIBRARY_DIRECTORY_ENVIRONMENT.items():
        value = environment.get(variable)
        if value:
            configured[group].add(value)
            if variable == "MUMPS_DIR":
                mumps_header_candidates.add(Path(value) / "include" / "dmumps_c.h")
    for token in _split_link_inputs(environment.get("USER_LIB")):
        link_inputs.add(token)
        configured[_library_group(token)].add(token)

    for profile, problem in units:
        executable = (
            root
            / "benchmarking"
            / "build"
            / profile
            / "EXEC"
            / f"{problem}_manufactured"
        )
        core_stamp = (
            root
            / "benchmarking"
            / "build"
            / profile
            / "_OBJ"
            / ".ads-core-build"
        )
        adapter_stamp = (
            root
            / "benchmarking"
            / "build"
            / profile
            / "benchmark_OBJ"
            / problem
            / ".ads-benchmark-adapter-build"
        )
        stamp, stamp_error = _parse_stamp(adapter_stamp)
        core_values, core_error = _parse_stamp(core_stamp)
        compiler_command = None
        compile_flags = None
        link_flags = None
        if stamp is not None:
            compiler_command = stamp.get("FF")
            compile_flags = stamp.get("FFLAGS")
            link_flags = stamp.get("USER_LIB")
        if compiler_command is None and core_values is not None:
            compiler_command = core_values.get("FF")
        if compile_flags is None and core_values is not None:
            compile_flags = core_values.get("FFLAGS")
        if compiler_command is None:
            compiler_command = (
                environment.get("FF")
                or environment.get("MPIFC")
                or environment.get("FC")
            )
        if compile_flags is None:
            compile_flags = environment.get("FFLAGS")
        if link_flags is None:
            link_flags = environment.get("USER_LIB")
        if compiler_command:
            compiler_commands.add(compiler_command)
        for token in _split_link_inputs(link_flags):
            link_inputs.add(token)
            group = _library_group(token)
            configured[group].add(token)
            path = Path(token)
            if group == "mumps":
                try:
                    prefix = path.parent.parent
                    mumps_header_candidates.add(prefix / "include" / "dmumps_c.h")
                except IndexError:
                    pass
        for flag_text in (compiler_command, compile_flags):
            if not flag_text:
                continue
            try:
                flag_tokens = shlex.split(flag_text)
            except ValueError:
                continue
            for position, token in enumerate(flag_tokens):
                include = None
                if token == "-I" and position + 1 < len(flag_tokens):
                    include = flag_tokens[position + 1]
                elif token.startswith("-I") and len(token) > 2:
                    include = token[2:]
                if include and "mumps" in include.lower():
                    mumps_header_candidates.add(Path(include) / "dmumps_c.h")

        unit_documents.append(
            {
                "profile": profile,
                "problem": problem,
                "executable": _file_identity(executable),
                "core_stamp": _file_identity(core_stamp),
                "adapter_stamp": _file_identity(adapter_stamp),
                "compiler_command": (
                    _available(compiler_command)
                    if compiler_command
                    else _unavailable(
                        stamp_error
                        or core_error
                        or "compiler command is not present in build stamps or environment"
                    )
                ),
                "compile_flags": (
                    _available(compile_flags)
                    if compile_flags is not None
                    else _unavailable(
                        stamp_error
                        or core_error
                        or "compile flags are not present in build stamps or environment"
                    )
                ),
                "link_flags": (
                    _available(link_flags)
                    if link_flags is not None
                    else _unavailable(
                        stamp_error
                        or "link flags are not present in adapter stamp or environment"
                    )
                ),
            }
        )

    compiler_versions: list[dict[str, object]] = []
    for command in sorted(compiler_commands):
        resolved = _resolved_executable(command, environment)
        compiler_versions.append(
            {
                "command": command,
                "resolved_executable": resolved,
                "version": _probe_version(resolved),
            }
        )

    if len(link_inputs) > _MAX_TRACKED_LIBRARIES:
        raise ProvenanceError(
            f"build references more than {_MAX_TRACKED_LIBRARIES} link libraries"
        )
    libraries = [_file_identity(Path(path)) for path in sorted(link_inputs)]
    libraries.sort(key=lambda item: str(item["path"]))

    mumps_version: dict[str, object] = _unavailable(
        "no readable MUMPS version header was discovered"
    )
    for header in sorted(mumps_header_candidates, key=str):
        content, error = _read_bounded_text(header, _MAX_STAMP_BYTES)
        if content is None:
            continue
        match = _MUMPS_VERSION.search(content)
        if match is not None:
            mumps_version = _available(
                {"version": match.group(1), "header": str(header.absolute())}
            )
            break
        if error is not None:
            mumps_version = _unavailable(error)

    configured_document = {
        name: (
            _available(sorted(values))
            if values
            else _unavailable(f"no configured {name} path was discovered")
        )
        for name, values in configured.items()
    }
    return {
        "profiles": sorted({profile for profile, _ in units}),
        "units": unit_documents,
        "compiler_versions": (
            _available(compiler_versions)
            if compiler_versions
            else _unavailable(
                "no compiler command was present in build stamps or environment"
            )
        ),
        "libraries": libraries,
        "configured_libraries": configured_document,
        "mumps_version": mumps_version,
    }


def _read_system_text(path: Path) -> str | None:
    content, _ = _read_bounded_text(path, _MAX_SYSTEM_TEXT_BYTES)
    return content


def _cpu_details() -> tuple[dict[str, object], dict[str, object]]:
    content = _read_system_text(Path("/proc/cpuinfo"))
    if content is None:
        processor = platform.processor().strip()
        return (
            _available(processor)
            if processor
            else _unavailable("CPU model is unavailable"),
            _unavailable("physical core topology is unavailable"),
        )
    model = None
    physical_pairs: set[tuple[str, str]] = set()
    for block in content.split("\n\n"):
        fields: dict[str, str] = {}
        for line in block.splitlines():
            if ":" in line:
                key, value = line.split(":", 1)
                fields[key.strip().lower()] = value.strip()
        if model is None:
            model = fields.get("model name") or fields.get("hardware")
        physical = fields.get("physical id")
        core = fields.get("core id")
        if physical is not None and core is not None:
            physical_pairs.add((physical, core))
    model_observation = (
        _available(model) if model else _unavailable("CPU model is unavailable")
    )
    physical_observation = (
        _available(len(physical_pairs))
        if physical_pairs
        else _unavailable("physical core topology is unavailable")
    )
    return model_observation, physical_observation


def _memory_bytes() -> dict[str, object]:
    try:
        page_size = os.sysconf("SC_PAGE_SIZE")
        pages = os.sysconf("SC_PHYS_PAGES")
        if type(page_size) is int and type(pages) is int and page_size > 0 and pages > 0:
            return _available(page_size * pages)
    except (OSError, ValueError):
        pass
    content = _read_system_text(Path("/proc/meminfo"))
    if content is not None:
        match = re.search(r"^MemTotal:\s+([0-9]+)\s+kB\s*$", content, re.MULTILINE)
        if match is not None:
            return _available(int(match.group(1)) * 1024)
    return _unavailable("physical memory size is unavailable")


def _collect_system() -> dict[str, object]:
    cpu_model, physical_cores = _cpu_details()
    logical = os.cpu_count()
    return {
        "hostname": (
            _available(socket.gethostname())
            if socket.gethostname()
            else _unavailable("hostname is unavailable")
        ),
        "operating_system": (
            _available(platform.system())
            if platform.system()
            else _unavailable("operating system is unavailable")
        ),
        "kernel": (
            _available(platform.release())
            if platform.release()
            else _unavailable("kernel release is unavailable")
        ),
        "architecture": (
            _available(platform.machine())
            if platform.machine()
            else _unavailable("machine architecture is unavailable")
        ),
        "cpu_model": cpu_model,
        "physical_cores": physical_cores,
        "logical_cores": (
            _available(logical)
            if type(logical) is int and logical > 0
            else _unavailable("logical CPU count is unavailable")
        ),
        "memory_bytes": _memory_bytes(),
    }


def _collect_environment(environment: Mapping[str, str]) -> dict[str, object]:
    names = sorted(
        _IMPORTANT_OMP_ENVIRONMENT
        | {
            name
            for name in environment
            if isinstance(name, str) and _OMP_NAME.fullmatch(name)
        }
    )
    if len(names) > _MAX_ENVIRONMENT_VARIABLES:
        raise ProvenanceError(
            f"more than {_MAX_ENVIRONMENT_VARIABLES} OMP/GOMP/KMP variables are set"
        )
    document: dict[str, object] = {}
    for name in names:
        value = environment.get(name)
        if value is None:
            document[name] = _unavailable("environment variable is not set")
        elif not isinstance(value, str):
            document[name] = _unavailable("environment value is not a string")
        elif len(value) > 4096:
            document[name] = _unavailable("environment value exceeds 4096 byte limit")
        else:
            document[name] = _available(value)
    return document


def _collect_launchers(
    descriptions: Iterable[Mapping[str, object]],
    environment: Mapping[str, str],
) -> dict[str, object]:
    records: list[dict[str, object]] = []
    for position, description in enumerate(descriptions):
        if position >= _MAX_LAUNCHERS:
            raise ProvenanceError(
                f"more than {_MAX_LAUNCHERS} launcher descriptions were supplied"
            )
        if not isinstance(description, Mapping):
            raise ProvenanceError(f"launcher description {position} is not an object")
        name = description.get("name")
        kind = description.get("kind")
        executable = description.get("executable")
        arguments = description.get("arguments")
        if not isinstance(name, str) or not name:
            raise ProvenanceError(f"launcher description {position}.name is invalid")
        if not isinstance(kind, str) or not kind:
            raise ProvenanceError(f"launcher description {position}.kind is invalid")
        if executable is not None and not isinstance(executable, str):
            raise ProvenanceError(
                f"launcher description {position}.executable is invalid"
            )
        if (
            not isinstance(arguments, list)
            or len(arguments) > _MAX_LAUNCHER_ARGUMENTS
            or any(
                not isinstance(argument, str)
                or not argument
                or "\0" in argument
                or len(argument) > 4096
                for argument in arguments
            )
        ):
            raise ProvenanceError(
                f"launcher description {position}.arguments is invalid"
            )
        executable_observation = (
            _available(executable)
            if executable is not None
            else _unavailable("launcher has no separate executable")
        )
        resolved = (
            _resolved_executable(executable, environment)
            if executable is not None
            else _unavailable("launcher has no separate executable")
        )
        records.append(
            {
                "name": name,
                "kind": kind,
                "executable": executable_observation,
                "resolved_executable": resolved,
                "arguments": list(arguments),
                "version": _probe_version(resolved),
                "available_mpi_slots": _optional_capacity(
                    description.get("available_mpi_slots"),
                    "available MPI slots were not declared",
                ),
                "available_cpu_slots": _optional_capacity(
                    description.get("available_cpu_slots"),
                    "available CPU slots were not declared",
                ),
                "allocated_threads_per_rank": _optional_capacity(
                    description.get("allocated_threads_per_rank"),
                    "allocated threads per rank were not declared",
                ),
            }
        )
    if not records:
        return _unavailable("launcher description was not supplied")
    records.sort(key=_canonical_json)
    return _available(records)


def _optional_capacity(value: object, reason: str) -> dict[str, object]:
    if value is None:
        return _unavailable(reason)
    if type(value) is not int or value <= 0:
        raise ProvenanceError("launcher capacity must be a positive integer")
    return _available(value)


def collect_execution_provenance(
    repository_root: Path,
    manifest: Mapping[str, object],
    *,
    launcher_descriptions: Iterable[Mapping[str, object]] = (),
    environment: Mapping[str, str] | None = None,
    recorded_at: datetime | None = None,
) -> dict[str, object]:
    """Collect one bounded, self-validating execution record.

    Collection never guesses absent versions.  Every unavailable scalar is an
    explicit ``{"value": null, "reason": ...}`` observation.  Launcher
    descriptions are deliberately supplied by the composition layer so this
    module does not depend on the component registry or execute a benchmark.
    """

    root = repository_root.resolve()
    runtime_environment: Mapping[str, str] = (
        dict(os.environ) if environment is None else dict(environment)
    )
    run_id, config_hash, units, rank_grids, bindings = _manifest_matrix(manifest)
    instant = (recorded_at or datetime.now().astimezone()).astimezone()
    offset = instant.utcoffset()
    zone_name = instant.tzname()
    timezone_document = {
        "name": (
            _available(zone_name)
            if zone_name
            else _unavailable("local timezone name is unavailable")
        ),
        "utc_offset_seconds": (
            _available(int(offset.total_seconds()))
            if offset is not None
            else _unavailable("local UTC offset is unavailable")
        ),
    }
    build = _collect_build(root, units, runtime_environment)
    launchers = _collect_launchers(launcher_descriptions, runtime_environment)
    document: dict[str, object] = {
        "schema_version": EXECUTION_SCHEMA_VERSION,
        "kind": EXECUTION_KIND,
        "run_id": run_id,
        "manifest_hash": manifest_digest(manifest),
        "config_hash": config_hash,
        "recorded_at": instant.isoformat(),
        "timezone": timezone_document,
        "system": _collect_system(),
        "topology": {"rank_grids": rank_grids, "bindings": bindings},
        "environment": _collect_environment(runtime_environment),
        "build": build,
        "launchers": launchers,
    }
    document["build_identity_hash"] = _digest(
        {"format": "ads-build-identity-v1", "build": build}
    )
    document["execution_compatibility_hash"] = _digest(
        {
            "format": "ads-execution-compatibility-v1",
            "build_identity_hash": document["build_identity_hash"],
            "launchers": launchers,
        }
    )
    document["record_hash"] = _digest(
        {"format": "ads-execution-record-v1", "record": document}
    )
    validate_execution_provenance(document, manifest=manifest)
    return document


def _exact_object(value: object, keys: set[str], field: str) -> Mapping[str, object]:
    if not isinstance(value, dict):
        raise ProvenanceError(f"{field} must be an object")
    missing = sorted(keys - set(value))
    extra = sorted(set(value) - keys)
    if missing or extra:
        details: list[str] = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if extra:
            details.append("unknown " + ", ".join(extra))
        raise ProvenanceError(f"{field} has invalid keys: {'; '.join(details)}")
    return value


def _text(value: object, field: str, *, allow_empty: bool = False) -> str:
    if (
        not isinstance(value, str)
        or (not allow_empty and not value)
        or "\0" in value
        or len(value) > _MAX_TEXT_LENGTH
    ):
        raise ProvenanceError(f"{field} must be bounded NUL-free text")
    return value


def _json_value(value: object, field: str, depth: int = 0) -> None:
    if depth > 12:
        raise ProvenanceError(f"{field} exceeds maximum nesting depth")
    if value is None or type(value) in {bool, int, float}:
        if isinstance(value, float) and (value != value or value in {float("inf"), float("-inf")}):
            raise ProvenanceError(f"{field} contains a non-finite number")
        return
    if isinstance(value, str):
        _text(value, field, allow_empty=True)
        return
    if isinstance(value, list):
        if len(value) > 4096:
            raise ProvenanceError(f"{field} contains too many entries")
        for index, item in enumerate(value):
            _json_value(item, f"{field}[{index}]", depth + 1)
        return
    if isinstance(value, dict):
        if len(value) > 4096:
            raise ProvenanceError(f"{field} contains too many keys")
        for key, item in value.items():
            _text(key, f"{field} key")
            _json_value(item, f"{field}.{key}", depth + 1)
        return
    raise ProvenanceError(f"{field} contains a non-JSON value")


def _observation(
    value: object,
    field: str,
    validator: Callable[[object, str], None] | None = None,
) -> Mapping[str, object]:
    observation = _exact_object(value, _OBSERVATION_KEYS, field)
    observed = observation["value"]
    reason = observation["reason"]
    if observed is None:
        _text(reason, f"{field}.reason")
    else:
        if reason is not None:
            raise ProvenanceError(f"{field}.reason must be null when value is available")
        if validator is None:
            _json_value(observed, f"{field}.value")
        else:
            validator(observed, f"{field}.value")
    return observation


def _string_value(value: object, field: str) -> None:
    _text(value, field, allow_empty=True)


def _positive_integer_value(value: object, field: str) -> None:
    if type(value) is not int or value <= 0:
        raise ProvenanceError(f"{field} must be a positive integer")


def _integer_value(value: object, field: str) -> None:
    if type(value) is not int:
        raise ProvenanceError(f"{field} must be an integer")


def _file_record(value: object, field: str) -> Mapping[str, object]:
    record = _exact_object(value, _FILE_KEYS, field)
    _text(record["path"], f"{field}.path")
    reason = record["reason"]
    size = record["size_bytes"]
    digest = record["sha256"]
    if reason is None:
        if type(size) is not int or size < 0:
            raise ProvenanceError(f"{field}.size_bytes must be nonnegative")
        if not isinstance(digest, str) or _SHA256.fullmatch(digest) is None:
            raise ProvenanceError(f"{field}.sha256 must be a SHA-256 digest")
    else:
        _text(reason, f"{field}.reason")
        if size is not None and (type(size) is not int or size < 0):
            raise ProvenanceError(f"{field}.size_bytes must be null or nonnegative")
        if digest is not None:
            raise ProvenanceError(f"{field}.sha256 must be null when unavailable")
    return record


def _validate_build(value: object) -> Mapping[str, object]:
    build = _exact_object(value, _BUILD_KEYS, "execution.build")
    profiles = build["profiles"]
    if not isinstance(profiles, list) or not profiles:
        raise ProvenanceError("execution.build.profiles must be a nonempty array")
    for index, profile in enumerate(profiles):
        _text(profile, f"execution.build.profiles[{index}]")
    if profiles != sorted(set(profiles)):
        raise ProvenanceError("execution.build.profiles must be sorted and unique")

    units = build["units"]
    if not isinstance(units, list) or not units or len(units) > _MAX_BUILD_UNITS:
        raise ProvenanceError("execution.build.units has invalid cardinality")
    coordinates: list[tuple[str, str]] = []
    for index, raw_unit in enumerate(units):
        field = f"execution.build.units[{index}]"
        unit = _exact_object(raw_unit, _BUILD_UNIT_KEYS, field)
        profile = _text(unit["profile"], f"{field}.profile")
        problem = _text(unit["problem"], f"{field}.problem")
        coordinates.append((profile, problem))
        for name in ("executable", "core_stamp", "adapter_stamp"):
            _file_record(unit[name], f"{field}.{name}")
        for name in ("compiler_command", "compile_flags", "link_flags"):
            _observation(unit[name], f"{field}.{name}", _string_value)
    if coordinates != sorted(set(coordinates)):
        raise ProvenanceError("execution.build.units must be sorted and unique")

    versions_observation = _observation(
        build["compiler_versions"], "execution.build.compiler_versions"
    )
    versions = versions_observation["value"]
    if versions is None:
        versions = []
    if not isinstance(versions, list) or len(versions) > _MAX_BUILD_UNITS:
        raise ProvenanceError(
            "execution.build.compiler_versions.value must be an array"
        )
    version_commands: list[str] = []
    for index, raw_version in enumerate(versions):
        field = f"execution.build.compiler_versions[{index}]"
        version = _exact_object(raw_version, _COMPILER_VERSION_KEYS, field)
        version_commands.append(_text(version["command"], f"{field}.command"))
        _observation(version["resolved_executable"], f"{field}.resolved_executable", _string_value)
        _observation(version["version"], f"{field}.version", _string_value)
    if version_commands != sorted(set(version_commands)):
        raise ProvenanceError("compiler version commands must be sorted and unique")

    libraries = build["libraries"]
    if not isinstance(libraries, list) or len(libraries) > _MAX_TRACKED_LIBRARIES:
        raise ProvenanceError("execution.build.libraries must be an array")
    library_paths: list[str] = []
    for index, library in enumerate(libraries):
        record = _file_record(library, f"execution.build.libraries[{index}]")
        assert isinstance(record["path"], str)
        library_paths.append(record["path"])
    if library_paths != sorted(set(library_paths)):
        raise ProvenanceError("execution.build.libraries must be sorted and unique")

    configured = _exact_object(
        build["configured_libraries"],
        set(_IMPORTANT_LIBRARY_KEYS),
        "execution.build.configured_libraries",
    )
    for name in _IMPORTANT_LIBRARY_KEYS:
        observation = _observation(
            configured[name], f"execution.build.configured_libraries.{name}"
        )
        paths = observation["value"]
        if paths is not None:
            if not isinstance(paths, list) or any(
                not isinstance(path, str) or not path for path in paths
            ):
                raise ProvenanceError(
                    f"execution.build.configured_libraries.{name}.value must be a path array"
                )
            if paths != sorted(set(paths)):
                raise ProvenanceError(
                    f"execution.build.configured_libraries.{name}.value must be sorted and unique"
                )
    _observation(build["mumps_version"], "execution.build.mumps_version")
    return build


def _validate_launchers(value: object) -> Mapping[str, object]:
    launchers = _observation(value, "execution.launchers")
    records = launchers["value"]
    if records is None:
        return launchers
    if (
        not isinstance(records, list)
        or not records
        or len(records) > _MAX_LAUNCHERS
    ):
        raise ProvenanceError("execution.launchers.value must be a nonempty array")
    for index, raw_record in enumerate(records):
        field = f"execution.launchers.value[{index}]"
        record = _exact_object(raw_record, _LAUNCHER_KEYS, field)
        _text(record["name"], f"{field}.name")
        _text(record["kind"], f"{field}.kind")
        _observation(record["executable"], f"{field}.executable", _string_value)
        _observation(record["resolved_executable"], f"{field}.resolved_executable", _string_value)
        arguments = record["arguments"]
        if (
            not isinstance(arguments, list)
            or len(arguments) > _MAX_LAUNCHER_ARGUMENTS
            or any(
                not isinstance(argument, str)
                or not argument
                or "\0" in argument
                or len(argument) > 4096
                for argument in arguments
            )
        ):
            raise ProvenanceError(f"{field}.arguments must be an argv array")
        _observation(record["version"], f"{field}.version", _string_value)
        for name in (
            "available_mpi_slots",
            "available_cpu_slots",
            "allocated_threads_per_rank",
        ):
            _observation(record[name], f"{field}.{name}", _positive_integer_value)
    if records != sorted(records, key=_canonical_json):
        raise ProvenanceError("execution.launchers.value must be sorted")
    return launchers


def validate_execution_provenance(
    document: Mapping[str, object],
    *,
    manifest: Mapping[str, object] | None = None,
) -> None:
    """Strictly validate an execution record and all cryptographic bindings."""

    execution = _exact_object(document, _EXECUTION_KEYS, "execution")
    if (
        type(execution["schema_version"]) is not int
        or execution["schema_version"] != EXECUTION_SCHEMA_VERSION
    ):
        raise ProvenanceError(
            f"execution.schema_version must be {EXECUTION_SCHEMA_VERSION}"
        )
    if execution["kind"] != EXECUTION_KIND:
        raise ProvenanceError(f"execution.kind must be {EXECUTION_KIND}")
    _text(execution["run_id"], "execution.run_id")
    for name in (
        "manifest_hash",
        "config_hash",
        "build_identity_hash",
        "execution_compatibility_hash",
        "record_hash",
    ):
        value = execution[name]
        if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
            raise ProvenanceError(f"execution.{name} must be a SHA-256 digest")
    recorded_at = _text(execution["recorded_at"], "execution.recorded_at")
    try:
        parsed_time = datetime.fromisoformat(recorded_at)
    except ValueError as error:
        raise ProvenanceError("execution.recorded_at is not ISO-8601") from error
    if parsed_time.tzinfo is None or parsed_time.utcoffset() is None:
        raise ProvenanceError("execution.recorded_at must include a UTC offset")

    timezone_document = _exact_object(
        execution["timezone"], _TIMEZONE_KEYS, "execution.timezone"
    )
    _observation(timezone_document["name"], "execution.timezone.name", _string_value)
    _observation(
        timezone_document["utc_offset_seconds"],
        "execution.timezone.utc_offset_seconds",
        _integer_value,
    )
    system = _exact_object(execution["system"], _SYSTEM_KEYS, "execution.system")
    for name in (
        "hostname",
        "operating_system",
        "kernel",
        "architecture",
        "cpu_model",
    ):
        _observation(system[name], f"execution.system.{name}", _string_value)
    for name in ("physical_cores", "logical_cores", "memory_bytes"):
        _observation(system[name], f"execution.system.{name}", _positive_integer_value)

    topology = _exact_object(
        execution["topology"], _TOPOLOGY_KEYS, "execution.topology"
    )
    rank_grids = topology["rank_grids"]
    if not isinstance(rank_grids, list) or not rank_grids:
        raise ProvenanceError("execution.topology.rank_grids must be nonempty")
    previous_grid: tuple[int, int, int, int] | None = None
    for index, raw_grid in enumerate(rank_grids):
        field = f"execution.topology.rank_grids[{index}]"
        grid = _exact_object(raw_grid, _RANK_GRID_KEYS, field)
        ranks = grid["ranks"]
        vector = grid["process_grid"]
        if type(ranks) is not int or ranks <= 0:
            raise ProvenanceError(f"{field}.ranks must be positive")
        if (
            not isinstance(vector, list)
            or len(vector) != 3
            or any(type(item) is not int or item <= 0 for item in vector)
        ):
            raise ProvenanceError(f"{field}.process_grid is invalid")
        coordinate = (ranks, vector[0], vector[1], vector[2])
        if previous_grid is not None and coordinate <= previous_grid:
            raise ProvenanceError("execution.topology.rank_grids must be sorted and unique")
        previous_grid = coordinate
    bindings = topology["bindings"]
    if not isinstance(bindings, list) or not bindings:
        raise ProvenanceError("execution.topology.bindings must be nonempty")
    canonical_bindings: list[str] = []
    for index, raw_binding in enumerate(bindings):
        field = f"execution.topology.bindings[{index}]"
        binding = _exact_object(raw_binding, _BINDING_KEYS, field)
        _positive_integer_value(binding["threads"], f"{field}.threads")
        if type(binding["dynamic"]) is not bool:
            raise ProvenanceError(f"{field}.dynamic must be boolean")
        _text(binding["proc_bind"], f"{field}.proc_bind")
        _text(binding["places"], f"{field}.places")
        canonical_bindings.append(_canonical_json(binding))
    if canonical_bindings != sorted(set(canonical_bindings)):
        raise ProvenanceError("execution.topology.bindings must be sorted and unique")

    environment = execution["environment"]
    if not isinstance(environment, dict) or len(environment) > _MAX_ENVIRONMENT_VARIABLES:
        raise ProvenanceError("execution.environment must be a bounded object")
    for name, observation in environment.items():
        if _OMP_NAME.fullmatch(name) is None:
            raise ProvenanceError(f"execution.environment has invalid key {name!r}")
        _observation(observation, f"execution.environment.{name}", _string_value)

    build = _validate_build(execution["build"])
    launchers = _validate_launchers(execution["launchers"])
    expected_build_hash = _digest(
        {"format": "ads-build-identity-v1", "build": build}
    )
    if execution["build_identity_hash"] != expected_build_hash:
        raise ProvenanceError("execution.build_identity_hash does not match build")
    expected_compatibility_hash = _digest(
        {
            "format": "ads-execution-compatibility-v1",
            "build_identity_hash": expected_build_hash,
            "launchers": launchers,
        }
    )
    if execution["execution_compatibility_hash"] != expected_compatibility_hash:
        raise ProvenanceError(
            "execution.execution_compatibility_hash does not match build and launchers"
        )
    without_record_hash = dict(execution)
    without_record_hash.pop("record_hash")
    expected_record_hash = _digest(
        {"format": "ads-execution-record-v1", "record": without_record_hash}
    )
    if execution["record_hash"] != expected_record_hash:
        raise ProvenanceError("execution.record_hash does not match record")

    if manifest is not None:
        run_id, config_hash, _, expected_grids, expected_bindings = _manifest_matrix(
            manifest
        )
        if execution["run_id"] != run_id:
            raise ProvenanceError("execution.run_id does not match manifest")
        if execution["config_hash"] != config_hash:
            raise ProvenanceError("execution.config_hash does not match manifest")
        if execution["manifest_hash"] != manifest_digest(manifest):
            raise ProvenanceError("execution.manifest_hash does not match manifest")
        if topology["rank_grids"] != expected_grids:
            raise ProvenanceError("execution rank grids do not match manifest")
        if topology["bindings"] != expected_bindings:
            raise ProvenanceError("execution bindings do not match manifest")


def _comparison_roots(
    execution: Mapping[str, object],
) -> tuple[tuple[Path, ...], tuple[Path, ...]]:
    """Infer owned repository/build prefixes only from standard artifact paths."""

    build = execution["build"]
    assert isinstance(build, Mapping)
    units = build["units"]
    assert isinstance(units, list)
    build_roots: set[Path] = set()
    repository_roots: set[Path] = set()
    for unit in units:
        assert isinstance(unit, Mapping)
        for field in ("executable", "core_stamp", "adapter_stamp"):
            identity = unit[field]
            assert isinstance(identity, Mapping)
            raw_path = identity["path"]
            assert isinstance(raw_path, str)
            path = Path(raw_path)
            if not path.is_absolute():
                continue
            for parent in path.parents:
                if parent.name == "build" and parent.parent.name == "benchmarking":
                    build_roots.add(parent)
                    repository_roots.add(parent.parent.parent)
                    break
    return (
        tuple(sorted(build_roots, key=lambda item: (-len(str(item)), str(item)))),
        tuple(
            sorted(repository_roots, key=lambda item: (-len(str(item)), str(item)))
        ),
    )


def _comparison_path(
    value: str,
    build_roots: Sequence[Path],
    repository_roots: Sequence[Path],
) -> str:
    path = Path(value)
    if not path.is_absolute():
        return value
    for root in build_roots:
        try:
            relative = path.relative_to(root)
        except ValueError:
            continue
        suffix = relative.as_posix()
        return (
            "<BENCHMARK_BUILD_ROOT>"
            if suffix == "."
            else f"<BENCHMARK_BUILD_ROOT>/{suffix}"
        )
    for root in repository_roots:
        try:
            relative = path.relative_to(root)
        except ValueError:
            continue
        suffix = relative.as_posix()
        return (
            "<REPOSITORY_ROOT>"
            if suffix == "."
            else f"<REPOSITORY_ROOT>/{suffix}"
        )
    return value


def _comparison_command(
    value: str,
    build_roots: Sequence[Path],
    repository_roots: Sequence[Path],
) -> str:
    """Canonicalize argv spelling and owned -I/-J/-L/-module paths."""

    try:
        arguments = shlex.split(value)
    except ValueError:
        # The strict record still preserves the unparsable value.  Keeping it
        # verbatim is safer than guessing how shell quoting was intended.
        return value
    normalized: list[str] = []
    index = 0
    while index < len(arguments):
        argument = arguments[index]
        if argument in {"-I", "-J", "-L", "-module"}:
            normalized.append(argument)
            if index + 1 < len(arguments):
                index += 1
                normalized.append(
                    _comparison_path(
                        arguments[index], build_roots, repository_roots
                    )
                )
        else:
            attached = False
            for prefix in ("-I", "-J", "-L"):
                if argument.startswith(prefix) and len(argument) > len(prefix):
                    normalized.append(
                        prefix
                        + _comparison_path(
                            argument[len(prefix) :],
                            build_roots,
                            repository_roots,
                        )
                    )
                    attached = True
                    break
            if not attached:
                normalized.append(
                    _comparison_path(argument, build_roots, repository_roots)
                )
        index += 1
    return shlex.join(normalized)


def _comparison_observation(
    value: Mapping[str, object],
    transform: Callable[[object], object],
) -> dict[str, object]:
    observed = value["value"]
    return {
        "value": None if observed is None else transform(observed),
        "reason": value["reason"],
    }


def _comparison_owned_texts(
    value: object,
    build_roots: Sequence[Path],
    repository_roots: Sequence[Path],
) -> object:
    """Remove owned absolute prefixes even when embedded in diagnostics."""

    if isinstance(value, str):
        result = value
        for root, token in (
            *((root, "<BENCHMARK_BUILD_ROOT>") for root in build_roots),
            *((root, "<REPOSITORY_ROOT>") for root in repository_roots),
        ):
            prefix = str(root)
            if result == prefix:
                result = token
            result = result.replace(prefix + os.sep, token + "/")
        return result
    if isinstance(value, list):
        return [
            _comparison_owned_texts(item, build_roots, repository_roots)
            for item in value
        ]
    if isinstance(value, dict):
        return {
            key: _comparison_owned_texts(item, build_roots, repository_roots)
            for key, item in value.items()
        }
    return value


def comparison_environment_identity(
    record: Mapping[str, object],
) -> dict[str, object]:
    """Return the path-neutral machine/toolchain identity used for A/B.

    This intentionally excludes run/time/manifest identity and the content of
    adapter executables and internal build stamps.  Two isolated worktrees can
    therefore be compared when their compiler, flags, external libraries,
    launcher, topology, OpenMP policy, and machine are equivalent.  External
    library paths and content hashes remain part of the identity.
    """

    validate_execution_provenance(record)
    execution = record
    build_roots, repository_roots = _comparison_roots(execution)
    raw_build = execution["build"]
    assert isinstance(raw_build, Mapping)

    units: list[dict[str, object]] = []
    raw_units = raw_build["units"]
    assert isinstance(raw_units, list)
    for raw_unit in raw_units:
        assert isinstance(raw_unit, Mapping)

        def normalized_command(value: object) -> object:
            assert isinstance(value, str)
            return _comparison_command(value, build_roots, repository_roots)

        units.append(
            {
                "profile": raw_unit["profile"],
                "problem": raw_unit["problem"],
                "compiler_command": _comparison_observation(
                    raw_unit["compiler_command"], normalized_command  # type: ignore[arg-type]
                ),
                "compile_flags": _comparison_observation(
                    raw_unit["compile_flags"], normalized_command  # type: ignore[arg-type]
                ),
                "link_flags": _comparison_observation(
                    raw_unit["link_flags"], normalized_command  # type: ignore[arg-type]
                ),
            }
        )

    def normalized_path(value: object) -> object:
        assert isinstance(value, str)
        return _comparison_path(value, build_roots, repository_roots)

    compiler_versions_observation = raw_build["compiler_versions"]
    assert isinstance(compiler_versions_observation, Mapping)
    raw_versions = compiler_versions_observation["value"]
    normalized_versions: list[dict[str, object]] | None = None
    if raw_versions is not None:
        assert isinstance(raw_versions, list)
        normalized_versions = []
        for version in raw_versions:
            assert isinstance(version, Mapping)
            command = version["command"]
            assert isinstance(command, str)
            resolved = version["resolved_executable"]
            assert isinstance(resolved, Mapping)
            normalized_versions.append(
                {
                    "command": _comparison_command(
                        command, build_roots, repository_roots
                    ),
                    "resolved_executable": _comparison_observation(
                        resolved, normalized_path
                    ),
                    "version": deepcopy(version["version"]),
                }
            )

    raw_libraries = raw_build["libraries"]
    assert isinstance(raw_libraries, list)
    libraries: list[dict[str, object]] = []
    for raw_library in raw_libraries:
        assert isinstance(raw_library, Mapping)
        path = raw_library["path"]
        assert isinstance(path, str)
        libraries.append(
            {
                "path": _comparison_path(
                    path, build_roots, repository_roots
                ),
                "size_bytes": raw_library["size_bytes"],
                "sha256": raw_library["sha256"],
                "reason": raw_library["reason"],
            }
        )
    libraries.sort(key=lambda item: str(item["path"]))

    raw_configured = raw_build["configured_libraries"]
    assert isinstance(raw_configured, Mapping)
    configured: dict[str, object] = {}
    for name in _IMPORTANT_LIBRARY_KEYS:
        observation = raw_configured[name]
        assert isinstance(observation, Mapping)

        def normalized_paths(value: object) -> object:
            assert isinstance(value, list)
            return sorted(
                _comparison_path(path, build_roots, repository_roots)
                for path in value
            )

        configured[name] = _comparison_observation(
            observation, normalized_paths
        )

    mumps_observation = raw_build["mumps_version"]
    assert isinstance(mumps_observation, Mapping)

    def normalized_mumps(value: object) -> object:
        if not isinstance(value, dict):
            return deepcopy(value)
        normalized = deepcopy(value)
        header = normalized.get("header")
        if isinstance(header, str):
            normalized["header"] = _comparison_path(
                header, build_roots, repository_roots
            )
        return normalized

    raw_launchers = execution["launchers"]
    assert isinstance(raw_launchers, Mapping)

    def normalized_launchers(value: object) -> object:
        assert isinstance(value, list)
        launchers: list[dict[str, object]] = []
        for raw_launcher in value:
            assert isinstance(raw_launcher, Mapping)
            launcher = deepcopy(dict(raw_launcher))
            for name in ("executable", "resolved_executable"):
                observation = launcher[name]
                assert isinstance(observation, Mapping)
                launcher[name] = _comparison_observation(
                    observation, normalized_path
                )
            arguments = launcher["arguments"]
            assert isinstance(arguments, list)
            launcher["arguments"] = [
                _comparison_path(argument, build_roots, repository_roots)
                for argument in arguments
            ]
            launchers.append(launcher)
        launchers.sort(key=_canonical_json)
        return launchers

    identity = {
        "schema_version": 1,
        "kind": "ads-benchmark-comparison-environment",
        "system": deepcopy(execution["system"]),
        "topology": deepcopy(execution["topology"]),
        "environment": deepcopy(execution["environment"]),
        "build": {
            "profiles": deepcopy(raw_build["profiles"]),
            "units": units,
            "compiler_versions": {
                "value": normalized_versions,
                "reason": compiler_versions_observation["reason"],
            },
            "libraries": libraries,
            "configured_libraries": configured,
            "mumps_version": _comparison_observation(
                mumps_observation, normalized_mumps
            ),
        },
        "launchers": _comparison_observation(
            raw_launchers, normalized_launchers
        ),
    }
    normalized_identity = _comparison_owned_texts(
        identity, build_roots, repository_roots
    )
    assert isinstance(normalized_identity, dict)
    return normalized_identity
