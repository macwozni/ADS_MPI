"""Contained, atomic storage for plans and per-case execution records."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from contextlib import contextmanager
import fcntl
import json
import os
from pathlib import Path
import re
import stat
import time

from .errors import StorageError


_SAFE_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$", re.ASCII)
_PLAN_KIND = "ads-benchmark-plan"
_PLAN_SCHEMA_VERSION = 1
_DIRECTORY_FLAGS = (
    os.O_RDONLY
    | getattr(os, "O_DIRECTORY", 0)
    | getattr(os, "O_NOFOLLOW", 0)
    | getattr(os, "O_CLOEXEC", 0)
)
_FILE_NOFOLLOW = getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_CLOEXEC", 0)
_MAX_JSON_BYTES = 64 * 1024 * 1024
_MAX_LOG_BYTES = 64 * 1024 * 1024
_MAX_GENERATED_ARTIFACT_BYTES = 64 * 1024 * 1024
_LOG_NAMES = frozenset({"stdout.log", "stderr.log"})
_GENERATED_ARTIFACTS = frozenset({"field_samples.csv"})
_ANALYSIS_ARTIFACTS = frozenset(
    {"analysis.json", "analysis.csv", "convergence.png", "strong-scaling.png"}
)


def validate_identifier(value: str, field: str) -> str:
    if (
        not isinstance(value, str)
        or value in {".", ".."}
        or not _SAFE_IDENTIFIER.fullmatch(value)
    ):
        raise StorageError(
            f"unsafe {field}: {value!r}; use 1-128 ASCII letters, digits, "
            "'.', '_', '-'"
        )
    return value


class ResultStore:
    """Own only verified run directories below ``repo/benchmarks``.

    Directory file descriptors plus ``O_NOFOLLOW`` keep every create/write
    operation anchored below the results root even if another process swaps a
    path component concurrently.  The shared results root itself is never
    marked or cleaned; ownership begins at an exclusively created run whose
    manifest has the expected kind, schema, and run ID.
    """

    def __init__(
        self, repository_root: Path, results_root: Path | None = None
    ) -> None:
        if not hasattr(os, "O_NOFOLLOW") or not hasattr(os, "O_DIRECTORY"):
            raise StorageError("safe result storage requires POSIX O_NOFOLLOW support")
        self.repository_root = repository_root.resolve()
        self.expected_root = self.repository_root / "benchmarks"
        candidate = results_root or self.expected_root
        if not candidate.is_absolute():
            candidate = self.repository_root / candidate
        if candidate.is_symlink() or candidate.resolve() != self.expected_root:
            raise StorageError(
                f"results root must be exactly {self.expected_root}, not {candidate}"
            )
        self.results_root = candidate
        if self.results_root.exists() and not self.results_root.is_dir():
            raise StorageError(f"results root is not a directory: {self.results_root}")

    def _run_path(self, run_id: str) -> Path:
        validate_identifier(run_id, "run_id")
        return self.results_root / run_id

    @staticmethod
    def _validate_manifest(
        manifest: Mapping[str, object], run_id: str
    ) -> None:
        if (
            type(manifest.get("schema_version")) is not int
            or manifest.get("schema_version") != _PLAN_SCHEMA_VERSION
            or manifest.get("kind") != _PLAN_KIND
            or manifest.get("run_id") != run_id
        ):
            raise StorageError(
                f"run {run_id} lacks a matching {_PLAN_KIND} ownership manifest"
            )

    def _results_fd(self, *, create: bool) -> int:
        try:
            repository_fd = os.open(self.repository_root, _DIRECTORY_FLAGS)
        except (OSError, UnicodeError) as error:
            raise StorageError(
                f"cannot open repository root {self.repository_root}: {error}"
            ) from error
        try:
            if create:
                try:
                    os.mkdir("benchmarks", mode=0o755, dir_fd=repository_fd)
                except FileExistsError:
                    pass
            try:
                return os.open("benchmarks", _DIRECTORY_FLAGS, dir_fd=repository_fd)
            except OSError as error:
                raise StorageError(
                    f"results root is missing or unsafe: {self.results_root}: {error}"
                ) from error
        finally:
            os.close(repository_fd)

    @staticmethod
    def _atomic_text_at(directory_fd: int, name: str, content: str) -> None:
        temporary = f".{name}.{os.getpid()}.{time.time_ns()}.tmp"
        descriptor: int | None = None
        try:
            descriptor = os.open(
                temporary,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | _FILE_NOFOLLOW,
                0o644,
                dir_fd=directory_fd,
            )
            with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
                descriptor = None
                stream.write(content)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(
                temporary,
                name,
                src_dir_fd=directory_fd,
                dst_dir_fd=directory_fd,
            )
            os.fsync(directory_fd)
        except (OSError, UnicodeError) as error:
            if descriptor is not None:
                os.close(descriptor)
            try:
                os.unlink(temporary, dir_fd=directory_fd)
            except OSError:
                pass
            raise StorageError(f"cannot write {name}: {error}") from error

    @staticmethod
    def _atomic_bytes_at(directory_fd: int, name: str, content: bytes) -> None:
        temporary = f".{name}.{os.getpid()}.{time.time_ns()}.tmp"
        descriptor: int | None = None
        try:
            descriptor = os.open(
                temporary,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | _FILE_NOFOLLOW,
                0o644,
                dir_fd=directory_fd,
            )
            with os.fdopen(descriptor, "wb") as stream:
                descriptor = None
                stream.write(content)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(
                temporary,
                name,
                src_dir_fd=directory_fd,
                dst_dir_fd=directory_fd,
            )
            os.fsync(directory_fd)
        except (OSError, UnicodeError) as error:
            if descriptor is not None:
                os.close(descriptor)
            try:
                os.unlink(temporary, dir_fd=directory_fd)
            except OSError:
                pass
            raise StorageError(f"cannot write {name}: {error}") from error

    @classmethod
    def _atomic_json_at(
        cls, directory_fd: int, name: str, document: Mapping[str, object]
    ) -> None:
        try:
            content = json.dumps(
                document,
                indent=2,
                sort_keys=True,
                ensure_ascii=True,
                allow_nan=False,
            ) + "\n"
        except (TypeError, ValueError) as error:
            raise StorageError(f"record for {name} is not valid JSON: {error}") from error
        cls._atomic_text_at(directory_fd, name, content)

    @staticmethod
    def _read_json_at(
        directory_fd: int,
        name: str,
        description: str,
        *,
        missing_ok: bool = False,
    ) -> dict[str, object] | None:
        descriptor: int | None = None
        try:
            descriptor = os.open(
                name,
                os.O_RDONLY | os.O_NONBLOCK | _FILE_NOFOLLOW,
                dir_fd=directory_fd,
            )
            metadata = os.fstat(descriptor)
            if (
                not stat.S_ISREG(metadata.st_mode)
                or metadata.st_size > _MAX_JSON_BYTES
            ):
                raise StorageError(
                    f"{description} must be a regular file below 64 MiB"
                )

            def unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
                document: dict[str, object] = {}
                for key, value in pairs:
                    if key in document:
                        raise StorageError(
                            f"{description} contains duplicate JSON key {key!r}"
                        )
                    document[key] = value
                return document

            def reject_constant(value: str) -> None:
                raise StorageError(
                    f"{description} contains non-finite JSON number {value}"
                )

            with os.fdopen(descriptor, "r", encoding="utf-8") as stream:
                descriptor = None
                document = json.load(
                    stream,
                    object_pairs_hook=unique_object,
                    parse_constant=reject_constant,
                )
        except StorageError:
            if descriptor is not None:
                os.close(descriptor)
            raise
        except FileNotFoundError:
            if descriptor is not None:
                os.close(descriptor)
            if missing_ok:
                return None
            raise StorageError(f"{description} is missing") from None
        except (
            OSError,
            UnicodeError,
            json.JSONDecodeError,
            ValueError,
            RecursionError,
        ) as error:
            if descriptor is not None:
                os.close(descriptor)
            raise StorageError(f"cannot read {description}: {error}") from error
        if not isinstance(document, dict):
            raise StorageError(f"{description} is not an object")
        return document

    @classmethod
    def _read_manifest_at(cls, run_fd: int, run_id: str) -> dict[str, object]:
        document = cls._read_json_at(
            run_fd,
            "manifest.json",
            f"run {run_id} ownership manifest",
        )
        assert document is not None
        return document

    @contextmanager
    def _owned_run_fds(self, run_id: str) -> Iterator[tuple[int, int]]:
        run_path = self._run_path(run_id)
        results_fd = self._results_fd(create=False)
        run_fd: int | None = None
        try:
            try:
                run_fd = os.open(run_id, _DIRECTORY_FLAGS, dir_fd=results_fd)
            except OSError as error:
                raise StorageError(f"run does not exist or is unsafe: {run_path}") from error
            manifest = self._read_manifest_at(run_fd, run_id)
            self._validate_manifest(manifest, run_id)
            yield results_fd, run_fd
        finally:
            if run_fd is not None:
                os.close(run_fd)
            os.close(results_fd)

    def create_run(
        self, run_id: str, manifest: Mapping[str, object]
    ) -> Path:
        run_path = self._run_path(run_id)
        self._validate_manifest(manifest, run_id)
        results_fd = self._results_fd(create=True)
        run_fd: int | None = None
        created = False
        try:
            try:
                os.mkdir(run_id, mode=0o755, dir_fd=results_fd)
                created = True
            except FileExistsError as error:
                raise StorageError(
                    f"run already exists; refusing overwrite: {run_path}"
                ) from error
            run_fd = os.open(run_id, _DIRECTORY_FLAGS, dir_fd=results_fd)
            self._atomic_json_at(run_fd, "manifest.json", manifest)
            return run_path
        except Exception:
            if created:
                if run_fd is not None:
                    try:
                        os.unlink("manifest.json", dir_fd=run_fd)
                    except OSError:
                        pass
                try:
                    os.rmdir(run_id, dir_fd=results_fd)
                except OSError:
                    pass
            raise
        finally:
            if run_fd is not None:
                os.close(run_fd)
            os.close(results_fd)

    def read_manifest(self, run_id: str) -> dict[str, object]:
        """Read an owned run's immutable execution manifest safely."""

        with self._owned_run_fds(run_id) as (_, run_fd):
            return self._read_manifest_at(run_fd, run_id)

    def write_analysis_artifact(
        self, run_id: str, name: str, content: str | bytes
    ) -> Path:
        """Atomically write one allowlisted artifact below an owned run."""

        if name not in _ANALYSIS_ARTIFACTS:
            raise StorageError(f"unsupported analysis artifact: {name}")
        if name.endswith(".png") and not isinstance(content, bytes):
            raise StorageError(f"{name} content must be bytes")
        if not name.endswith(".png") and not isinstance(content, str):
            raise StorageError(f"{name} content must be text")
        with self._owned_run_fds(run_id) as (_, run_fd):
            try:
                os.mkdir("analysis", mode=0o755, dir_fd=run_fd)
            except FileExistsError:
                pass
            try:
                analysis_fd = os.open("analysis", _DIRECTORY_FLAGS, dir_fd=run_fd)
            except OSError as error:
                raise StorageError(
                    f"unsafe analysis directory in run {run_id}"
                ) from error
            try:
                if isinstance(content, bytes):
                    self._atomic_bytes_at(analysis_fd, name, content)
                else:
                    self._atomic_text_at(analysis_fd, name, content)
            finally:
                os.close(analysis_fd)
        return self._run_path(run_id) / "analysis" / name

    def remove_analysis_artifact(self, run_id: str, name: str) -> None:
        """Safely remove one stale allowlisted analysis artifact."""

        if name not in _ANALYSIS_ARTIFACTS:
            raise StorageError(f"unsupported analysis artifact: {name}")
        with self._owned_run_fds(run_id) as (_, run_fd):
            try:
                analysis_fd = os.open("analysis", _DIRECTORY_FLAGS, dir_fd=run_fd)
            except FileNotFoundError:
                return
            except OSError as error:
                raise StorageError(
                    f"unsafe analysis directory in run {run_id}"
                ) from error
            try:
                try:
                    os.unlink(name, dir_fd=analysis_fd)
                    os.fsync(analysis_fd)
                except FileNotFoundError:
                    pass
                except OSError as error:
                    raise StorageError(
                        f"cannot remove analysis artifact {name}: {error}"
                    ) from error
            finally:
                os.close(analysis_fd)

    @contextmanager
    def execution_lock(self, run_id: str) -> Iterator[None]:
        """Prevent two runner processes from appending to one run concurrently."""

        descriptor: int | None = None
        with self._owned_run_fds(run_id) as (_, run_fd):
            try:
                descriptor = os.open(
                    ".execution.lock",
                    os.O_RDWR | os.O_CREAT | _FILE_NOFOLLOW,
                    0o644,
                    dir_fd=run_fd,
                )
                metadata = os.fstat(descriptor)
                if not stat.S_ISREG(metadata.st_mode):
                    raise StorageError(
                        f"run {run_id} execution lock is not a regular file"
                    )
                try:
                    fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError as error:
                    raise StorageError(
                        f"run {run_id} is already being executed or resumed"
                    ) from error
                yield
            except StorageError:
                raise
            except OSError as error:
                raise StorageError(
                    f"cannot lock run {run_id} for execution: {error}"
                ) from error
            finally:
                if descriptor is not None:
                    try:
                        fcntl.flock(descriptor, fcntl.LOCK_UN)
                    finally:
                        os.close(descriptor)

    def case_directory_path(self, run_id: str, case_id: str) -> Path:
        """Return the validated lexical path for a case without creating it."""

        run_path = self._run_path(run_id)
        validate_identifier(case_id, "case_id")
        return run_path / "cases" / case_id

    def create_case_directory(self, run_id: str, case_id: str) -> Path:
        run_path = self._run_path(run_id)
        validate_identifier(case_id, "case_id")
        with self._owned_run_fds(run_id) as (_, run_fd):
            try:
                os.mkdir("cases", mode=0o755, dir_fd=run_fd)
            except FileExistsError:
                pass
            try:
                cases_fd = os.open("cases", _DIRECTORY_FLAGS, dir_fd=run_fd)
            except OSError as error:
                raise StorageError(f"unsafe cases directory in run {run_id}") from error
            try:
                try:
                    os.mkdir(case_id, mode=0o755, dir_fd=cases_fd)
                except FileExistsError as error:
                    raise StorageError(
                        f"case already exists; refusing overwrite: {case_id}"
                    ) from error
            finally:
                os.close(cases_fd)
        return run_path / "cases" / case_id

    def ensure_case_directory(self, run_id: str, case_id: str) -> tuple[Path, bool]:
        """Return a safely anchored case directory, creating it when absent."""

        run_path = self._run_path(run_id)
        validate_identifier(case_id, "case_id")
        created = False
        with self._owned_run_fds(run_id) as (_, run_fd):
            try:
                os.mkdir("cases", mode=0o755, dir_fd=run_fd)
            except FileExistsError:
                pass
            try:
                cases_fd = os.open("cases", _DIRECTORY_FLAGS, dir_fd=run_fd)
            except OSError as error:
                raise StorageError(f"unsafe cases directory in run {run_id}") from error
            case_fd: int | None = None
            try:
                try:
                    os.mkdir(case_id, mode=0o755, dir_fd=cases_fd)
                    created = True
                except FileExistsError:
                    pass
                try:
                    case_fd = os.open(case_id, _DIRECTORY_FLAGS, dir_fd=cases_fd)
                except OSError as error:
                    raise StorageError(
                        f"case path is missing or unsafe: {case_id}"
                    ) from error
            finally:
                if case_fd is not None:
                    os.close(case_fd)
                os.close(cases_fd)
        return run_path / "cases" / case_id, created

    def _case_coordinates(self, case_directory: Path) -> tuple[str, str]:
        if not case_directory.is_absolute():
            raise StorageError(f"case path must be absolute: {case_directory}")
        cases_path = case_directory.parent
        run_path = cases_path.parent
        run_id = run_path.name
        case_id = case_directory.name
        validate_identifier(run_id, "run_id")
        validate_identifier(case_id, "case_id")
        expected = self.results_root / run_id / "cases" / case_id
        if cases_path.name != "cases" or run_path.parent != self.results_root or case_directory != expected:
            raise StorageError(f"case path has invalid result layout: {case_directory}")
        return run_id, case_id

    @contextmanager
    def _case_fd(self, case_directory: Path) -> Iterator[int]:
        run_id, case_id = self._case_coordinates(case_directory)
        with self._owned_run_fds(run_id) as (_, run_fd):
            try:
                cases_fd = os.open("cases", _DIRECTORY_FLAGS, dir_fd=run_fd)
            except OSError as error:
                raise StorageError(f"unsafe cases directory in run {run_id}") from error
            case_fd: int | None = None
            try:
                try:
                    case_fd = os.open(case_id, _DIRECTORY_FLAGS, dir_fd=cases_fd)
                except OSError as error:
                    raise StorageError(
                        f"case path is missing or unsafe: {case_directory}"
                    ) from error
                yield case_fd
            finally:
                if case_fd is not None:
                    os.close(case_fd)
                os.close(cases_fd)

    @contextmanager
    def open_case_directory(self, case_directory: Path) -> Iterator[int]:
        """Keep a verified case-directory inode open across child launch."""

        with self._case_fd(case_directory) as case_fd:
            yield case_fd

    def write_status(self, case_directory: Path, document: Mapping[str, object]) -> None:
        with self._case_fd(case_directory) as case_fd:
            self._atomic_json_at(case_fd, "status.json", document)

    def write_result(self, case_directory: Path, document: Mapping[str, object]) -> None:
        with self._case_fd(case_directory) as case_fd:
            self._atomic_json_at(case_fd, "result.json", document)

    def write_log(self, case_directory: Path, name: str, content: str) -> None:
        if name not in _LOG_NAMES:
            raise StorageError(f"unsupported log name: {name}")
        with self._case_fd(case_directory) as case_fd:
            self._atomic_text_at(case_fd, name, content)

    def read_case_json(
        self,
        case_directory: Path,
        name: str,
        *,
        missing_ok: bool = False,
    ) -> dict[str, object] | None:
        if name not in {"status.json", "result.json"}:
            raise StorageError(f"unsupported case JSON name: {name}")
        _, case_id = self._case_coordinates(case_directory)
        with self._case_fd(case_directory) as case_fd:
            return self._read_json_at(
                case_fd,
                name,
                f"case {case_id} {name}",
                missing_ok=missing_ok,
            )

    def case_log_is_complete(self, case_directory: Path, name: str) -> bool:
        if name not in _LOG_NAMES:
            raise StorageError(f"unsupported log name: {name}")
        return self.case_file_is_regular(case_directory, name)

    def case_file_is_regular(self, case_directory: Path, name: str) -> bool:
        if name not in _LOG_NAMES | _GENERATED_ARTIFACTS:
            raise StorageError(f"unsupported case file name: {name}")
        _, case_id = self._case_coordinates(case_directory)
        descriptor: int | None = None
        with self._case_fd(case_directory) as case_fd:
            try:
                descriptor = os.open(
                    name,
                    os.O_RDONLY | os.O_NONBLOCK | _FILE_NOFOLLOW,
                    dir_fd=case_fd,
                )
                return stat.S_ISREG(os.fstat(descriptor).st_mode)
            except FileNotFoundError:
                return False
            except OSError as error:
                raise StorageError(
                    f"cannot inspect case {case_id} {name}: {error}"
                ) from error
            finally:
                if descriptor is not None:
                    os.close(descriptor)

    def read_case_log(self, case_directory: Path, name: str) -> str:
        if name not in _LOG_NAMES:
            raise StorageError(f"unsupported log name: {name}")
        _, case_id = self._case_coordinates(case_directory)
        descriptor: int | None = None
        try:
            with self._case_fd(case_directory) as case_fd:
                descriptor = os.open(
                    name,
                    os.O_RDONLY | os.O_NONBLOCK | _FILE_NOFOLLOW,
                    dir_fd=case_fd,
                )
                metadata = os.fstat(descriptor)
                if (
                    not stat.S_ISREG(metadata.st_mode)
                    or metadata.st_size > _MAX_LOG_BYTES
                ):
                    raise StorageError(
                        f"case {case_id} {name} must be a regular file below 64 MiB"
                    )
                with os.fdopen(descriptor, "r", encoding="utf-8") as stream:
                    descriptor = None
                    return stream.read()
        except StorageError:
            raise
        except (OSError, UnicodeError) as error:
            raise StorageError(
                f"cannot read case {case_id} {name}: {error}"
            ) from error
        finally:
            if descriptor is not None:
                os.close(descriptor)

    def read_generated_artifact(self, case_directory: Path, name: str) -> str:
        """Read one allowlisted, bounded text artifact from an owned case."""

        if name not in _GENERATED_ARTIFACTS:
            raise StorageError(f"unsupported generated artifact: {name}")
        _, case_id = self._case_coordinates(case_directory)
        descriptor: int | None = None
        try:
            with self._case_fd(case_directory) as case_fd:
                descriptor = os.open(
                    name,
                    os.O_RDONLY | os.O_NONBLOCK | _FILE_NOFOLLOW,
                    dir_fd=case_fd,
                )
                metadata = os.fstat(descriptor)
                if (
                    not stat.S_ISREG(metadata.st_mode)
                    or metadata.st_size > _MAX_GENERATED_ARTIFACT_BYTES
                ):
                    raise StorageError(
                        f"case {case_id} {name} must be a regular file below 64 MiB"
                    )
                with os.fdopen(descriptor, "rb") as stream:
                    descriptor = None
                    content = stream.read(_MAX_GENERATED_ARTIFACT_BYTES + 1)
                if len(content) > _MAX_GENERATED_ARTIFACT_BYTES:
                    raise StorageError(
                        f"case {case_id} {name} must be a regular file below 64 MiB"
                    )
                return content.decode("utf-8")
        except StorageError:
            raise
        except (OSError, UnicodeError) as error:
            raise StorageError(
                f"cannot read case {case_id} generated artifact {name}: {error}"
            ) from error
        finally:
            if descriptor is not None:
                os.close(descriptor)

    def remove_result(self, case_directory: Path) -> None:
        """Remove a stale result before a retry; missing results are harmless."""

        with self._case_fd(case_directory) as case_fd:
            try:
                os.unlink("result.json", dir_fd=case_fd)
                os.fsync(case_fd)
            except FileNotFoundError:
                pass
            except OSError as error:
                raise StorageError(
                    f"cannot clear stale result in {case_directory}: {error}"
                ) from error

    def remove_generated_artifact(self, case_directory: Path, name: str) -> None:
        if name not in _GENERATED_ARTIFACTS:
            raise StorageError(f"unsupported generated artifact: {name}")
        with self._case_fd(case_directory) as case_fd:
            try:
                os.unlink(name, dir_fd=case_fd)
                os.fsync(case_fd)
            except FileNotFoundError:
                pass
            except OSError as error:
                raise StorageError(
                    f"cannot clear stale {name} in {case_directory}: {error}"
                ) from error
