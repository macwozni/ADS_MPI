"""Contained, atomic storage for plans and per-case execution records."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from contextlib import contextmanager
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
    def _read_manifest_at(run_fd: int, run_id: str) -> dict[str, object]:
        descriptor: int | None = None
        try:
            descriptor = os.open(
                "manifest.json",
                os.O_RDONLY | os.O_NONBLOCK | _FILE_NOFOLLOW,
                dir_fd=run_fd,
            )
            metadata = os.fstat(descriptor)
            if (
                not stat.S_ISREG(metadata.st_mode)
                or metadata.st_size > 64 * 1024 * 1024
            ):
                raise StorageError(
                    f"run {run_id} ownership manifest must be a regular file below 64 MiB"
                )
            with os.fdopen(descriptor, "r", encoding="utf-8") as stream:
                descriptor = None
                document = json.load(stream)
        except StorageError:
            if descriptor is not None:
                os.close(descriptor)
            raise
        except (OSError, UnicodeError, json.JSONDecodeError) as error:
            if descriptor is not None:
                os.close(descriptor)
            raise StorageError(
                f"run {run_id} has no readable ownership manifest: {error}"
            ) from error
        if not isinstance(document, dict):
            raise StorageError(f"run {run_id} ownership manifest is not an object")
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

    def write_status(self, case_directory: Path, document: Mapping[str, object]) -> None:
        with self._case_fd(case_directory) as case_fd:
            self._atomic_json_at(case_fd, "status.json", document)

    def write_result(self, case_directory: Path, document: Mapping[str, object]) -> None:
        with self._case_fd(case_directory) as case_fd:
            self._atomic_json_at(case_fd, "result.json", document)

    def write_log(self, case_directory: Path, name: str, content: str) -> None:
        if name not in {"stdout.log", "stderr.log"}:
            raise StorageError(f"unsupported log name: {name}")
        with self._case_fd(case_directory) as case_fd:
            self._atomic_text_at(case_fd, name, content)
