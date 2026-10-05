"""Contained, atomic storage for plans and per-case execution records."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping
from contextlib import contextmanager
import fcntl
import json
import os
from pathlib import Path
import re
import stat
import time

from .errors import ProvenanceError, StorageError
from .provenance import validate_execution_provenance


_SAFE_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$", re.ASCII)
_SHARD_FILE = re.compile(r"^shard-([0-9]{6})\.json$", re.ASCII)
_SHARD_LOCK = ".shards.lock"
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
    {
        "analysis.json",
        "analysis.csv",
        "comparison.json",
        "comparison.csv",
        "ab-schedule.json",
        "convergence.png",
        "strong-scaling.png",
        "weak-scaling.png",
    }
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
    def _json_text(cls, name: str, document: Mapping[str, object]) -> str:
        try:
            return json.dumps(
                document,
                indent=2,
                sort_keys=True,
                ensure_ascii=True,
                allow_nan=False,
            ) + "\n"
        except (TypeError, ValueError) as error:
            raise StorageError(f"record for {name} is not valid JSON: {error}") from error

    @classmethod
    def _atomic_json_at(
        cls, directory_fd: int, name: str, document: Mapping[str, object]
    ) -> None:
        content = cls._json_text(name, document)
        cls._atomic_text_at(directory_fd, name, content)

    @classmethod
    def _exclusive_json_at(
        cls, directory_fd: int, name: str, document: Mapping[str, object]
    ) -> None:
        """Durably publish JSON without ever replacing an existing name."""

        content = cls._json_text(name, document)
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
            try:
                # A hard-link publish is the portable POSIX no-replace
                # primitive.  Unlike a preceding stat followed by replace,
                # two writers cannot both claim the same immutable name.
                os.link(
                    temporary,
                    name,
                    src_dir_fd=directory_fd,
                    dst_dir_fd=directory_fd,
                    follow_symlinks=False,
                )
            except FileExistsError as error:
                raise StorageError(
                    f"file already exists; refusing overwrite: {name}"
                ) from error
            os.unlink(temporary, dir_fd=directory_fd)
            os.fsync(directory_fd)
        except StorageError:
            if descriptor is not None:
                os.close(descriptor)
            try:
                os.unlink(temporary, dir_fd=directory_fd)
            except OSError:
                pass
            raise
        except (OSError, UnicodeError) as error:
            if descriptor is not None:
                os.close(descriptor)
            try:
                os.unlink(temporary, dir_fd=directory_fd)
            except OSError:
                pass
            raise StorageError(f"cannot write {name}: {error}") from error

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
        self,
        run_id: str,
        manifest: Mapping[str, object],
        *,
        execution_record: Mapping[str, object] | None = None,
    ) -> Path:
        run_path = self._run_path(run_id)
        self._validate_manifest(manifest, run_id)
        if execution_record is not None:
            try:
                validate_execution_provenance(
                    execution_record, manifest=manifest
                )
            except ProvenanceError as error:
                raise StorageError(f"invalid execution record: {error}") from error
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
            if execution_record is not None:
                self._exclusive_json_at(
                    run_fd, "execution.json", execution_record
                )
            return run_path
        except Exception:
            if created:
                if run_fd is not None:
                    try:
                        os.unlink("manifest.json", dir_fd=run_fd)
                    except OSError:
                        pass
                    try:
                        os.unlink("execution.json", dir_fd=run_fd)
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

    def write_execution_record(
        self, run_id: str, document: Mapping[str, object]
    ) -> Path:
        """Exclusively publish one manifest-bound ``execution.json`` record."""

        with self._owned_run_fds(run_id) as (_, run_fd):
            manifest = self._read_manifest_at(run_fd, run_id)
            try:
                validate_execution_provenance(document, manifest=manifest)
            except ProvenanceError as error:
                raise StorageError(f"invalid execution record: {error}") from error
            self._exclusive_json_at(run_fd, "execution.json", document)
        return self._run_path(run_id) / "execution.json"

    def read_execution_record(
        self, run_id: str, *, missing_ok: bool = False
    ) -> dict[str, object] | None:
        """Read and strictly verify an owned run's ``execution.json`` record."""

        with self._owned_run_fds(run_id) as (_, run_fd):
            manifest = self._read_manifest_at(run_fd, run_id)
            document = self._read_json_at(
                run_fd,
                "execution.json",
                f"run {run_id} execution record",
                missing_ok=missing_ok,
            )
            if document is None:
                return None
            try:
                validate_execution_provenance(document, manifest=manifest)
            except ProvenanceError as error:
                raise StorageError(f"invalid execution record: {error}") from error
            return document

    @staticmethod
    def _shard_name(index: int) -> str:
        if type(index) is not int or index < 0 or index > 999_999:
            raise StorageError("shard index must be between 0 and 999999")
        return f"shard-{index:06d}.json"

    @contextmanager
    def _shard_manifest_lock(
        self, run_id: str, *, exclusive: bool
    ) -> Iterator[int]:
        """Anchor and lock one parent's shard namespace for a whole operation."""

        descriptor: int | None = None
        locked = False
        with self._owned_run_fds(run_id) as (_, run_fd):
            try:
                descriptor = os.open(
                    _SHARD_LOCK,
                    os.O_RDWR | os.O_CREAT | _FILE_NOFOLLOW,
                    0o644,
                    dir_fd=run_fd,
                )
                if not stat.S_ISREG(os.fstat(descriptor).st_mode):
                    raise StorageError(
                        f"run {run_id} shard lock is not a regular file"
                    )
                operation = fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH
                try:
                    fcntl.flock(descriptor, operation | fcntl.LOCK_NB)
                except BlockingIOError as error:
                    raise StorageError(
                        f"run {run_id} shard manifests are currently in use"
                    ) from error
                locked = True
                yield run_fd
            except StorageError:
                raise
            except OSError as error:
                raise StorageError(
                    f"cannot lock shard manifests for run {run_id}: {error}"
                ) from error
            finally:
                if descriptor is not None:
                    try:
                        if locked:
                            fcntl.flock(descriptor, fcntl.LOCK_UN)
                    finally:
                        os.close(descriptor)

    @staticmethod
    def _open_shards_at(run_fd: int, run_id: str, *, create: bool) -> int:
        if create:
            try:
                os.mkdir("shards", mode=0o755, dir_fd=run_fd)
            except FileExistsError:
                pass
        try:
            return os.open("shards", _DIRECTORY_FLAGS, dir_fd=run_fd)
        except OSError as error:
            qualifier = "unsafe" if create else "no safe"
            raise StorageError(
                f"run {run_id} has {qualifier} shards directory"
            ) from error

    def _write_shard_manifest_at(
        self,
        run_fd: int,
        run_id: str,
        index: int,
        document: Mapping[str, object],
    ) -> Path:
        name = self._shard_name(index)
        if (
            document.get("kind") != "ads-benchmark-plan-shard"
            or type(document.get("shard_index")) is not int
            or document.get("shard_index") != index
        ):
            raise StorageError("shard document identity does not match its index")
        shards_fd = self._open_shards_at(run_fd, run_id, create=True)
        try:
            self._exclusive_json_at(shards_fd, name, document)
        finally:
            os.close(shards_fd)
        return self._run_path(run_id) / "shards" / name

    def write_shard_manifest(
        self,
        run_id: str,
        index: int,
        document: Mapping[str, object],
    ) -> Path:
        """Atomically write one indexed shard below an owned parent plan."""

        with self._shard_manifest_lock(run_id, exclusive=True) as run_fd:
            return self._write_shard_manifest_at(
                run_fd, run_id, index, document
            )

    def write_shard_manifests(
        self,
        run_id: str,
        documents: Iterable[Mapping[str, object]],
    ) -> tuple[Path, ...]:
        """Publish a complete shard set under one exclusive namespace lock."""

        prepared: list[tuple[int, Mapping[str, object]]] = []
        seen: set[int] = set()
        for position, document in enumerate(documents):
            index = document.get("shard_index")
            if type(index) is not int:
                raise StorageError(
                    f"shard document {position} has an invalid shard_index"
                )
            self._shard_name(index)
            if index in seen:
                raise StorageError(f"duplicate shard index in write set: {index}")
            seen.add(index)
            if document.get("kind") != "ads-benchmark-plan-shard":
                raise StorageError(
                    f"shard document {position} has an invalid kind"
                )
            # Reject serialization failures before publishing the first member.
            self._json_text(self._shard_name(index), document)
            prepared.append((index, document))
        if not prepared:
            raise StorageError("cannot write an empty shard manifest set")

        with self._shard_manifest_lock(run_id, exclusive=True) as run_fd:
            return tuple(
                self._write_shard_manifest_at(
                    run_fd, run_id, index, document
                )
                for index, document in prepared
            )

    def read_shard_manifest(self, run_id: str, index: int) -> dict[str, object]:
        """Read one bounded, regular shard document from an owned plan."""

        name = self._shard_name(index)
        with self._shard_manifest_lock(run_id, exclusive=False) as run_fd:
            shards_fd = self._open_shards_at(run_fd, run_id, create=False)
            try:
                document = self._read_json_at(
                    shards_fd,
                    name,
                    f"run {run_id} shard {index}",
                )
            finally:
                os.close(shards_fd)
        assert document is not None
        return document

    def read_all_shard_manifests(
        self, run_id: str
    ) -> tuple[dict[str, object], ...]:
        """Read every deterministically named shard from an owned plan."""

        with self._shard_manifest_lock(run_id, exclusive=False) as run_fd:
            shards_fd = self._open_shards_at(run_fd, run_id, create=False)
            try:
                names = sorted(
                    name
                    for name in os.listdir(shards_fd)
                    if _SHARD_FILE.fullmatch(name)
                )
                if not names:
                    raise StorageError(f"run {run_id} contains no shard manifests")
                documents = []
                for name in names:
                    document = self._read_json_at(
                        shards_fd,
                        name,
                        f"run {run_id} shard manifest {name}",
                    )
                    assert document is not None
                    documents.append(document)
                return tuple(documents)
            finally:
                os.close(shards_fd)

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

    def seal_generated_artifact(
        self, case_directory: Path, name: str
    ) -> Path:
        """Validate and atomically republish one subprocess-created artifact.

        Benchmark executables necessarily create their field CSV directly.
        Once the child has exited, the engine bounds and decodes that file,
        then replaces it through the same fsync-and-rename path used for every
        other persisted result artifact.
        """

        content = self.read_generated_artifact(case_directory, name)
        _, case_id = self._case_coordinates(case_directory)
        try:
            with self._case_fd(case_directory) as case_fd:
                self._atomic_text_at(case_fd, name, content)
        except StorageError:
            raise
        except (OSError, UnicodeError) as error:
            raise StorageError(
                f"cannot seal case {case_id} {name}: {error}"
            ) from error
        return case_directory / name

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
