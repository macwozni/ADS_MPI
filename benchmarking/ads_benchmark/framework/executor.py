"""Generic subprocess execution with timeout, logs, and status transitions."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import signal
import shutil
import subprocess
import time

from .errors import ExecutionError, StorageError
from .model import ExecutionContext, PlannedCase
from .registry import Catalog
from .storage import ResultStore
from .validation import MAX_TIMEOUT_SECONDS, decimal_value


RESULT_SCHEMA_VERSION = 1
CASE_STATES = frozenset({"planned", "running", "passed", "failed", "timeout"})


def _timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def _validated_argv(values: object, owner: str) -> tuple[str, ...]:
    if isinstance(values, (str, bytes)):
        raise ExecutionError(f"{owner} must return an argv sequence, not a string")
    try:
        arguments = tuple(values)  # type: ignore[arg-type]
    except TypeError as error:
        raise ExecutionError(f"{owner} returned a non-iterable command") from error
    if not arguments:
        raise ExecutionError(f"{owner} returned an empty command")
    if any(
        not isinstance(argument, str) or not argument or "\0" in argument
        for argument in arguments
    ):
        raise ExecutionError(
            f"{owner} command entries must be nonempty NUL-free strings"
        )
    return arguments


def _merge_environment(
    destination: dict[str, str], values: object, owner: str
) -> None:
    if not isinstance(values, Mapping):
        raise ExecutionError(f"{owner} environment must be a mapping")
    for key, value in values.items():
        if (
            not isinstance(key, str)
            or not key
            or "=" in key
            or "\0" in key
            or not isinstance(value, str)
            or "\0" in value
        ):
            raise ExecutionError(
                f"{owner} environment requires valid NUL-free string keys and values"
            )
        destination[key] = value


def _strict_json(value: object) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )


def _timeout_text(value: str | bytes | None) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value


@dataclass(frozen=True)
class ExecutionOutcome:
    case: PlannedCase
    status: Mapping[str, object]
    skipped: bool


@dataclass(frozen=True)
class ExecutionSummary:
    outcomes: tuple[ExecutionOutcome, ...]

    @property
    def passed(self) -> int:
        return sum(
            outcome.status.get("state") == "passed" and not outcome.skipped
            for outcome in self.outcomes
        )

    @property
    def skipped(self) -> int:
        return sum(outcome.skipped for outcome in self.outcomes)

    @property
    def failed(self) -> int:
        return sum(
            outcome.status.get("state") != "passed"
            for outcome in self.outcomes
        )


class Executor:
    """Execute any registered adapter without branches on adapter names."""

    def __init__(
        self,
        catalog: Catalog,
        store: ResultStore,
        repository_root: Path,
    ) -> None:
        self.catalog = catalog
        self.store = store
        self.repository_root = repository_root.resolve()

    @staticmethod
    def _require_supported_measurement(case: PlannedCase) -> None:
        measurement = case.spec.measurement
        if measurement.warmups != 0 or measurement.samples != 1:
            raise ExecutionError(
                "execution currently supports only warmups=0 and samples=1; "
                f"case requests warmups={measurement.warmups} and "
                f"samples={measurement.samples}"
            )

    def _current_command(
        self, case: PlannedCase, case_directory: Path
    ) -> tuple[str, ...]:
        adapter = self.catalog.adapters.get(case.spec.problem)
        payload = _validated_argv(
            adapter.build_payload_command(
                case.spec,
                ExecutionContext(
                    repository_root=self.repository_root,
                    case_directory=case_directory,
                ),
            ),
            f"adapter {adapter.name}",
        )
        launcher = self.catalog.launchers.get(case.spec.launcher)
        return _validated_argv(
            launcher.command(payload, case.spec),
            f"launcher {launcher.name}",
        )

    @staticmethod
    def _require_executable(
        executable: str,
        *,
        working_directory: Path,
        environment: Mapping[str, str],
        owner: str,
    ) -> None:
        if os.sep in executable or (os.altsep and os.altsep in executable):
            candidate = Path(executable)
            if not candidate.is_absolute():
                candidate = working_directory / candidate
            available = candidate.is_file() and os.access(candidate, os.X_OK)
            resolved = str(candidate)
        else:
            located = shutil.which(executable, path=environment.get("PATH"))
            available = located is not None
            resolved = executable
        if not available:
            raise ExecutionError(
                f"execution preflight cannot find executable for {owner}: {resolved}"
            )

    def preflight(
        self,
        cases: Sequence[PlannedCase],
        *,
        environment: Mapping[str, str] | None = None,
    ) -> None:
        """Validate commands and executable availability without result writes."""

        checked: set[tuple[str, str, str]] = set()
        for case in cases:
            context = ExecutionContext(
                repository_root=self.repository_root,
                case_directory=(
                    self.store.results_root
                    / ".preflight"
                    / case.case_id
                ),
            )
            try:
                self._require_supported_measurement(case)
                adapter = self.catalog.adapters.get(case.spec.problem)
                if not adapter.execution_ready:
                    raise ExecutionError(
                        f"adapter {adapter.name} is not execution-ready"
                    )
                adapter.validate_case(case.spec)
                launcher = self.catalog.launchers.get(case.spec.launcher)
                launcher.validate_case(case.spec)
                payload = _validated_argv(
                    adapter.build_payload_command(case.spec, context),
                    f"adapter {adapter.name}",
                )
                command = _validated_argv(
                    launcher.command(payload, case.spec),
                    f"launcher {launcher.name}",
                )
                runtime_environment = os.environ.copy()
                if environment:
                    _merge_environment(runtime_environment, environment, "caller")
                _merge_environment(
                    runtime_environment,
                    launcher.environment(case.spec),
                    f"launcher {launcher.name}",
                )
                for owner, executable in (
                    (f"adapter {adapter.name}", payload[0]),
                    (f"launcher {launcher.name}", command[0]),
                ):
                    key = (owner, executable, runtime_environment.get("PATH", ""))
                    if key in checked:
                        continue
                    self._require_executable(
                        executable,
                        working_directory=context.case_directory,
                        environment=runtime_environment,
                        owner=owner,
                    )
                    checked.add(key)
            except ExecutionError as error:
                raise ExecutionError(
                    f"execution preflight failed for {case.case_id}: {error}"
                ) from error
            except Exception as error:
                raise ExecutionError(
                    f"execution preflight failed for {case.case_id}: {error}"
                ) from error

    def _planned_status(self, case: PlannedCase) -> dict[str, object]:
        return {
            "schema_version": RESULT_SCHEMA_VERSION,
            "case_id": case.case_id,
            "state": "planned",
            "planned_at": _timestamp(),
        }

    def _prepare_case(
        self, run_id: str, case: PlannedCase, *, allow_existing: bool
    ) -> Path:
        case_directory, created = self.store.ensure_case_directory(
            run_id, case.case_id
        )
        if not created and not allow_existing:
            raise StorageError(
                f"case already exists; refusing overwrite: {case.case_id}"
            )
        if created:
            self.store.write_status(case_directory, self._planned_status(case))
        return case_directory

    def _verified_completion(
        self, case_directory: Path, case: PlannedCase
    ) -> tuple[Mapping[str, object] | None, Mapping[str, object] | None]:
        """Accept only a complete result that still satisfies its adapter contract."""

        try:
            status = self.store.read_case_json(
                case_directory, "status.json", missing_ok=True
            )
            result = self.store.read_case_json(
                case_directory, "result.json", missing_ok=True
            )
        except StorageError:
            return None, None
        if status is None or result is None:
            return None, status
        if set(status) != {
            "schema_version",
            "case_id",
            "state",
            "started_at",
            "finished_at",
            "duration_seconds",
            "return_code",
            "command",
        } or (
            type(status.get("schema_version")) is not int
            or status.get("schema_version") != RESULT_SCHEMA_VERSION
            or status.get("case_id") != case.case_id
            or status.get("state") != "passed"
            or type(status.get("return_code")) is not int
            or status.get("return_code") != 0
        ):
            return None, status
        if any(
            not isinstance(status.get(field), str) or not status.get(field)
            for field in ("started_at", "finished_at")
        ):
            return None, status
        duration = status.get("duration_seconds")
        if (
            isinstance(duration, bool)
            or not isinstance(duration, (int, float))
            or not math.isfinite(duration)
            or duration < 0
        ):
            return None, status
        command = status.get("command")
        if (
            not isinstance(command, list)
            or not command
            or any(
                not isinstance(argument, str) or not argument or "\0" in argument
                for argument in command
            )
        ):
            return None, status
        try:
            self._require_supported_measurement(case)
            adapter = self.catalog.adapters.get(case.spec.problem)
            expected_payload = _validated_argv(
                adapter.build_payload_command(
                    case.spec,
                    ExecutionContext(
                        repository_root=self.repository_root,
                        case_directory=case_directory,
                    ),
                ),
                f"adapter {adapter.name}",
            )
        except Exception:
            return None, status
        # Launcher executables and flags may be selected from the environment
        # (for example MPIEXEC), so offline analysis cannot reconstruct the
        # complete historical argv reliably.  The adapter payload is
        # deterministic from the frozen case and must remain its exact suffix.
        if (
            len(command) < len(expected_payload)
            or tuple(command[-len(expected_payload) :]) != expected_payload
        ):
            return None, status
        try:
            stdout = self.store.read_case_log(case_directory, "stdout.log")
            stderr = self.store.read_case_log(case_directory, "stderr.log")
        except StorageError:
            return None, status
        if set(result) != {
            "schema_version",
            "kind",
            "case_id",
            "status",
            "configuration",
            "timing",
            "domain_result",
        }:
            return None, status
        if (
            type(result.get("schema_version")) is not int
            or result.get("schema_version") != RESULT_SCHEMA_VERSION
            or result.get("kind") != "ads-benchmark-case-result"
            or result.get("case_id") != case.case_id
            or result.get("status") != "passed"
            or result.get("configuration") != case.spec.to_dict()
        ):
            return None, status
        timing = result.get("timing")
        if not isinstance(timing, dict) or set(timing) != {"wall_seconds"}:
            return None, status
        wall_seconds = timing.get("wall_seconds")
        if (
            isinstance(wall_seconds, bool)
            or not isinstance(wall_seconds, (int, float))
            or not math.isfinite(wall_seconds)
            or wall_seconds < 0
            or wall_seconds != duration
        ):
            return None, status
        domain_result = result.get("domain_result")
        if not isinstance(domain_result, Mapping):
            return None, status
        try:
            reparsed = adapter.parse_result(stdout, stderr)
            if not isinstance(reparsed, Mapping):
                return None, status
            adapter.validate_result(case.spec, reparsed)
            if _strict_json(dict(reparsed)) != _strict_json(dict(domain_result)):
                return None, status
            if (
                case.spec.sampling.write_samples
                and not self.store.case_file_is_regular(
                    case_directory, "field_samples.csv"
                )
            ):
                return None, status
        except Exception:
            return None, status
        return result, status

    def verified_result(
        self, case_directory: Path, case: PlannedCase
    ) -> Mapping[str, object] | None:
        """Return the fully verified persisted result, or ``None`` for retry."""

        result, _ = self._verified_completion(case_directory, case)
        return result

    def execute_frozen(
        self,
        run_id: str,
        cases: Sequence[PlannedCase],
        *,
        resume: bool = False,
        environment: Mapping[str, str] | None = None,
        observer: Callable[[int, int, ExecutionOutcome], None] | None = None,
        preflight: bool = True,
    ) -> ExecutionSummary:
        """Execute a manifest's cases, with verified result reuse on resume."""

        if preflight:
            self.preflight(cases, environment=environment)
        outcomes: list[ExecutionOutcome] = []
        with self.store.execution_lock(run_id):
            verified_cases: dict[
                str,
                tuple[
                    Mapping[str, object] | None,
                    Mapping[str, object] | None,
                ],
            ] = {}
            if resume:
                # This first path pass is intentionally read-only.  Missing
                # case directories must remain missing if launcher
                # compatibility later refuses the resume.
                directories = {
                    case.case_id: self.store.case_directory_path(
                        run_id, case.case_id
                    )
                    for case in cases
                }
                # Validate every reusable case before retrying any incomplete
                # one.  Payload-only verification remains sufficient for
                # offline analysis, but one resumed run must never mix launcher
                # prefixes selected by different configurations/environments.
                for case in cases:
                    case_directory = directories[case.case_id]
                    completion = self._verified_completion(case_directory, case)
                    verified_cases[case.case_id] = completion
                    verified, previous_status = completion
                    if verified is None:
                        continue
                    assert previous_status is not None
                    try:
                        current_command = self._current_command(
                            case, case_directory
                        )
                    except Exception as error:
                        raise ExecutionError(
                            "cannot reconstruct current resume command for "
                            f"{case.case_id}: {error}"
                        ) from error
                    if previous_status.get("command") != list(current_command):
                        raise ExecutionError(
                            "resume launcher command mismatch for "
                            f"{case.case_id}; refusing to mix launcher "
                            "configurations in one run"
                        )
                directories = {
                    case.case_id: self._prepare_case(
                        run_id, case, allow_existing=True
                    )
                    for case in cases
                }
            else:
                directories = {
                    case.case_id: self._prepare_case(
                        run_id, case, allow_existing=False
                    )
                    for case in cases
                }
            total = len(cases)
            for position, case in enumerate(cases, start=1):
                case_directory = directories[case.case_id]
                if resume:
                    verified, previous_status = verified_cases[case.case_id]
                else:
                    verified, previous_status = self._verified_completion(
                        case_directory, case
                    )
                if resume and verified is not None:
                    assert previous_status is not None
                    outcome = ExecutionOutcome(
                        case=case, status=previous_status, skipped=True
                    )
                else:
                    if resume:
                        self.store.write_status(
                            case_directory, self._planned_status(case)
                        )
                    try:
                        status = self.execute(
                            run_id,
                            case,
                            environment=environment,
                            reuse_case_directory=True,
                        )
                    except ExecutionError:
                        recorded = self.store.read_case_json(
                            case_directory, "status.json"
                        )
                        assert recorded is not None
                        status = recorded
                    outcome = ExecutionOutcome(
                        case=case, status=status, skipped=False
                    )
                outcomes.append(outcome)
                if observer is not None:
                    observer(position, total, outcome)
        return ExecutionSummary(tuple(outcomes))

    @staticmethod
    def _terminate(process: subprocess.Popen[str]) -> None:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass

        # Reap the direct child while allowing every process in its new
        # session a short grace period.  A child may outlive its parent and
        # keep stdout/stderr pipes open, so checking only process.wait() is not
        # sufficient.
        deadline = time.monotonic() + 0.25
        while time.monotonic() < deadline:
            process.poll()
            try:
                os.killpg(process.pid, 0)
            except ProcessLookupError:
                break
            time.sleep(0.01)
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        try:
            process.wait(timeout=1)
        except subprocess.TimeoutExpired:
            process.kill()
            try:
                process.wait(timeout=1)
            except subprocess.TimeoutExpired:
                pass

    @staticmethod
    def _close_pipes(process: subprocess.Popen[str]) -> None:
        for stream in (process.stdout, process.stderr):
            if stream is not None:
                try:
                    stream.close()
                except OSError:
                    pass

    def execute(
        self,
        run_id: str,
        case: PlannedCase,
        *,
        environment: Mapping[str, str] | None = None,
        reuse_case_directory: bool = False,
    ) -> dict[str, object]:
        case_directory = self._prepare_case(
            run_id, case, allow_existing=reuse_case_directory
        )
        context = ExecutionContext(
            repository_root=self.repository_root,
            case_directory=case_directory,
        )
        try:
            self._require_supported_measurement(case)
            adapter = self.catalog.adapters.get(case.spec.problem)
            if not adapter.execution_ready:
                raise ExecutionError(
                    f"adapter {adapter.name} is planning-only in this implementation stage"
                )
            launcher = self.catalog.launchers.get(case.spec.launcher)
            timeout_decimal = decimal_value(
                case.spec.measurement.timeout_seconds, "timeout_seconds"
            )
            if timeout_decimal <= 0 or timeout_decimal > MAX_TIMEOUT_SECONDS:
                raise ExecutionError("timeout_seconds is outside runtime range")
            try:
                timeout = float(timeout_decimal)
            except (OverflowError, ValueError) as error:
                raise ExecutionError(
                    "timeout_seconds is outside runtime range"
                ) from error
            if not math.isfinite(timeout) or timeout <= 0:
                raise ExecutionError("timeout_seconds is outside runtime range")
            self.store.remove_result(case_directory)
            self.store.remove_generated_artifact(
                case_directory, "field_samples.csv"
            )
            self.store.write_log(case_directory, "stdout.log", "")
            self.store.write_log(case_directory, "stderr.log", "")
            payload = _validated_argv(
                adapter.build_payload_command(case.spec, context),
                f"adapter {adapter.name}",
            )
            command = _validated_argv(
                launcher.command(payload, case.spec), f"launcher {launcher.name}"
            )
            runtime_environment = os.environ.copy()
            if environment:
                _merge_environment(runtime_environment, environment, "caller")
            _merge_environment(
                runtime_environment,
                launcher.environment(case.spec),
                f"launcher {launcher.name}",
            )
            # Case semantics win over ambient/caller values so a recorded plan
            # cannot silently run with a different OpenMP layout.
            runtime_environment.update(
                {
                    "OMP_NUM_THREADS": str(case.spec.openmp_threads),
                    "OMP_DYNAMIC": "FALSE",
                    "OMP_PROC_BIND": "close",
                    "OMP_PLACES": "cores",
                }
            )
        except Exception as error:
            # Construction may fail before stale artifacts were cleared.
            self.store.remove_result(case_directory)
            self.store.remove_generated_artifact(
                case_directory, "field_samples.csv"
            )
            self.store.write_log(case_directory, "stdout.log", "")
            self.store.write_log(case_directory, "stderr.log", "")
            self.store.write_status(
                case_directory,
                {
                    "schema_version": RESULT_SCHEMA_VERSION,
                    "case_id": case.case_id,
                    "state": "failed",
                    "finished_at": _timestamp(),
                    "error": f"command construction failed: {error}",
                },
            )
            if isinstance(error, ExecutionError):
                raise
            raise ExecutionError(f"cannot construct command: {error}") from error

        started_at = _timestamp()
        self.store.write_status(
            case_directory,
            {
                "schema_version": RESULT_SCHEMA_VERSION,
                "case_id": case.case_id,
                "state": "running",
                "started_at": started_at,
                "command": list(command),
            },
        )
        start = time.monotonic()
        stdout = ""
        stderr = ""
        return_code: int | None = None
        state = "failed"
        error_message: str | None = None
        process: subprocess.Popen[str] | None = None
        try:
            with self.store.open_case_directory(case_directory) as case_fd:
                anchored_cwd = f"/proc/self/fd/{case_fd}"
                if not os.path.isdir(anchored_cwd):
                    raise ExecutionError(
                        "safe child cwd requires a mounted /proc/self/fd"
                    )
                process = subprocess.Popen(
                    command,
                    # The inherited descriptor anchors cwd to the verified
                    # inode even if the visible case path is swapped later.
                    cwd=anchored_cwd,
                    pass_fds=(case_fd,),
                    env=runtime_environment,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    start_new_session=True,
                )
                try:
                    stdout, stderr = process.communicate(timeout=timeout)
                    return_code = process.returncode
                except subprocess.TimeoutExpired as error:
                    stdout = _timeout_text(error.stdout)
                    stderr = _timeout_text(error.stderr)
                    self._terminate(process)
                    # A detached descendant may retain the inherited pipe writers.
                    # Closing our readers guarantees timeout handling never waits
                    # indefinitely for EOF from a process outside the owned group.
                    self._close_pipes(process)
                    return_code = process.returncode
                    state = "timeout"
                    error_message = (
                        f"process exceeded timeout of {timeout:g} seconds"
                    )
        except Exception as error:
            if process is not None:
                self._terminate(process)
                self._close_pipes(process)
                return_code = process.returncode
            error_message = f"process execution failed: {error}"
        except BaseException:
            if process is not None:
                self._terminate(process)
                self._close_pipes(process)
            raise

        duration = time.monotonic() - start
        self.store.write_log(case_directory, "stdout.log", stdout)
        self.store.write_log(case_directory, "stderr.log", stderr)

        parsed: Mapping[str, object] | None = None
        if error_message is None and return_code == 0:
            try:
                parsed = adapter.parse_result(stdout, stderr)
                if not isinstance(parsed, Mapping):
                    raise TypeError("adapter result must be a mapping")
            except Exception as error:
                error_message = f"result parser failed: {error}"
            if error_message is None and parsed is not None:
                try:
                    adapter.validate_result(case.spec, parsed)
                    if (
                        case.spec.sampling.write_samples
                        and not self.store.case_file_is_regular(
                            case_directory, "field_samples.csv"
                        )
                    ):
                        raise ExecutionError(
                            "requested field_samples.csv is missing or not regular"
                        )
                    state = "passed"
                except Exception as error:
                    error_message = f"result validation failed: {error}"
        elif error_message is None:
            error_message = f"process exited with status {return_code}"

        finished_at = _timestamp()
        if state == "passed" and parsed is not None:
            result = {
                "schema_version": RESULT_SCHEMA_VERSION,
                "kind": "ads-benchmark-case-result",
                "case_id": case.case_id,
                "status": "passed",
                "configuration": case.spec.to_dict(),
                "timing": {"wall_seconds": duration},
                "domain_result": dict(parsed),
            }
            try:
                self.store.write_result(case_directory, result)
            except Exception as error:
                state = "failed"
                error_message = f"result serialization failed: {error}"

        status: dict[str, object] = {
            "schema_version": RESULT_SCHEMA_VERSION,
            "case_id": case.case_id,
            "state": state,
            "started_at": started_at,
            "finished_at": finished_at,
            "duration_seconds": duration,
            "return_code": return_code,
            "command": list(command),
        }
        if error_message is not None:
            status["error"] = error_message
        self.store.write_status(case_directory, status)
        return status
