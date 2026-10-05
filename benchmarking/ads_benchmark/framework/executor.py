"""Generic subprocess execution with timeout, logs, and status transitions."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
import errno
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
FAILURE_KINDS = frozenset(
    {"numerical", "mpi", "timeout", "resource", "configuration"}
)
RETRYABLE_FAILURE_KINDS = frozenset({"mpi", "timeout", "resource"})
_RESOURCE_ERRNOS = frozenset(
    {errno.ENOMEM, errno.ENOSPC, errno.EMFILE, errno.ENFILE, errno.EAGAIN}
)
_REPEATED_TIMING_KEYS = {
    "wall_seconds",
    "metric",
    "warmup_samples",
    "measured_samples",
    "warmup_process_wall_seconds",
    "measured_process_wall_seconds",
    "minimum_reliable_seconds",
    "reliable",
    "openmp_environment",
}
_OPENMP_ENVIRONMENT_KEYS = {
    "OMP_NUM_THREADS",
    "OMP_DYNAMIC",
    "OMP_PROC_BIND",
    "OMP_PLACES",
}


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


def _is_resource_exhaustion(error: BaseException) -> bool:
    """Recognize bounded OS resource failures through wrapper exceptions."""

    seen: set[int] = set()
    current: BaseException | None = error
    while current is not None and id(current) not in seen:
        seen.add(id(current))
        if isinstance(current, MemoryError):
            return True
        if isinstance(current, OSError) and current.errno in _RESOURCE_ERRNOS:
            return True
        current = current.__cause__ or current.__context__
    return False


def _failure_kind_for_exception(
    error: BaseException, *, default: str
) -> str:
    if default not in FAILURE_KINDS:
        raise ValueError(f"invalid default failure kind: {default}")
    return "resource" if _is_resource_exhaustion(error) else default


def _attempt_record(
    status: Mapping[str, object], attempt: int
) -> dict[str, object]:
    """Keep retry history small; stdout/stderr remain the final attempt logs."""

    state = status.get("state")
    failure_kind = status.get("failure_kind")
    error = status.get("error")
    return_code = status.get("return_code")
    return {
        "attempt": attempt,
        "state": state,
        "failure_kind": failure_kind if state != "passed" else None,
        "error": error if state != "passed" else None,
        "return_code": return_code,
    }


def _valid_attempt_history(
    value: object, final_status: Mapping[str, object]
) -> bool:
    if not isinstance(value, list) or not value:
        return False
    expected_keys = {
        "attempt",
        "state",
        "failure_kind",
        "error",
        "return_code",
    }
    for index, record in enumerate(value, start=1):
        if not isinstance(record, dict) or set(record) != expected_keys:
            return False
        if (
            type(record.get("attempt")) is not int
            or record.get("attempt") != index
        ):
            return False
        state = record.get("state")
        if state not in {"passed", "failed", "timeout"}:
            return False
        return_code = record.get("return_code")
        if return_code is not None and type(return_code) is not int:
            return False
        if state == "passed":
            if (
                record.get("failure_kind") is not None
                or record.get("error") is not None
            ):
                return False
        else:
            if record.get("failure_kind") not in FAILURE_KINDS:
                return False
            if not isinstance(record.get("error"), str) or not record.get("error"):
                return False
        if (
            index < len(value)
            and record.get("failure_kind") not in RETRYABLE_FAILURE_KINDS
        ):
            return False
    final = value[-1]
    return (
        final.get("state") == final_status.get("state")
        and final.get("return_code") == final_status.get("return_code")
    )


def _uses_repetitions(case: PlannedCase) -> bool:
    measurement = case.spec.measurement
    return (
        measurement.warmups != 0
        or measurement.samples != 1
        or measurement.minimum_sample_seconds is not None
        or case.spec.family in {"strong", "weak"}
    )


def _minimum_sample_seconds(case: PlannedCase) -> float:
    text = case.spec.measurement.minimum_sample_seconds
    if text is None:
        return 0.0
    value = decimal_value(text, "minimum_sample_seconds")
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise ExecutionError(
            "minimum_sample_seconds is outside runtime range"
        ) from error
    if not math.isfinite(result) or result <= 0.0:
        raise ExecutionError("minimum_sample_seconds is outside runtime range")
    return result


def _openmp_environment(case: PlannedCase) -> dict[str, str]:
    dynamic = case.spec.openmp_dynamic
    return {
        "OMP_NUM_THREADS": str(case.spec.openmp_threads),
        "OMP_DYNAMIC": "TRUE" if dynamic is True else "FALSE",
        "OMP_PROC_BIND": case.spec.openmp_proc_bind or "close",
        "OMP_PLACES": case.spec.openmp_places or "cores",
    }


def _physical_step_seconds(result: Mapping[str, object]) -> float:
    value = result.get("physical_step_wall_seconds")
    seconds = _nonnegative_seconds(value)
    if seconds is None:
        raise ExecutionError(
            "result physical_step_wall_seconds must be a finite nonnegative number"
        )
    return seconds


def _nonnegative_seconds(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    try:
        seconds = float(value)
    except (OverflowError, ValueError):
        return None
    if not math.isfinite(seconds) or seconds < 0.0:
        return None
    return seconds


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


@dataclass(frozen=True)
class _ProcessAttempt:
    stdout: str
    stderr: str
    return_code: int | None
    duration_seconds: float
    state: str
    error: str | None
    failure_kind: str | None


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
                adapter = self.catalog.adapters.get(case.spec.problem)
                if not adapter.execution_ready:
                    raise ExecutionError(
                        f"adapter {adapter.name} is not execution-ready"
                    )
                adapter.validate_case(case.spec)
                launcher = self.catalog.launchers.get(case.spec.launcher)
                launcher.validate_case(case.spec)
                validate_resources = getattr(launcher, "validate_resources", None)
                if callable(validate_resources):
                    validate_resources(case.spec)
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
        passed_status_keys = {
            "schema_version",
            "case_id",
            "state",
            "started_at",
            "finished_at",
            "duration_seconds",
            "return_code",
            "command",
        }
        if set(status) not in (
            passed_status_keys,
            passed_status_keys | {"attempts"},
        ) or (
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
        duration = _nonnegative_seconds(status.get("duration_seconds"))
        if duration is None:
            return None, status
        if "attempts" in status and not _valid_attempt_history(
            status.get("attempts"), status
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
        historical_payload = tuple(command[-len(expected_payload) :])
        # Detached-ref orchestration deliberately removes its owned worktrees
        # after a successful comparison.  The persisted execution record binds
        # the historical executable by digest, while offline verification can
        # still require the adapter basename and every semantic payload
        # argument without requiring that old absolute build path to exist.
        if (
            len(command) < len(expected_payload)
            or Path(historical_payload[0]).name
            != Path(expected_payload[0]).name
            or historical_payload[1:] != expected_payload[1:]
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
        repeated = _uses_repetitions(case)
        timing = result.get("timing")
        expected_timing_keys = (
            _REPEATED_TIMING_KEYS if repeated else {"wall_seconds"}
        )
        if not isinstance(timing, dict) or set(timing) != expected_timing_keys:
            return None, status
        wall_seconds = _nonnegative_seconds(timing.get("wall_seconds"))
        if wall_seconds is None or wall_seconds != duration:
            return None, status

        measured_samples: list[float] | None = None
        if repeated:
            if timing.get("metric") != "physical_step_wall_seconds":
                return None, status

            def sample_series(field: str, count: int) -> list[float] | None:
                raw_values = timing.get(field)
                if not isinstance(raw_values, list) or len(raw_values) != count:
                    return None
                values: list[float] = []
                for raw_value in raw_values:
                    value = _nonnegative_seconds(raw_value)
                    if value is None:
                        return None
                    values.append(value)
                return values

            warmup_samples = sample_series(
                "warmup_samples", case.spec.measurement.warmups
            )
            measured_samples = sample_series(
                "measured_samples", case.spec.measurement.samples
            )
            warmup_process_wall = sample_series(
                "warmup_process_wall_seconds", case.spec.measurement.warmups
            )
            measured_process_wall = sample_series(
                "measured_process_wall_seconds", case.spec.measurement.samples
            )
            if any(
                values is None
                for values in (
                    warmup_samples,
                    measured_samples,
                    warmup_process_wall,
                    measured_process_wall,
                )
            ):
                return None, status

            minimum = _nonnegative_seconds(
                timing.get("minimum_reliable_seconds")
            )
            try:
                expected_minimum = _minimum_sample_seconds(case)
            except Exception:
                return None, status
            if minimum is None or minimum != expected_minimum:
                return None, status
            assert measured_samples is not None
            expected_reliable = min(measured_samples) >= expected_minimum
            if (
                type(timing.get("reliable")) is not bool
                or timing.get("reliable") is not expected_reliable
            ):
                return None, status
            openmp_environment = timing.get("openmp_environment")
            if (
                not isinstance(openmp_environment, dict)
                or set(openmp_environment) != _OPENMP_ENVIRONMENT_KEYS
                or any(
                    not isinstance(value, str)
                    for value in openmp_environment.values()
                )
                or openmp_environment != _openmp_environment(case)
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
                repeated
                and measured_samples is not None
                and _physical_step_seconds(reparsed) != measured_samples[-1]
            ):
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

    def _execute_with_retries(
        self,
        run_id: str,
        case: PlannedCase,
        case_directory: Path,
        *,
        environment: Mapping[str, str] | None,
        max_retries: int,
    ) -> Mapping[str, object]:
        attempts: list[dict[str, object]] = []
        for attempt in range(1, max_retries + 2):
            try:
                status: Mapping[str, object] = self.execute(
                    run_id,
                    case,
                    environment=environment,
                    reuse_case_directory=True,
                )
            except ExecutionError:
                recorded = self.store.read_case_json(
                    case_directory, "status.json"
                )
                if recorded is None:
                    raise
                status = recorded

            state = status.get("state")
            if state not in {"passed", "failed", "timeout"}:
                raise ExecutionError(
                    f"case {case.case_id} returned invalid terminal state {state!r}"
                )
            failure_kind = status.get("failure_kind")
            if state != "passed" and failure_kind not in FAILURE_KINDS:
                raise ExecutionError(
                    f"case {case.case_id} failure lacks a valid failure_kind"
                )

            attempts.append(_attempt_record(status, attempt))
            final_status = dict(status)
            if max_retries > 0:
                final_status["attempts"] = list(attempts)
                self.store.write_status(case_directory, final_status)

            if state == "passed":
                return final_status
            if (
                failure_kind not in RETRYABLE_FAILURE_KINDS
                or attempt > max_retries
            ):
                return final_status

        raise AssertionError("retry loop exhausted without a terminal status")

    def execute_new_case_locked(
        self,
        run_id: str,
        case: PlannedCase,
        *,
        environment: Mapping[str, str] | None = None,
        max_retries: int = 0,
    ) -> Mapping[str, object]:
        """Execute one new case while the caller owns the run execution lock.

        The A/B orchestrator must hold two run locks at once so that it can
        alternate matched cases.  Keeping this small entry point in the
        executor lets that orchestrator reuse the ordinary process, logging,
        status, parsing, and retry implementation instead of growing a second
        execution engine.
        """

        if type(max_retries) is not int or max_retries < 0:
            raise ExecutionError("max_retries must be a nonnegative integer")
        case_directory = self._prepare_case(
            run_id, case, allow_existing=False
        )
        return self._execute_with_retries(
            run_id,
            case,
            case_directory,
            environment=environment,
            max_retries=max_retries,
        )

    def execute_frozen(
        self,
        run_id: str,
        cases: Sequence[PlannedCase],
        *,
        resume: bool = False,
        environment: Mapping[str, str] | None = None,
        observer: Callable[[int, int, ExecutionOutcome], None] | None = None,
        preflight: bool = True,
        max_retries: int = 0,
    ) -> ExecutionSummary:
        """Execute a manifest's cases, with verified result reuse on resume."""

        if type(max_retries) is not int or max_retries < 0:
            raise ExecutionError("max_retries must be a nonnegative integer")
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
                    status = self._execute_with_retries(
                        run_id,
                        case,
                        case_directory,
                        environment=environment,
                        max_retries=max_retries,
                    )
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

    def _run_process_attempt(
        self,
        command: Sequence[str],
        case_directory: Path,
        runtime_environment: Mapping[str, str],
        timeout: float,
        *,
        nonzero_failure_kind: str,
    ) -> _ProcessAttempt:
        start = time.monotonic()
        stdout = ""
        stderr = ""
        return_code: int | None = None
        state = "failed"
        error_message: str | None = None
        failure_kind: str | None = None
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
                    self._close_pipes(process)
                    return_code = process.returncode
                    state = "timeout"
                    failure_kind = "timeout"
                    error_message = (
                        f"process exceeded timeout of {timeout:g} seconds"
                    )
        except Exception as error:
            if process is not None:
                self._terminate(process)
                self._close_pipes(process)
                return_code = process.returncode
            failure_kind = _failure_kind_for_exception(
                error, default="configuration"
            )
            error_message = f"process execution failed: {error}"
        except BaseException:
            if process is not None:
                self._terminate(process)
                self._close_pipes(process)
            raise

        duration = time.monotonic() - start
        if error_message is None and return_code == 0:
            state = "passed"
        elif error_message is None:
            failure_kind = nonzero_failure_kind
            error_message = f"process exited with status {return_code}"
        return _ProcessAttempt(
            stdout=stdout,
            stderr=stderr,
            return_code=return_code,
            duration_seconds=duration,
            state=state,
            error=error_message,
            failure_kind=failure_kind,
        )

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
            minimum_sample_seconds = _minimum_sample_seconds(case)
            repeated = _uses_repetitions(case)
            openmp_environment = _openmp_environment(case)
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
            runtime_environment.update(openmp_environment)
        except Exception as error:
            # Construction may fail before stale artifacts were cleared.
            failure_kind = _failure_kind_for_exception(
                error, default="configuration"
            )
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
                    "failure_kind": failure_kind,
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
        state = "passed"
        error_message: str | None = None
        failure_kind: str | None = None
        final_parsed: Mapping[str, object] | None = None
        warmup_samples: list[float] = []
        measured_samples: list[float] = []
        warmup_process_wall_seconds: list[float] = []
        measured_process_wall_seconds: list[float] = []

        phases = (
            ("warmup", case.spec.measurement.warmups),
            ("measured", case.spec.measurement.samples),
        )
        for phase, count in phases:
            if error_message is not None:
                break
            for index in range(1, count + 1):
                try:
                    # Every repetition is an independent, identical process.
                    # Removing the artifact first makes validation specific to
                    # this attempt instead of accepting a stale earlier file.
                    self.store.remove_generated_artifact(
                        case_directory, "field_samples.csv"
                    )
                except Exception as error:
                    state = "failed"
                    failure_kind = _failure_kind_for_exception(
                        error, default="configuration"
                    )
                    error_message = (
                        f"{phase} {index}: artifact preparation failed: {error}"
                    )
                    break

                attempt = self._run_process_attempt(
                    command,
                    case_directory,
                    runtime_environment,
                    timeout,
                    nonzero_failure_kind=(
                        "mpi" if launcher.name == "mpi" else "numerical"
                    ),
                )
                stdout = attempt.stdout
                stderr = attempt.stderr
                return_code = attempt.return_code
                if attempt.state != "passed":
                    state = attempt.state
                    failure_kind = attempt.failure_kind
                    error_message = f"{phase} {index}: {attempt.error}"
                    break

                try:
                    parsed = adapter.parse_result(stdout, stderr)
                    if not isinstance(parsed, Mapping):
                        raise TypeError("adapter result must be a mapping")
                except Exception as error:
                    state = "failed"
                    failure_kind = "numerical"
                    error_message = f"result parser failed: {error}"
                    break
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
                    if case.spec.sampling.write_samples:
                        self.store.seal_generated_artifact(
                            case_directory, "field_samples.csv"
                        )
                except Exception as error:
                    state = "failed"
                    failure_kind = "numerical"
                    error_message = f"result validation failed: {error}"
                    break

                if repeated:
                    try:
                        sample_seconds = _physical_step_seconds(parsed)
                    except Exception as error:
                        state = "failed"
                        failure_kind = "numerical"
                        error_message = f"result measurement failed: {error}"
                        break
                    if phase == "warmup":
                        warmup_samples.append(sample_seconds)
                        warmup_process_wall_seconds.append(
                            attempt.duration_seconds
                        )
                    else:
                        measured_samples.append(sample_seconds)
                        measured_process_wall_seconds.append(
                            attempt.duration_seconds
                        )
                if phase == "measured":
                    final_parsed = parsed

            if error_message is not None:
                break

        if state == "passed" and final_parsed is None:
            state = "failed"
            failure_kind = "numerical"
            error_message = "execution produced no measured result"

        duration = time.monotonic() - start
        self.store.write_log(case_directory, "stdout.log", stdout)
        self.store.write_log(case_directory, "stderr.log", stderr)

        finished_at = _timestamp()
        if state == "passed" and final_parsed is not None:
            if repeated:
                timing: dict[str, object] = {
                    "wall_seconds": duration,
                    "metric": "physical_step_wall_seconds",
                    "warmup_samples": warmup_samples,
                    "measured_samples": measured_samples,
                    "warmup_process_wall_seconds": (
                        warmup_process_wall_seconds
                    ),
                    "measured_process_wall_seconds": (
                        measured_process_wall_seconds
                    ),
                    "minimum_reliable_seconds": minimum_sample_seconds,
                    "reliable": (
                        min(measured_samples) >= minimum_sample_seconds
                    ),
                    "openmp_environment": dict(openmp_environment),
                }
            else:
                timing = {"wall_seconds": duration}
            result = {
                "schema_version": RESULT_SCHEMA_VERSION,
                "kind": "ads-benchmark-case-result",
                "case_id": case.case_id,
                "status": "passed",
                "configuration": case.spec.to_dict(),
                "timing": timing,
                "domain_result": dict(final_parsed),
            }
            try:
                self.store.write_result(case_directory, result)
            except Exception as error:
                state = "failed"
                failure_kind = _failure_kind_for_exception(
                    error, default="numerical"
                )
                error_message = f"result serialization failed: {error}"
                try:
                    self.store.remove_result(case_directory)
                except Exception:
                    pass

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
            status["failure_kind"] = failure_kind or "numerical"
            status["error"] = error_message
        self.store.write_status(case_directory, status)
        return status
