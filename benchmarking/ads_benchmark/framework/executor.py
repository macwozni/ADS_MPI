"""Generic subprocess execution with timeout, logs, and status transitions."""

from __future__ import annotations

from collections.abc import Mapping
from datetime import datetime, timezone
import math
import os
from pathlib import Path
import signal
import subprocess
import time

from .errors import ExecutionError
from .model import ExecutionContext, PlannedCase
from .registry import Catalog
from .storage import ResultStore
from .validation import MAX_TIMEOUT_SECONDS, decimal_value


RESULT_SCHEMA_VERSION = 1


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
            process.wait()

    def execute(
        self,
        run_id: str,
        case: PlannedCase,
        *,
        environment: Mapping[str, str] | None = None,
    ) -> dict[str, object]:
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
            raise ExecutionError("timeout_seconds is outside runtime range") from error
        if not math.isfinite(timeout) or timeout <= 0:
            raise ExecutionError("timeout_seconds is outside runtime range")
        case_directory = self.store.create_case_directory(run_id, case.case_id)
        context = ExecutionContext(
            repository_root=self.repository_root,
            case_directory=case_directory,
        )
        try:
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
            process = subprocess.Popen(
                command,
                # Every payload gets an isolated, owned working directory.
                # Relative solver outputs therefore stay with this case
                # instead of leaking into the repository checkout.
                cwd=case_directory,
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
            except subprocess.TimeoutExpired:
                self._terminate(process)
                # communicate() retains already-read pipe data after a timeout,
                # so retry once after termination instead of concatenating a
                # potentially byte-typed TimeoutExpired payload.
                stdout, stderr = process.communicate()
                return_code = process.returncode
                state = "timeout"
                error_message = f"process exceeded timeout of {timeout:g} seconds"
        except Exception as error:
            if process is not None:
                self._terminate(process)
                try:
                    stdout, stderr = process.communicate()
                except Exception:
                    stdout = stdout or ""
                    stderr = stderr or ""
                return_code = process.returncode
            error_message = f"process execution failed: {error}"
        except BaseException:
            if process is not None:
                self._terminate(process)
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
                state = "passed"
            except Exception as error:
                error_message = f"result parser failed: {error}"
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
