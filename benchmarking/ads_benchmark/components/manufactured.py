"""Executable adapters for the shared manufactured transient harness."""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import json
import math
from typing import Any, Mapping, Sequence

from ..framework.errors import ExecutionError
from ..framework.model import CaseSpec, ExecutionContext


RESULT_PREFIX = "ADS_BENCHMARK_RESULT "
RESULT_SCHEMA_VERSION = 1
RESULT_KIND = "ads-manufactured-transient-result"
SUPPORTED_EXACT_CASES = {
    "temporal-polynomial": 3,
    "spatial-cosine": 1,
}
_SUPPORTED_SCHEMES = frozenset({"dg", "pr", "be"})
_RESULT_KEYS = {
    "schema_version",
    "kind",
    "exact_case",
    "problem",
    "scheme",
    "requested_final_time",
    "actual_final_time",
    "time_step",
    "steps",
    "initial_l2_error",
    "initial_linf_error",
    "l2_error",
    "linf_error",
    "solution_l2_norm",
    "field_checksum",
    "sample_points_per_axis",
    "field_samples_written",
    "physical_step_wall_seconds",
    "solver_status",
}
_NONNEGATIVE_REAL_FIELDS = (
    "initial_l2_error",
    "initial_linf_error",
    "l2_error",
    "linf_error",
    "solution_l2_norm",
    "physical_step_wall_seconds",
)


def _reject_constant(value: str) -> None:
    raise ExecutionError(f"non-finite JSON number is not allowed: {value}")


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ExecutionError(f"duplicate result key: {key}")
        result[key] = value
    return result


def _binary64_value(value: str, field: str) -> float:
    try:
        parsed = float(Fraction(value))
    except (OverflowError, TypeError, ValueError, ZeroDivisionError) as error:
        raise ExecutionError(f"{field} cannot be represented by the solver") from error
    if not math.isfinite(parsed) or parsed <= 0.0:
        raise ExecutionError(f"{field} cannot be represented by the solver")
    return parsed


def _binary64_argument(value: str, field: str) -> str:
    """Convert an exact model value deliberately to the solver precision."""

    parsed = _binary64_value(value, field)
    if "/" not in value:
        return value
    return format(parsed, ".17g")


def _integer(document: Mapping[str, Any], field: str, *, minimum: int) -> int:
    value = document[field]
    if type(value) is not int or value < minimum:
        qualifier = "nonnegative" if minimum == 0 else "positive"
        raise ExecutionError(f"result {field} must be a {qualifier} integer")
    return value


def _real(
    document: Mapping[str, Any], field: str, *, strictly_positive: bool = False
) -> float:
    value = document[field]
    if isinstance(value, bool) or not isinstance(value, (int, float, Decimal)):
        raise ExecutionError(f"result {field} must be a finite number")
    try:
        normalized = float(value)
    except (OverflowError, ValueError) as error:
        raise ExecutionError(f"result {field} must be a finite number") from error
    if not math.isfinite(normalized):
        raise ExecutionError(f"result {field} must be a finite number")
    if strictly_positive and normalized <= 0.0:
        raise ExecutionError(f"result {field} must be positive")
    return normalized


@dataclass(frozen=True)
class ManufacturedTransientAdapter:
    """Thin Python side of one shared Fortran manufactured adapter."""

    name: str
    execution_ready: bool = True

    def validate_case(self, case: CaseSpec) -> None:
        if case.problem != self.name:
            raise ValueError(f"adapter {self.name} received problem {case.problem}")
        if case.exact_case not in SUPPORTED_EXACT_CASES:
            raise ValueError(
                f"adapter {self.name} does not support exact case "
                f"{case.exact_case}"
            )
        required_case = {
            "temporal": "temporal-polynomial",
            "h": "spatial-cosine",
            "p": "spatial-cosine",
            "validation": "spatial-cosine",
        }.get(case.family)
        if required_case is not None and case.exact_case != required_case:
            raise ValueError(
                f"{case.family} family requires exact case {required_case}"
            )
        if case.scheme not in _SUPPORTED_SCHEMES:
            raise ValueError(f"adapter {self.name} does not support {case.scheme}")
        minimum_degree = SUPPORTED_EXACT_CASES[case.exact_case]
        if any(degree < minimum_degree for degree in case.trial_degree):
            raise ValueError(
                f"adapter {self.name} requires trial degree at least "
                f"{minimum_degree} on every axis for {case.exact_case}"
            )

    def build_payload_command(
        self, case: CaseSpec, context: ExecutionContext
    ) -> Sequence[str]:
        self.validate_case(case)
        executable = (
            context.repository_root
            / "benchmarking"
            / "build"
            / case.build_profile
            / "EXEC"
            / f"{self.name}_manufactured"
        )
        payload = (
            str(executable),
            case.scheme,
            _binary64_argument(case.time.final_time, "final_time"),
            str(case.time.steps),
            *(str(value) for value in case.mesh),
            *(str(value) for value in case.test_degree),
            *(str(value) for value in case.trial_degree),
            *(str(value) for value in case.mpi.process_grid),
            str(case.sampling.points_per_axis),
            "1" if case.sampling.write_samples else "0",
        )
        # Stage 3 stored this exact temporal argv, so it must remain verifiable.
        if case.exact_case == "temporal-polynomial":
            return payload
        return (*payload, case.exact_case)

    def parse_result(self, stdout: str, stderr: str) -> Mapping[str, object]:
        stderr_records = [
            line for line in stderr.splitlines() if line.startswith(RESULT_PREFIX)
        ]
        if stderr_records:
            raise ExecutionError("tagged benchmark result must be written to stdout")
        records = [
            line for line in stdout.splitlines() if line.startswith(RESULT_PREFIX)
        ]
        if len(records) != 1:
            raise ExecutionError(
                "expected exactly one tagged manufactured benchmark result"
            )
        payload = records[0][len(RESULT_PREFIX) :]
        try:
            document = json.loads(
                payload,
                parse_float=Decimal,
                parse_constant=_reject_constant,
                object_pairs_hook=_unique_object,
            )
        except ExecutionError:
            raise
        except (json.JSONDecodeError, TypeError, ValueError) as error:
            raise ExecutionError(f"invalid tagged benchmark JSON: {error}") from error
        if not isinstance(document, dict):
            raise ExecutionError("tagged benchmark result must be a JSON object")

        missing = sorted(_RESULT_KEYS - set(document))
        extra = sorted(set(document) - _RESULT_KEYS)
        if missing or extra:
            details: list[str] = []
            if missing:
                details.append("missing " + ", ".join(missing))
            if extra:
                details.append("unknown " + ", ".join(extra))
            raise ExecutionError(
                "tagged benchmark result has invalid keys: " + "; ".join(details)
            )

        if (
            type(document["schema_version"]) is not int
            or document["schema_version"] != RESULT_SCHEMA_VERSION
        ):
            raise ExecutionError(
                f"result schema_version must be {RESULT_SCHEMA_VERSION}"
            )
        expected_strings = {"kind": RESULT_KIND, "problem": self.name}
        for field, expected in expected_strings.items():
            if document[field] != expected or not isinstance(document[field], str):
                raise ExecutionError(f"result {field} must be {expected!r}")
        if (
            not isinstance(document["exact_case"], str)
            or document["exact_case"] not in SUPPORTED_EXACT_CASES
        ):
            raise ExecutionError(
                "result exact_case must be a registered manufactured case"
            )
        if (
            not isinstance(document["scheme"], str)
            or document["scheme"] not in _SUPPORTED_SCHEMES
        ):
            raise ExecutionError("result scheme must be one of be, dg, pr")

        requested_time = _real(
            document, "requested_final_time", strictly_positive=True
        )
        actual_time = _real(document, "actual_final_time", strictly_positive=True)
        time_step = _real(document, "time_step", strictly_positive=True)
        steps = _integer(document, "steps", minimum=1)
        normalized: dict[str, object] = dict(document)
        normalized.update(
            {
                "requested_final_time": requested_time,
                "actual_final_time": actual_time,
                "time_step": time_step,
                "steps": steps,
            }
        )
        for field in _NONNEGATIVE_REAL_FIELDS:
            value = _real(document, field)
            if value < 0.0:
                raise ExecutionError(f"result {field} must be nonnegative")
            normalized[field] = value
        normalized["field_checksum"] = _real(document, "field_checksum")

        sample_points = _integer(document, "sample_points_per_axis", minimum=1)
        if sample_points < 2 or sample_points > 257:
            raise ExecutionError(
                "result sample_points_per_axis must be between 2 and 257"
            )
        normalized["sample_points_per_axis"] = sample_points
        if type(document["field_samples_written"]) is not bool:
            raise ExecutionError("result field_samples_written must be a boolean")
        solver_status = _integer(document, "solver_status", minimum=0)
        if solver_status != 0:
            raise ExecutionError("result solver_status must be zero")
        normalized["solver_status"] = solver_status

        tolerance = 128.0 * math.ulp(1.0)
        if not math.isclose(
            time_step * steps,
            requested_time,
            rel_tol=tolerance,
            abs_tol=tolerance * max(1.0, requested_time),
        ):
            raise ExecutionError("result time_step*steps does not equal final time")
        if not math.isclose(
            actual_time,
            requested_time,
            rel_tol=tolerance,
            abs_tol=tolerance * max(1.0, requested_time),
        ):
            raise ExecutionError(
                "result actual_final_time does not equal requested time"
            )
        return normalized

    def validate_result(
        self, case: CaseSpec, result: Mapping[str, object]
    ) -> None:
        """Bind a valid manufactured record to the exact planned case."""

        self.validate_case(case)
        expected_strings = {
            "problem": case.problem,
            "scheme": case.scheme,
            "exact_case": case.exact_case,
        }
        for field, expected in expected_strings.items():
            if not isinstance(result.get(field), str) or result[field] != expected:
                raise ExecutionError(
                    f"result {field} does not match the planned case"
                )

        exact_integer_fields = {
            "steps": case.time.steps,
            "sample_points_per_axis": case.sampling.points_per_axis,
        }
        for field, expected in exact_integer_fields.items():
            if type(result.get(field)) is not int or result[field] != expected:
                raise ExecutionError(
                    f"result {field} does not match the planned case"
                )
        if (
            type(result.get("field_samples_written")) is not bool
            or result["field_samples_written"] is not case.sampling.write_samples
        ):
            raise ExecutionError(
                "result field_samples_written does not match the planned case"
            )

        expected_final_time = _binary64_value(
            case.time.final_time, "final_time"
        )
        expected_time_step = _binary64_value(case.time.time_step, "time_step")
        for field, expected in (
            ("requested_final_time", expected_final_time),
            ("actual_final_time", expected_final_time),
            ("time_step", expected_time_step),
        ):
            actual = _real(result, field, strictly_positive=True)
            tolerance = 8.0 * max(math.ulp(actual), math.ulp(expected))
            if abs(actual - expected) > tolerance:
                raise ExecutionError(
                    f"result {field} does not match the planned case"
                )
