"""Shared, strict result loading for strong- and weak-scaling analyses."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
import json
import math

from .base import AnalysisError
from .statistics import SampleStatistics, summarize_samples
from .validation import FieldResultPoint, extract_field_result


REPEATED_TIMING_KEYS = {
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
OPENMP_ENVIRONMENT_KEYS = {
    "OMP_NUM_THREADS",
    "OMP_DYNAMIC",
    "OMP_PROC_BIND",
    "OMP_PLACES",
}


def canonical(value: object) -> str:
    try:
        return json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
    except (TypeError, ValueError) as error:
        raise AnalysisError(f"configuration is not strict JSON: {error}") from error


def mapping(value: object, field: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise AnalysisError(f"{field} must be an object")
    return value


def positive_integer(value: object, field: str) -> int:
    if type(value) is not int or value <= 0:
        raise AnalysisError(f"{field} must be a positive integer")
    return value


def _nonnegative_integer(value: object, field: str) -> int:
    if type(value) is not int or value < 0:
        raise AnalysisError(f"{field} must be a nonnegative integer")
    return value


def _finite(value: object, field: str, *, nonnegative: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, Decimal)):
        raise AnalysisError(f"{field} must be a finite number")
    result = float(value)
    if not math.isfinite(result) or (nonnegative and result < 0.0):
        raise AnalysisError(f"{field} must be finite and nonnegative")
    return result


def _seconds_list(value: object, field: str) -> tuple[float, ...]:
    if not isinstance(value, list):
        raise AnalysisError(f"{field} must be an array")
    return tuple(
        _finite(item, f"{field}[{index}]", nonnegative=True)
        for index, item in enumerate(value)
    )


def _decimal_string(value: object, field: str) -> float:
    if not isinstance(value, str) or not value:
        raise AnalysisError(f"{field} must be a nonempty decimal string")
    try:
        parsed = Decimal(value)
    except (InvalidOperation, ValueError) as error:
        raise AnalysisError(f"{field} must be a decimal string") from error
    if not parsed.is_finite() or parsed <= 0:
        raise AnalysisError(f"{field} must be positive and finite")
    return float(parsed)


@dataclass(frozen=True)
class ScalingPoint:
    field: FieldResultPoint
    warmup_samples: tuple[float, ...]
    measured_samples: tuple[float, ...]
    warmup_process_wall_seconds: tuple[float, ...]
    measured_process_wall_seconds: tuple[float, ...]
    statistics: SampleStatistics
    minimum_reliable_seconds: float
    reliable: bool
    openmp_environment: Mapping[str, str]

    @property
    def resources(self) -> int:
        return self.field.mpi_ranks * self.field.openmp_threads


def extract_scaling_point(
    document: Mapping[str, object], *, expected_family: str
) -> ScalingPoint:
    """Validate common repeated timings, OpenMP policy, and full field data."""

    point = extract_field_result(document, expected_family=expected_family)
    case_id = point.case_id
    configuration = point.configuration
    measurement = mapping(
        configuration.get("measurement"), f"{case_id}.measurement"
    )
    warmup_count = _nonnegative_integer(
        measurement.get("warmups"), f"{case_id}.measurement.warmups"
    )
    sample_count = positive_integer(
        measurement.get("samples"), f"{case_id}.measurement.samples"
    )
    requested_minimum = _decimal_string(
        measurement.get("minimum_sample_seconds"),
        f"{case_id}.measurement.minimum_sample_seconds",
    )

    timing = mapping(document.get("timing"), f"{case_id}.timing")
    if set(timing) != REPEATED_TIMING_KEYS:
        missing = sorted(REPEATED_TIMING_KEYS - set(timing))
        extra = sorted(set(timing) - REPEATED_TIMING_KEYS)
        raise AnalysisError(
            f"{case_id}.timing has invalid keys: missing={missing}, extra={extra}"
        )
    _finite(timing.get("wall_seconds"), f"{case_id}.timing.wall_seconds", nonnegative=True)
    if timing.get("metric") != "physical_step_wall_seconds":
        raise AnalysisError(
            f"{case_id}: timing metric must be physical_step_wall_seconds"
        )
    warmups = _seconds_list(
        timing.get("warmup_samples"), f"{case_id}.timing.warmup_samples"
    )
    measured = _seconds_list(
        timing.get("measured_samples"), f"{case_id}.timing.measured_samples"
    )
    warmup_process = _seconds_list(
        timing.get("warmup_process_wall_seconds"),
        f"{case_id}.timing.warmup_process_wall_seconds",
    )
    measured_process = _seconds_list(
        timing.get("measured_process_wall_seconds"),
        f"{case_id}.timing.measured_process_wall_seconds",
    )
    if len(warmups) != warmup_count or len(warmup_process) != warmup_count:
        raise AnalysisError(f"{case_id}: persisted warmup count differs from the plan")
    if len(measured) != sample_count or len(measured_process) != sample_count:
        raise AnalysisError(f"{case_id}: persisted measured count differs from the plan")
    domain = mapping(document.get("domain_result"), f"{case_id}.domain_result")
    domain_sample = _finite(
        domain.get("physical_step_wall_seconds"),
        f"{case_id}.domain_result.physical_step_wall_seconds",
        nonnegative=True,
    )
    if domain_sample != measured[-1]:
        raise AnalysisError(
            f"{case_id}: final domain timing differs from the last measured sample"
        )

    minimum = _finite(
        timing.get("minimum_reliable_seconds"),
        f"{case_id}.timing.minimum_reliable_seconds",
        nonnegative=True,
    )
    if not math.isclose(minimum, requested_minimum, rel_tol=0.0, abs_tol=0.0):
        raise AnalysisError(f"{case_id}: reliability threshold differs from the plan")
    expected_reliability = bool(measured) and min(measured) >= minimum
    if type(timing.get("reliable")) is not bool:
        raise AnalysisError(f"{case_id}.timing.reliable must be a boolean")
    if timing["reliable"] is not expected_reliability:
        raise AnalysisError(f"{case_id}: inconsistent timing reliability flag")

    openmp = mapping(configuration.get("openmp"), f"{case_id}.openmp")
    if openmp.get("dynamic") is not False:
        raise AnalysisError(
            f"{case_id}: {expected_family} scaling requires OMP_DYNAMIC=FALSE"
        )
    proc_bind = openmp.get("proc_bind")
    places = openmp.get("places")
    if not isinstance(proc_bind, str) or not proc_bind:
        raise AnalysisError(f"{case_id}.openmp.proc_bind must be a nonempty string")
    if not isinstance(places, str) or not places:
        raise AnalysisError(f"{case_id}.openmp.places must be a nonempty string")
    environment = mapping(
        timing.get("openmp_environment"),
        f"{case_id}.timing.openmp_environment",
    )
    if set(environment) != OPENMP_ENVIRONMENT_KEYS or any(
        not isinstance(value, str) for value in environment.values()
    ):
        raise AnalysisError(f"{case_id}: invalid persisted OpenMP environment")
    expected_environment = {
        "OMP_NUM_THREADS": str(point.openmp_threads),
        "OMP_DYNAMIC": "FALSE",
        "OMP_PROC_BIND": proc_bind,
        "OMP_PLACES": places,
    }
    if dict(environment) != expected_environment:
        raise AnalysisError(f"{case_id}: persisted OpenMP environment differs from the plan")

    return ScalingPoint(
        field=point,
        warmup_samples=warmups,
        measured_samples=measured,
        warmup_process_wall_seconds=warmup_process,
        measured_process_wall_seconds=measured_process,
        statistics=summarize_samples(measured),
        minimum_reliable_seconds=minimum,
        reliable=expected_reliability,
        openmp_environment=expected_environment,
    )


__all__ = [
    "ScalingPoint",
    "canonical",
    "extract_scaling_point",
    "mapping",
    "positive_integer",
]
