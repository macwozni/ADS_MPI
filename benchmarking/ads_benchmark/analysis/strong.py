"""Strong-scaling analysis gated by complete-field validation."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from decimal import Decimal, InvalidOperation
import hashlib
import json
import math

from ..framework.model import PlannedCase
from ..validation import FieldComparisonError, compare_fields
from .base import AnalysisError, AnalysisReport
from .statistics import SampleStatistics, scaling_metrics, summarize_samples
from .validation import (
    PARALLEL_ABSOLUTE_TOLERANCE,
    PARALLEL_RELATIVE_TOLERANCE,
    FieldResultPoint,
    extract_field_result,
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


def _canonical(value: object) -> str:
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


def _mapping(value: object, field: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise AnalysisError(f"{field} must be an object")
    return value


def _positive_integer(value: object, field: str) -> int:
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
class _StrongPoint:
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


def _extract_strong_point(document: Mapping[str, object]) -> _StrongPoint:
    point = extract_field_result(document, expected_family="strong")
    case_id = point.case_id
    configuration = point.configuration
    measurement = _mapping(
        configuration.get("measurement"), f"{case_id}.measurement"
    )
    warmup_count = _nonnegative_integer(
        measurement.get("warmups"), f"{case_id}.measurement.warmups"
    )
    sample_count = _positive_integer(
        measurement.get("samples"), f"{case_id}.measurement.samples"
    )
    requested_minimum = _decimal_string(
        measurement.get("minimum_sample_seconds"),
        f"{case_id}.measurement.minimum_sample_seconds",
    )

    timing = _mapping(document.get("timing"), f"{case_id}.timing")
    if set(timing) != _REPEATED_TIMING_KEYS:
        missing = sorted(_REPEATED_TIMING_KEYS - set(timing))
        extra = sorted(set(timing) - _REPEATED_TIMING_KEYS)
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
    domain = _mapping(document.get("domain_result"), f"{case_id}.domain_result")
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

    openmp = _mapping(configuration.get("openmp"), f"{case_id}.openmp")
    if openmp.get("dynamic") is not False:
        raise AnalysisError(f"{case_id}: strong scaling requires OMP_DYNAMIC=FALSE")
    proc_bind = openmp.get("proc_bind")
    places = openmp.get("places")
    if not isinstance(proc_bind, str) or not proc_bind:
        raise AnalysisError(f"{case_id}.openmp.proc_bind must be a nonempty string")
    if not isinstance(places, str) or not places:
        raise AnalysisError(f"{case_id}.openmp.places must be a nonempty string")
    environment = _mapping(
        timing.get("openmp_environment"),
        f"{case_id}.timing.openmp_environment",
    )
    if set(environment) != _OPENMP_ENVIRONMENT_KEYS or any(
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

    return _StrongPoint(
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


def _configuration_without_layout(
    configuration: Mapping[str, object],
) -> dict[str, object]:
    result = dict(configuration)
    result.pop("mpi", None)
    openmp = dict(_mapping(result.get("openmp"), "configuration.openmp"))
    openmp.pop("threads", None)
    result["openmp"] = openmp
    return result


def _is_field_reference(point: _StrongPoint) -> bool:
    return (
        point.field.mpi_ranks == 1
        and point.field.mpi_grid == (1, 1, 1)
        and point.field.openmp_threads == 1
    )


def _point_order(point: _StrongPoint) -> tuple[object, ...]:
    return (
        point.resources,
        point.field.mpi_ranks,
        point.field.mpi_grid,
        point.field.openmp_threads,
        point.field.case_id,
    )


@dataclass(frozen=True)
class StrongScalingAnalyzer:
    """Compute scaling only after each field and timing sample is qualified."""

    expected_configurations: tuple[str, ...] | None = None
    name: str = "strong-scaling"
    family: str = "strong"

    def validate_plan(self, planned_cases: Sequence[PlannedCase]) -> None:
        """Require one full-field reference for every selected physics group."""

        if not planned_cases:
            raise AnalysisError("cannot validate an empty strong-scaling plan")
        reference_counts: dict[str, int] = {}
        for case in planned_cases:
            if case.spec.family != self.family:
                raise AnalysisError(
                    f"analyzer {self.name} cannot analyze family {case.spec.family}"
                )
            physics_key = _canonical(
                _configuration_without_layout(case.spec.to_dict())
            )
            reference_counts.setdefault(physics_key, 0)
            if (
                case.spec.mpi.ranks == 1
                and case.spec.mpi.process_grid == (1, 1, 1)
                and case.spec.openmp_threads == 1
            ):
                reference_counts[physics_key] += 1
        invalid = [
            hashlib.sha256(key.encode("utf-8")).hexdigest()[:16]
            for key, count in reference_counts.items()
            if count != 1
        ]
        if invalid:
            raise AnalysisError(
                "strong-scaling plan requires exactly one MPI=1, OMP=1 "
                "full-field reference per physics group; invalid groups="
                + ",".join(sorted(invalid))
            )

    def configure(
        self, planned_cases: Sequence[PlannedCase]
    ) -> "StrongScalingAnalyzer":
        self.validate_plan(planned_cases)
        configurations: set[str] = set()
        for case in planned_cases:
            if case.spec.family != self.family:
                raise AnalysisError(
                    f"analyzer {self.name} cannot analyze family {case.spec.family}"
                )
            encoded = _canonical(case.spec.to_dict())
            if encoded in configurations:
                raise AnalysisError("frozen strong-scaling plan contains a duplicate case")
            configurations.add(encoded)
        return replace(self, expected_configurations=tuple(sorted(configurations)))

    def analyze(
        self,
        results: Iterable[Mapping[str, object]],
        *,
        source_run: str | None = None,
    ) -> AnalysisReport:
        documents = tuple(results)
        if not documents:
            raise AnalysisError("strong-scaling analysis received no results")

        seen_ids: set[str] = set()
        seen_configurations: set[str] = set()
        physics_groups: dict[str, list[Mapping[str, object]]] = {}
        physics_documents: dict[str, Mapping[str, object]] = {}
        for index, document in enumerate(documents):
            if not isinstance(document, Mapping):
                raise AnalysisError(f"result {index} is not an object")
            case_id = document.get("case_id")
            if not isinstance(case_id, str) or not case_id:
                raise AnalysisError(f"result {index}.case_id must be a nonempty string")
            if case_id in seen_ids:
                raise AnalysisError(f"duplicate case_id: {case_id}")
            seen_ids.add(case_id)
            configuration = _mapping(
                document.get("configuration"), f"{case_id}.configuration"
            )
            if configuration.get("family") != self.family:
                raise AnalysisError(
                    f"{case_id}: expected {self.family} experiment family"
                )
            encoded = _canonical(configuration)
            if encoded in seen_configurations:
                raise AnalysisError(
                    f"duplicate strong configuration at case {case_id}"
                )
            seen_configurations.add(encoded)
            physics = _configuration_without_layout(configuration)
            physics_key = _canonical(physics)
            physics_groups.setdefault(physics_key, []).append(document)
            physics_documents[physics_key] = physics

        if self.expected_configurations is not None:
            expected = set(self.expected_configurations)
            if seen_configurations != expected:
                def labels(values: set[str]) -> list[str]:
                    return sorted(
                        hashlib.sha256(value.encode("utf-8")).hexdigest()[:16]
                        for value in values
                    )

                raise AnalysisError(
                    "results do not match frozen strong-scaling plan; "
                    f"missing={labels(expected - seen_configurations)}, "
                    f"unexpected={labels(seen_configurations - expected)}"
                )

        series_documents: list[dict[str, object]] = []
        diagnostics: list[str] = []
        invalid_timing_count = 0
        unreliable_count = 0
        field_failure_count = 0
        failed_series = 0

        for physics_key in sorted(physics_groups):
            # Generated field artifacts are lazy references.  Materialize only
            # one otherwise-identical physics group at a time so the complete
            # cluster matrix does not retain every sampled 3-D field in RAM.
            group = [
                _extract_strong_point(document)
                for document in physics_groups[physics_key]
            ]
            references = [point for point in group if _is_field_reference(point)]
            if len(references) != 1:
                group_id = hashlib.sha256(physics_key.encode("utf-8")).hexdigest()[:16]
                raise AnalysisError(
                    f"{group_id}: strong scaling requires exactly one MPI=1, OMP=1 "
                    "full-field reference"
                )
            field_reference = references[0]

            qualified: dict[str, tuple[bool, Mapping[str, object]]] = {}
            for point in group:
                try:
                    comparison = compare_fields(
                        field_reference.field.numerical,
                        point.field.numerical,
                        absolute_tolerance=PARALLEL_ABSOLUTE_TOLERANCE,
                        relative_tolerance=PARALLEL_RELATIVE_TOLERANCE,
                    )
                except FieldComparisonError as error:
                    raise AnalysisError(
                        f"{point.field.case_id}: full-field comparison failed: {error}"
                    ) from error
                field_valid = comparison.passed
                timing_valid = field_valid and point.reliable and (
                    point.statistics.median > 0.0
                )
                qualified[point.field.case_id] = (timing_valid, comparison.to_dict())
                if not field_valid:
                    field_failure_count += 1
                    diagnostics.append(
                        f"{point.field.case_id}: field differs from MPI=1,OMP=1; "
                        f"{comparison.diagnostic()}"
                    )
                if not point.reliable:
                    unreliable_count += 1
                    diagnostics.append(
                        f"{point.field.case_id}: physical-step sample below "
                        f"{point.minimum_reliable_seconds:.8g} s"
                    )
                if not timing_valid:
                    invalid_timing_count += 1

            eligible = sorted(
                (
                    point
                    for point in group
                    if qualified[point.field.case_id][0]
                ),
                key=_point_order,
            )
            baseline = eligible[0] if eligible else None
            if baseline is None:
                diagnostics.append(
                    "no valid timing baseline for physics group "
                    + hashlib.sha256(physics_key.encode("utf-8")).hexdigest()[:16]
                )

            by_mpi: dict[tuple[int, tuple[int, int, int]], list[_StrongPoint]] = {}
            for point in group:
                key = (point.field.mpi_ranks, point.field.mpi_grid)
                by_mpi.setdefault(key, []).append(point)

            for mpi_key in sorted(by_mpi):
                mpi_points = sorted(
                    by_mpi[mpi_key], key=lambda item: item.field.openmp_threads
                )
                thread_counts = [item.field.openmp_threads for item in mpi_points]
                if len(thread_counts) != len(set(thread_counts)):
                    raise AnalysisError(
                        f"duplicate OpenMP count for MPI layout {mpi_key}"
                    )
                series_key = _canonical(
                    {
                        "physics": physics_documents[physics_key],
                        "mpi_ranks": mpi_key[0],
                        "mpi_grid": list(mpi_key[1]),
                    }
                )
                series_id = hashlib.sha256(series_key.encode("utf-8")).hexdigest()[:16]
                point_documents: list[dict[str, object]] = []
                series_failed = baseline is None
                for point in mpi_points:
                    timing_valid, comparison = qualified[point.field.case_id]
                    metrics = None
                    if timing_valid and baseline is not None:
                        metrics = scaling_metrics(
                            baseline_seconds=baseline.statistics.median,
                            measured_seconds=point.statistics.median,
                            baseline_resources=baseline.resources,
                            resources=point.resources,
                        )
                    else:
                        series_failed = True
                    point_documents.append(
                        {
                            "case_id": point.field.case_id,
                            "mpi_ranks": point.field.mpi_ranks,
                            "mpi_grid": list(point.field.mpi_grid),
                            "openmp_threads": point.field.openmp_threads,
                            "resources": point.resources,
                            "openmp_environment": dict(point.openmp_environment),
                            "warmup_samples": list(point.warmup_samples),
                            "measured_samples": list(point.measured_samples),
                            "warmup_process_wall_seconds": list(
                                point.warmup_process_wall_seconds
                            ),
                            "measured_process_wall_seconds": list(
                                point.measured_process_wall_seconds
                            ),
                            "statistics": point.statistics.to_dict(),
                            "minimum_reliable_seconds": (
                                point.minimum_reliable_seconds
                            ),
                            "measurement_reliable": point.reliable,
                            "field_valid": bool(comparison["passed"]),
                            "timing_valid": timing_valid,
                            "field_comparison": comparison,
                            "speedup": metrics.speedup if metrics else None,
                            "efficiency": metrics.efficiency if metrics else None,
                        }
                    )
                if series_failed:
                    failed_series += 1
                representative = mpi_points[0].field
                series_documents.append(
                    {
                        "group_id": series_id,
                        "physics_group_id": hashlib.sha256(
                            physics_key.encode("utf-8")
                        ).hexdigest()[:16],
                        "series_kind": "strong-scaling",
                        "status": "failed" if series_failed else "passed",
                        "problem": representative.problem,
                        "scheme": representative.scheme,
                        "final_time": float(representative.final_time),
                        "steps": representative.steps,
                        "mesh": list(
                            _mapping(
                                representative.configuration.get("mesh"),
                                "configuration.mesh",
                            )["elements"]
                        ),
                        "test_degree": list(
                            _mapping(
                                representative.configuration.get("spaces"),
                                "configuration.spaces",
                            )["test_degree"]
                        ),
                        "trial_degree": list(
                            _mapping(
                                representative.configuration.get("spaces"),
                                "configuration.spaces",
                            )["trial_degree"]
                        ),
                        "mpi_ranks": mpi_key[0],
                        "mpi_grid": list(mpi_key[1]),
                        "field_reference_case_id": field_reference.field.case_id,
                        "timing_reference_case_id": (
                            baseline.field.case_id if baseline is not None else None
                        ),
                        "timing_reference_resources": (
                            baseline.resources if baseline is not None else None
                        ),
                        "timing_reference_median_seconds": (
                            baseline.statistics.median if baseline is not None else None
                        ),
                        "controlled_configuration": dict(
                            physics_documents[physics_key]
                        ),
                        "points": point_documents,
                    }
                )

        series_documents.sort(
            key=lambda item: (
                str(item["problem"]),
                str(item["scheme"]),
                tuple(item["mesh"]),
                tuple(item["test_degree"]),
                tuple(item["trial_degree"]),
                tuple(item["mpi_grid"]),
            )
        )
        return AnalysisReport(
            analyzer=self.name,
            family=self.family,
            status="passed" if failed_series == 0 else "failed",
            source_run=source_run,
            summary={
                "series_count": len(series_documents),
                "passed_series": len(series_documents) - failed_series,
                "failed_series": failed_series,
                "case_count": len(documents),
                "invalid_timing_count": invalid_timing_count,
                "unreliable_measurement_count": unreliable_count,
                "field_failure_count": field_failure_count,
                "parallel_absolute_tolerance": PARALLEL_ABSOLUTE_TOLERANCE,
                "parallel_relative_tolerance": PARALLEL_RELATIVE_TOLERANCE,
                "primary_metric": "physical_step_wall_seconds",
                "spread_metric": "median_absolute_deviation",
                "timing_reference_policy": (
                    "minimum valid ranks*threads; ties by ranks, MPI grid, "
                    "threads, case_id"
                ),
            },
            series=tuple(series_documents),
            diagnostics=tuple(diagnostics),
        )


__all__ = ["StrongScalingAnalyzer"]
