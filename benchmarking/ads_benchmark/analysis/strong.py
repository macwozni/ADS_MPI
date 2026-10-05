"""Strong-scaling analysis gated by complete-field validation."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
import hashlib

from ..framework.model import PlannedCase
from ..validation import FieldComparisonError, compare_fields
from .base import AnalysisError, AnalysisReport
from .scaling import (
    ScalingPoint,
    canonical,
    extract_scaling_point,
    mapping,
)
from .statistics import scaling_metrics
from .validation import (
    PARALLEL_ABSOLUTE_TOLERANCE,
    PARALLEL_RELATIVE_TOLERANCE,
)


def _configuration_without_layout(
    configuration: Mapping[str, object],
) -> dict[str, object]:
    result = dict(configuration)
    result.pop("mpi", None)
    openmp = dict(mapping(result.get("openmp"), "configuration.openmp"))
    openmp.pop("threads", None)
    result["openmp"] = openmp
    return result


def _is_field_reference(point: ScalingPoint) -> bool:
    return (
        point.field.mpi_ranks == 1
        and point.field.mpi_grid == (1, 1, 1)
        and point.field.openmp_threads == 1
    )


def _point_order(point: ScalingPoint) -> tuple[object, ...]:
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
            physics_key = canonical(
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
            encoded = canonical(case.spec.to_dict())
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
            configuration = mapping(
                document.get("configuration"), f"{case_id}.configuration"
            )
            if configuration.get("family") != self.family:
                raise AnalysisError(
                    f"{case_id}: expected {self.family} experiment family"
                )
            encoded = canonical(configuration)
            if encoded in seen_configurations:
                raise AnalysisError(
                    f"duplicate strong configuration at case {case_id}"
                )
            seen_configurations.add(encoded)
            physics = _configuration_without_layout(configuration)
            physics_key = canonical(physics)
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
                extract_scaling_point(document, expected_family=self.family)
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

            by_mpi: dict[tuple[int, tuple[int, int, int]], list[ScalingPoint]] = {}
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
                series_key = canonical(
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
                            mapping(
                                representative.configuration.get("mesh"),
                                "configuration.mesh",
                            )["elements"]
                        ),
                        "test_degree": list(
                            mapping(
                                representative.configuration.get("spaces"),
                                "configuration.spaces",
                            )["test_degree"]
                        ),
                        "trial_degree": list(
                            mapping(
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
