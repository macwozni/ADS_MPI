"""Weak-scaling analysis gated by same-mesh complete-field validation."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from fractions import Fraction
import hashlib
import math

from ..framework.model import PlannedCase
from ..validation import FieldComparisonError, compare_fields
from .base import AnalysisError, AnalysisReport
from .scaling import (
    ScalingPoint,
    canonical,
    extract_scaling_point,
    mapping,
    positive_integer,
)
from .statistics import weak_scaling_efficiency
from .validation import (
    PARALLEL_ABSOLUTE_TOLERANCE,
    PARALLEL_RELATIVE_TOLERANCE,
)


Vector3 = tuple[int, int, int]
_WEAK_SCALING_KEYS = {"local_elements", "workload_basis", "role"}
_WEAK_ROLES = {"measurement", "field-reference"}


def _vector3(value: object, field: str) -> Vector3:
    if not isinstance(value, list) or len(value) != 3:
        raise AnalysisError(f"{field} must contain exactly three integers")
    result = tuple(
        positive_integer(component, f"{field}[{index}]")
        for index, component in enumerate(value)
    )
    return result  # type: ignore[return-value]


@dataclass(frozen=True)
class _WeakMetadata:
    local_elements: Vector3
    workload_basis: str
    role: str


def _weak_metadata(
    configuration: Mapping[str, object], case_id: str
) -> _WeakMetadata:
    weak = mapping(
        configuration.get("weak_scaling"), f"{case_id}.weak_scaling"
    )
    if set(weak) != _WEAK_SCALING_KEYS:
        missing = sorted(_WEAK_SCALING_KEYS - set(weak))
        extra = sorted(set(weak) - _WEAK_SCALING_KEYS)
        raise AnalysisError(
            f"{case_id}.weak_scaling has invalid keys: "
            f"missing={missing}, extra={extra}"
        )
    local_elements = _vector3(
        weak.get("local_elements"), f"{case_id}.weak_scaling.local_elements"
    )
    workload_basis = weak.get("workload_basis")
    if workload_basis != "per-rank":
        raise AnalysisError(
            f"{case_id}: unsupported weak-scaling workload basis "
            f"{workload_basis!r}; expected 'per-rank'"
        )
    role = weak.get("role")
    if role not in _WEAK_ROLES:
        raise AnalysisError(
            f"{case_id}.weak_scaling.role must be 'measurement' or "
            "'field-reference'"
        )
    return _WeakMetadata(local_elements, workload_basis, str(role))


def _validate_weak_configuration(
    configuration: Mapping[str, object], case_id: str
) -> _WeakMetadata:
    metadata = _weak_metadata(configuration, case_id)
    mesh = _vector3(
        mapping(configuration.get("mesh"), f"{case_id}.mesh").get("elements"),
        f"{case_id}.mesh.elements",
    )
    mpi = mapping(configuration.get("mpi"), f"{case_id}.mpi")
    ranks = positive_integer(mpi.get("ranks"), f"{case_id}.mpi.ranks")
    grid = _vector3(
        mpi.get("process_grid"), f"{case_id}.mpi.process_grid"
    )
    if math.prod(grid) != ranks:
        raise AnalysisError(f"{case_id}: MPI ranks differ from process-grid product")
    expected_mesh = tuple(
        local * processes
        for local, processes in zip(metadata.local_elements, grid, strict=True)
    )
    if mesh != expected_mesh:
        raise AnalysisError(
            f"{case_id}: global mesh {mesh} does not equal per-rank local "
            f"elements {metadata.local_elements} times process grid {grid}"
        )
    openmp = mapping(configuration.get("openmp"), f"{case_id}.openmp")
    threads = positive_integer(openmp.get("threads"), f"{case_id}.openmp.threads")
    if metadata.role == "field-reference" and (
        ranks != 1 or grid != (1, 1, 1) or threads != 1
    ):
        raise AnalysisError(
            f"{case_id}: weak field-reference must use MPI=1, grid=1x1x1, OMP=1"
        )
    return metadata


def _configuration_without_resource_layout(
    configuration: Mapping[str, object],
) -> dict[str, object]:
    """Identity of one per-rank weak series; OMP count stays controlled."""

    result = dict(configuration)
    result.pop("mpi", None)
    result.pop("mesh", None)
    weak = dict(mapping(result.get("weak_scaling"), "configuration.weak_scaling"))
    weak.pop("role", None)
    result["weak_scaling"] = weak
    return result


def _configuration_for_field_reference(
    configuration: Mapping[str, object],
) -> dict[str, object]:
    """Identity of a discrete field, independent of its parallel layout."""

    result = dict(configuration)
    result.pop("mpi", None)
    result.pop("weak_scaling", None)
    openmp = dict(mapping(result.get("openmp"), "configuration.openmp"))
    openmp.pop("threads", None)
    result["openmp"] = openmp
    return result


def _is_serial_field_reference(point: ScalingPoint) -> bool:
    return (
        point.field.mpi_ranks == 1
        and point.field.mpi_grid == (1, 1, 1)
        and point.field.openmp_threads == 1
    )


def _measurement_order(point: Mapping[str, object]) -> tuple[object, ...]:
    return (
        point["resources"],
        point["mpi_ranks"],
        tuple(point["mpi_grid"]),  # type: ignore[arg-type]
        point["case_id"],
    )


@dataclass(frozen=True)
class WeakScalingAnalyzer:
    """Compute per-rank weak efficiency after complete-field qualification."""

    expected_configurations: tuple[str, ...] | None = None
    name: str = "weak-scaling"
    family: str = "weak"

    def validate_plan(self, planned_cases: Sequence[PlannedCase]) -> None:
        if not planned_cases:
            raise AnalysisError("cannot validate an empty weak-scaling plan")

        field_reference_counts: dict[str, int] = {}
        timing_groups: dict[str, list[Mapping[str, object]]] = {}
        configurations: set[str] = set()
        for case in planned_cases:
            if case.spec.family != self.family:
                raise AnalysisError(
                    f"analyzer {self.name} cannot analyze family {case.spec.family}"
                )
            configuration = case.spec.to_dict()
            metadata = _validate_weak_configuration(configuration, case.case_id)
            encoded = canonical(configuration)
            if encoded in configurations:
                raise AnalysisError("frozen weak-scaling plan contains a duplicate case")
            configurations.add(encoded)

            field_key = canonical(_configuration_for_field_reference(configuration))
            field_reference_counts.setdefault(field_key, 0)
            if (
                case.spec.mpi.ranks == 1
                and case.spec.mpi.process_grid == (1, 1, 1)
                and case.spec.openmp_threads == 1
            ):
                field_reference_counts[field_key] += 1
            if metadata.role == "measurement":
                timing_key = canonical(
                    _configuration_without_resource_layout(configuration)
                )
                timing_groups.setdefault(timing_key, []).append(configuration)

        invalid_fields = [
            hashlib.sha256(key.encode("utf-8")).hexdigest()[:16]
            for key, count in field_reference_counts.items()
            if count != 1
        ]
        if invalid_fields:
            raise AnalysisError(
                "weak-scaling plan requires exactly one same-mesh MPI=1, OMP=1 "
                "full-field reference per discrete configuration; invalid groups="
                + ",".join(sorted(invalid_fields))
            )
        if not timing_groups:
            raise AnalysisError("weak-scaling plan contains no measurement cases")
        missing_baselines: list[str] = []
        for key, configurations_in_group in timing_groups.items():
            has_baseline = any(
                mapping(configuration["mpi"], "configuration.mpi").get("ranks")
                == 1
                for configuration in configurations_in_group
            )
            if not has_baseline:
                missing_baselines.append(
                    hashlib.sha256(key.encode("utf-8")).hexdigest()[:16]
                )
        if missing_baselines:
            raise AnalysisError(
                "weak-scaling plan requires an MPI=1 measurement baseline in "
                "every local-work/OMP series; missing groups="
                + ",".join(sorted(missing_baselines))
            )

    def configure(
        self, planned_cases: Sequence[PlannedCase]
    ) -> "WeakScalingAnalyzer":
        self.validate_plan(planned_cases)
        configurations = tuple(
            sorted(canonical(case.spec.to_dict()) for case in planned_cases)
        )
        return replace(self, expected_configurations=configurations)

    def analyze(
        self,
        results: Iterable[Mapping[str, object]],
        *,
        source_run: str | None = None,
    ) -> AnalysisReport:
        documents = tuple(results)
        if not documents:
            raise AnalysisError("weak-scaling analysis received no results")

        seen_ids: set[str] = set()
        seen_configurations: set[str] = set()
        field_groups: dict[str, list[Mapping[str, object]]] = {}
        metadata_by_case: dict[str, _WeakMetadata] = {}
        timing_key_by_case: dict[str, str] = {}
        timing_configurations: dict[str, Mapping[str, object]] = {}
        measurement_count = 0
        reference_count = 0

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
            metadata = _validate_weak_configuration(configuration, case_id)
            metadata_by_case[case_id] = metadata
            encoded = canonical(configuration)
            if encoded in seen_configurations:
                raise AnalysisError(f"duplicate weak configuration at case {case_id}")
            seen_configurations.add(encoded)

            field_configuration = _configuration_for_field_reference(configuration)
            field_key = canonical(field_configuration)
            field_groups.setdefault(field_key, []).append(document)
            if metadata.role == "measurement":
                measurement_count += 1
                timing_configuration = _configuration_without_resource_layout(
                    configuration
                )
                timing_key = canonical(timing_configuration)
                timing_key_by_case[case_id] = timing_key
                timing_configurations[timing_key] = timing_configuration
            else:
                reference_count += 1

        if self.expected_configurations is not None:
            expected = set(self.expected_configurations)
            if seen_configurations != expected:
                def labels(values: set[str]) -> list[str]:
                    return sorted(
                        hashlib.sha256(value.encode("utf-8")).hexdigest()[:16]
                        for value in values
                    )

                raise AnalysisError(
                    "results do not match frozen weak-scaling plan; "
                    f"missing={labels(expected - seen_configurations)}, "
                    f"unexpected={labels(seen_configurations - expected)}"
                )
        if measurement_count == 0:
            raise AnalysisError("weak-scaling results contain no measurement cases")

        timing_points: dict[str, list[dict[str, object]]] = {}
        diagnostics: list[str] = []
        invalid_timing_count = 0
        unreliable_count = 0
        field_failure_count = 0

        for field_key in sorted(field_groups):
            # Materialize one global-mesh field group at a time.  Persist only
            # scalar comparison/timing documents, never all sampled 3-D fields.
            points = [
                extract_scaling_point(document, expected_family=self.family)
                for document in field_groups[field_key]
            ]
            references = [point for point in points if _is_serial_field_reference(point)]
            if len(references) != 1:
                group_id = hashlib.sha256(field_key.encode("utf-8")).hexdigest()[:16]
                raise AnalysisError(
                    f"{group_id}: weak scaling requires exactly one same-mesh "
                    "MPI=1, OMP=1 full-field reference"
                )
            field_reference = references[0]

            for point in points:
                metadata = metadata_by_case[point.field.case_id]
                if metadata.role != "measurement":
                    continue
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
                if not field_valid:
                    field_failure_count += 1
                    diagnostics.append(
                        f"{point.field.case_id}: field differs from same-mesh "
                        f"MPI=1,OMP=1 reference; {comparison.diagnostic()}"
                    )
                if not point.reliable:
                    unreliable_count += 1
                    diagnostics.append(
                        f"{point.field.case_id}: physical-step sample below "
                        f"{point.minimum_reliable_seconds:.8g} s"
                    )
                if not timing_valid:
                    invalid_timing_count += 1

                mesh = _vector3(
                    mapping(
                        point.field.configuration.get("mesh"),
                        f"{point.field.case_id}.mesh",
                    ).get("elements"),
                    f"{point.field.case_id}.mesh.elements",
                )
                point_document: dict[str, object] = {
                    "case_id": point.field.case_id,
                    "mpi_ranks": point.field.mpi_ranks,
                    "mpi_grid": list(point.field.mpi_grid),
                    "openmp_threads": point.field.openmp_threads,
                    "resources": point.resources,
                    "global_mesh": list(mesh),
                    "local_elements": list(metadata.local_elements),
                    "workload_basis": metadata.workload_basis,
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
                    "minimum_reliable_seconds": point.minimum_reliable_seconds,
                    "measurement_reliable": point.reliable,
                    "field_valid": bool(comparison.passed),
                    "timing_valid": timing_valid,
                    "field_comparison": comparison.to_dict(),
                    "field_reference_case_id": field_reference.field.case_id,
                    "weak_scaling_efficiency": None,
                }
                timing_key = timing_key_by_case[point.field.case_id]
                timing_points.setdefault(timing_key, []).append(point_document)

        series_documents: list[dict[str, object]] = []
        failed_series = 0
        for timing_key in sorted(timing_points):
            points = sorted(timing_points[timing_key], key=_measurement_order)
            resources = [int(point["resources"]) for point in points]
            if len(resources) != len(set(resources)):
                raise AnalysisError(
                    "duplicate resource count in weak-scaling local-work/OMP series "
                    + hashlib.sha256(timing_key.encode("utf-8")).hexdigest()[:16]
                )
            baseline_points = [
                point for point in points if point["mpi_ranks"] == 1
            ]
            if len(baseline_points) != 1:
                raise AnalysisError(
                    "weak-scaling local-work/OMP series requires exactly one "
                    "MPI=1 timing baseline "
                    + hashlib.sha256(timing_key.encode("utf-8")).hexdigest()[:16]
                )
            baseline = (
                baseline_points[0]
                if baseline_points[0]["timing_valid"]
                else None
            )
            series_failed = baseline is None
            if baseline is None:
                diagnostics.append(
                    "no valid MPI=1 timing baseline for weak-scaling series "
                    + hashlib.sha256(timing_key.encode("utf-8")).hexdigest()[:16]
                )
            else:
                baseline_seconds = float(
                    mapping(baseline["statistics"], "statistics")["median_seconds"]
                )
                for point in points:
                    if point["timing_valid"]:
                        measured_seconds = float(
                            mapping(point["statistics"], "statistics")[
                                "median_seconds"
                            ]
                        )
                        point["weak_scaling_efficiency"] = weak_scaling_efficiency(
                            baseline_seconds=baseline_seconds,
                            measured_seconds=measured_seconds,
                        )
                    else:
                        series_failed = True
            if series_failed:
                failed_series += 1

            configuration = timing_configurations[timing_key]
            weak = mapping(configuration["weak_scaling"], "weak_scaling")
            spaces = mapping(configuration["spaces"], "spaces")
            time = mapping(configuration["time"], "time")
            openmp = mapping(configuration["openmp"], "openmp")
            series_id = hashlib.sha256(timing_key.encode("utf-8")).hexdigest()[:16]
            series_documents.append(
                {
                    "group_id": series_id,
                    "series_kind": "weak-scaling",
                    "status": "failed" if series_failed else "passed",
                    "problem": configuration["problem"],
                    "scheme": configuration["scheme"],
                    "final_time": float(Fraction(str(time["final_time"]))),
                    "steps": time["steps"],
                    "test_degree": list(spaces["test_degree"]),
                    "trial_degree": list(spaces["trial_degree"]),
                    "local_elements": list(weak["local_elements"]),
                    "workload_basis": weak["workload_basis"],
                    "openmp_threads": openmp["threads"],
                    "timing_reference_case_id": (
                        baseline["case_id"] if baseline is not None else None
                    ),
                    "timing_reference_resources": (
                        baseline["resources"] if baseline is not None else None
                    ),
                    "timing_reference_median_seconds": (
                        mapping(baseline["statistics"], "statistics")[
                            "median_seconds"
                        ]
                        if baseline is not None
                        else None
                    ),
                    "controlled_configuration": dict(configuration),
                    "points": points,
                }
            )

        series_documents.sort(
            key=lambda item: (
                str(item["problem"]),
                str(item["scheme"]),
                tuple(item["local_elements"]),
                tuple(item["test_degree"]),
                tuple(item["trial_degree"]),
                int(item["openmp_threads"]),
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
                "measurement_case_count": measurement_count,
                "field_reference_case_count": reference_count,
                "invalid_timing_count": invalid_timing_count,
                "unreliable_measurement_count": unreliable_count,
                "field_failure_count": field_failure_count,
                "parallel_absolute_tolerance": PARALLEL_ABSOLUTE_TOLERANCE,
                "parallel_relative_tolerance": PARALLEL_RELATIVE_TOLERANCE,
                "primary_metric": "physical_step_wall_seconds",
                "spread_metric": "median_absolute_deviation",
                "efficiency_metric": "weak_scaling_efficiency",
                "workload_policy": "constant elements per MPI rank",
                "timing_reference_policy": (
                    "mandatory valid MPI=1 measurement within identical "
                    "physics, per-rank local elements, OpenMP count and "
                    "binding policy"
                ),
                "field_reference_policy": (
                    "same global mesh and physics at MPI=1, OMP=1"
                ),
            },
            series=tuple(series_documents),
            diagnostics=tuple(diagnostics),
        )


__all__ = ["WeakScalingAnalyzer"]
