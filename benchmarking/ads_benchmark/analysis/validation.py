"""Full-field numerical validation across MPI/OpenMP layouts and schemes.

The analyzer treats a serial ``MPI=1, OMP=1`` result as the parallelism
reference for one otherwise identical configuration.  It compares every
sample before making the corresponding timing eligible.  Serial references
also provide the independent analytical and cross-scheme audits.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from fractions import Fraction
import hashlib
import itertools
import json
import math

from ..framework.model import PlannedCase
from ..validation import (
    FieldComparison,
    FieldComparisonError,
    RegularGridField,
    compare_fields,
    parse_regular_grid_csv,
    scalar_component,
)
from .base import AnalysisError, AnalysisReport


VALIDATION_EXACT_CASE = "spatial-cosine"

# Parallel layouts execute the same discrete problem.  These limits cover
# binary64 reduction/reconstruction noise without hiding a changed solution.
PARALLEL_ABSOLUTE_TOLERANCE = 1.0e-11
PARALLEL_RELATIVE_TOLERANCE = 1.0e-10

# The frozen validation workload uses the unit-amplitude spatial-cosine case
# on a deliberately small mesh.  Its declared final analytic accuracy is an
# absolute five percent of that unit amplitude.  Keeping the relative term at
# zero makes the policy exactly 0.05 rather than an accidental abs+rel 0.10
# allowance near extrema.
ANALYTIC_ABSOLUTE_TOLERANCE = 5.0e-2
ANALYTIC_RELATIVE_TOLERANCE = 0.0

# DG/PR/BE must agree pointwise to one percent of unit amplitude at the finest
# planned dt.  Neither pairwise norm may grow by more than 0.1% from the
# coarsest level, and at least one must contract by a measurable 0.1% (unless
# the final fields are already bitwise identical).
SCHEME_ABSOLUTE_TOLERANCE = 1.0e-2
SCHEME_RELATIVE_TOLERANCE = 0.0
SCHEME_TREND_RELAXATION = 1.0e-3
SCHEME_TREND_MINIMUM_REDUCTION = 1.0e-3


Vector3 = tuple[int, int, int]


@dataclass(frozen=True)
class _ValidationPoint:
    case_id: str
    configuration: Mapping[str, object]
    problem: str
    scheme: str
    final_time: Fraction
    time_step: Fraction
    steps: int
    mpi_ranks: int
    mpi_grid: Vector3
    openmp_threads: int
    wall_seconds: float
    reported_l2_error: float
    reported_linf_error: float
    reported_solution_l2_norm: float
    reported_field_checksum: float
    field: RegularGridField
    numerical: RegularGridField
    exact: RegularGridField


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise AnalysisError(f"{name} must be an object")
    return value


def _string(value: object, name: str) -> str:
    if not isinstance(value, str) or not value:
        raise AnalysisError(f"{name} must be a nonempty string")
    return value


def _integer(value: object, name: str, *, positive: bool = False) -> int:
    if type(value) is not int or (positive and value <= 0):
        qualifier = "positive " if positive else ""
        raise AnalysisError(f"{name} must be a {qualifier}integer")
    return value


def _finite(value: object, name: str, *, nonnegative: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise AnalysisError(f"{name} must be a finite number")
    result = float(value)
    if not math.isfinite(result) or (nonnegative and result < 0.0):
        qualifier = " nonnegative" if nonnegative else ""
        raise AnalysisError(f"{name} must be a finite{qualifier} number")
    return result


def _fraction(value: object, name: str) -> Fraction:
    if not isinstance(value, str):
        raise AnalysisError(f"{name} must be an exact-number string")
    try:
        result = Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise AnalysisError(f"{name} is not an exact number") from error
    if result <= 0:
        raise AnalysisError(f"{name} must be positive")
    return result


def _vector3(value: object, name: str) -> Vector3:
    if not isinstance(value, list) or len(value) != 3:
        raise AnalysisError(f"{name} must contain exactly three integers")
    result = tuple(
        _integer(component, f"{name}[{index}]", positive=True)
        for index, component in enumerate(value)
    )
    return result  # type: ignore[return-value]


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


def _matches_exact_float(actual: float, expected: Fraction) -> bool:
    expected_float = float(expected)
    tolerance = max(64.0 * math.ulp(expected_float), 1.0e-14 * expected_float)
    return abs(actual - expected_float) <= tolerance


def _artifact(document: Mapping[str, object], case_id: str) -> object:
    artifacts = _mapping(
        document.get("analysis_artifacts"), f"{case_id}.analysis_artifacts"
    )
    artifact = artifacts.get("field_samples_csv")
    if isinstance(artifact, str) or callable(getattr(artifact, "read_text", None)):
        return artifact
    raise AnalysisError(f"{case_id}: field_samples.csv is unavailable")


def _validate_sample_contract(
    point_field: RegularGridField,
    *,
    case_id: str,
    final_time: Fraction,
    domain: Mapping[str, object],
) -> None:
    if point_field.component_names != ("numerical", "exact", "error"):
        raise AnalysisError(
            f"{case_id}: field_samples.csv header must be "
            "x,y,z,numerical,exact,error"
        )
    for axis_name, axis in zip(("x", "y", "z"), point_field.axes, strict=True):
        if not math.isclose(axis[0], 0.0, rel_tol=0.0, abs_tol=2.0e-14) or not math.isclose(
            axis[-1], 1.0, rel_tol=0.0, abs_tol=2.0e-14
        ):
            raise AnalysisError(
                f"{case_id}: sampled {axis_name} axis must cover the physical [0,1] domain"
            )

    numerical = point_field.component_values("numerical")
    exact = point_field.component_values("exact")
    errors = point_field.component_values("error")
    physical_time = float(final_time)
    observed_linf = 0.0
    sample_sum = 0.0
    sample_compensation = 0.0
    weighted_sum = 0.0
    weighted_compensation = 0.0
    for index, (actual, expected, error_value) in enumerate(
        zip(numerical, exact, errors, strict=True)
    ):
        x, y, z = point_field.coordinates(index)
        analytical = (
            math.exp(-physical_time)
            * math.cos(math.pi * x)
            * math.cos(math.pi * y)
            * math.cos(math.pi * z)
        )
        if not math.isclose(
            expected, analytical, rel_tol=2.0e-13, abs_tol=2.0e-14
        ):
            raise AnalysisError(
                f"{case_id}: wrong exact field value at coordinate {(x, y, z)}"
            )
        if not math.isclose(
            error_value,
            actual - expected,
            rel_tol=2.0e-11,
            abs_tol=5.0e-14,
        ):
            raise AnalysisError(
                f"{case_id}: inconsistent error field at coordinate {(x, y, z)}"
            )
        observed_linf = max(observed_linf, abs(error_value))
        corrected = actual - sample_compensation
        updated = sample_sum + corrected
        sample_compensation = (updated - sample_sum) - corrected
        sample_sum = updated
        corrected = (index + 1) * actual - weighted_compensation
        updated = weighted_sum + corrected
        weighted_compensation = (updated - weighted_sum) - corrected
        weighted_sum = updated

    reported_linf = _finite(
        domain.get("linf_error"), f"{case_id}.domain.linf_error", nonnegative=True
    )
    if not math.isclose(
        observed_linf, reported_linf, rel_tol=2.0e-11, abs_tol=5.0e-14
    ):
        raise AnalysisError(
            f"{case_id}: sampled Linf error differs from the domain result"
        )
    observed_checksum = sample_sum + weighted_sum / (point_field.point_count + 1)
    reported_checksum = _finite(
        domain.get("field_checksum"), f"{case_id}.domain.field_checksum"
    )
    if not math.isclose(
        observed_checksum,
        reported_checksum,
        rel_tol=2.0e-11,
        abs_tol=2.0e-12,
    ):
        raise AnalysisError(
            f"{case_id}: field sample checksum differs from the domain result"
        )


def _extract_point(document: Mapping[str, object]) -> _ValidationPoint:
    case_id = _string(document.get("case_id"), "case_id")
    if document.get("kind") != "ads-benchmark-case-result":
        raise AnalysisError(f"{case_id}: unsupported result kind")
    if document.get("status") != "passed":
        raise AnalysisError(f"{case_id}: result status is not passed")
    configuration = _mapping(
        document.get("configuration"), f"{case_id}.configuration"
    )
    if configuration.get("family") != "validation":
        raise AnalysisError(f"{case_id}: expected validation experiment family")
    problem = _string(configuration.get("problem"), f"{case_id}.problem")
    scheme = _string(configuration.get("scheme"), f"{case_id}.scheme").lower()
    exact_case = _string(
        configuration.get("exact_case"), f"{case_id}.exact_case"
    )
    if exact_case != VALIDATION_EXACT_CASE:
        raise AnalysisError(
            f"{case_id}: validation requires exact case {VALIDATION_EXACT_CASE!r}"
        )

    time = _mapping(configuration.get("time"), f"{case_id}.time")
    final_time = _fraction(time.get("final_time"), f"{case_id}.final_time")
    time_step = _fraction(time.get("time_step"), f"{case_id}.time_step")
    steps = _integer(time.get("steps"), f"{case_id}.steps", positive=True)
    if time_step * steps != final_time:
        raise AnalysisError(f"{case_id}: inconsistent T, dt, and step count")

    mpi = _mapping(configuration.get("mpi"), f"{case_id}.mpi")
    mpi_ranks = _integer(mpi.get("ranks"), f"{case_id}.mpi.ranks", positive=True)
    mpi_grid = _vector3(mpi.get("process_grid"), f"{case_id}.mpi.process_grid")
    if math.prod(mpi_grid) != mpi_ranks:
        raise AnalysisError(f"{case_id}: MPI ranks differ from the process-grid product")
    openmp = _mapping(configuration.get("openmp"), f"{case_id}.openmp")
    openmp_threads = _integer(
        openmp.get("threads"), f"{case_id}.openmp.threads", positive=True
    )
    sampling = _mapping(configuration.get("sampling"), f"{case_id}.sampling")
    points = _integer(
        sampling.get("points_per_axis"),
        f"{case_id}.sampling.points_per_axis",
        positive=True,
    )
    if points < 2 or sampling.get("write_samples") is not True:
        raise AnalysisError(
            f"{case_id}: full-field validation requires write_samples=true and at least two points"
        )

    domain = _mapping(document.get("domain_result"), f"{case_id}.domain_result")
    if domain.get("problem") != problem or domain.get("scheme") != scheme:
        raise AnalysisError(f"{case_id}: domain problem/scheme differs from the plan")
    if domain.get("exact_case") != exact_case:
        raise AnalysisError(f"{case_id}: domain exact case differs from the plan")
    if _integer(domain.get("solver_status"), f"{case_id}.solver_status") != 0:
        raise AnalysisError(f"{case_id}: solver_status is nonzero")
    if _integer(domain.get("steps"), f"{case_id}.domain.steps") != steps:
        raise AnalysisError(f"{case_id}: domain step count differs from the plan")
    if _integer(
        domain.get("sample_points_per_axis"),
        f"{case_id}.domain.sample_points_per_axis",
    ) != points:
        raise AnalysisError(f"{case_id}: domain sample grid differs from the plan")
    if domain.get("field_samples_written") is not True:
        raise AnalysisError(f"{case_id}: domain did not write field samples")
    if not _matches_exact_float(
        _finite(domain.get("actual_final_time"), f"{case_id}.actual_final_time"),
        final_time,
    ) or not _matches_exact_float(
        _finite(
            domain.get("requested_final_time"),
            f"{case_id}.requested_final_time",
        ),
        final_time,
    ):
        raise AnalysisError(f"{case_id}: domain final time differs from the plan")
    if not _matches_exact_float(
        _finite(domain.get("time_step"), f"{case_id}.domain.time_step"),
        time_step,
    ):
        raise AnalysisError(f"{case_id}: domain time step differs from the plan")

    try:
        field = parse_regular_grid_csv(
            _artifact(document, case_id), shape=points, label=case_id
        )
        _validate_sample_contract(
            field, case_id=case_id, final_time=final_time, domain=domain
        )
        numerical = scalar_component(
            field, "numerical", name="value", label=f"{case_id}:numerical"
        )
        exact = scalar_component(
            field, "exact", name="value", label=f"{case_id}:exact"
        )
    except FieldComparisonError as error:
        raise AnalysisError(f"{case_id}: {error}") from error

    timing = _mapping(document.get("timing"), f"{case_id}.timing")
    return _ValidationPoint(
        case_id=case_id,
        configuration=configuration,
        problem=problem,
        scheme=scheme,
        final_time=final_time,
        time_step=time_step,
        steps=steps,
        mpi_ranks=mpi_ranks,
        mpi_grid=mpi_grid,
        openmp_threads=openmp_threads,
        wall_seconds=_finite(
            timing.get("wall_seconds"), f"{case_id}.timing.wall_seconds", nonnegative=True
        ),
        reported_l2_error=_finite(
            domain.get("l2_error"), f"{case_id}.l2_error", nonnegative=True
        ),
        reported_linf_error=_finite(
            domain.get("linf_error"), f"{case_id}.linf_error", nonnegative=True
        ),
        reported_solution_l2_norm=_finite(
            domain.get("solution_l2_norm"),
            f"{case_id}.solution_l2_norm",
            nonnegative=True,
        ),
        reported_field_checksum=_finite(
            domain.get("field_checksum"), f"{case_id}.field_checksum"
        ),
        field=field,
        numerical=numerical,
        exact=exact,
    )


def _configuration_without_layout(configuration: Mapping[str, object]) -> dict[str, object]:
    result = dict(configuration)
    result.pop("mpi", None)
    result.pop("openmp", None)
    return result


def _physics_configuration(configuration: Mapping[str, object]) -> dict[str, object]:
    result = _configuration_without_layout(configuration)
    result.pop("scheme", None)
    time = dict(_mapping(result.get("time"), "configuration.time"))
    time.pop("steps", None)
    time.pop("time_step", None)
    result["time"] = time
    return result


def _layout_signature(point: _ValidationPoint) -> tuple[int, Vector3, int]:
    return point.mpi_ranks, point.mpi_grid, point.openmp_threads


def _is_reference(point: _ValidationPoint) -> bool:
    return (
        point.mpi_ranks == 1
        and point.mpi_grid == (1, 1, 1)
        and point.openmp_threads == 1
    )


def _comparison_document(comparison: FieldComparison) -> dict[str, object]:
    return comparison.to_dict()


def _case_document(point: _ValidationPoint, *, timing_valid: bool) -> dict[str, object]:
    return {
        "case_id": point.case_id,
        "mpi_ranks": point.mpi_ranks,
        "mpi_grid": list(point.mpi_grid),
        "openmp_threads": point.openmp_threads,
        "wall_seconds": point.wall_seconds,
        "timing_valid": timing_valid,
        "reported_l2_error": point.reported_l2_error,
        "reported_linf_error": point.reported_linf_error,
        "reported_solution_l2_norm": point.reported_solution_l2_norm,
        "reported_field_checksum": point.reported_field_checksum,
    }


@dataclass(frozen=True)
class FieldValidationAnalyzer:
    """Validate complete fields before accepting any parallel timing."""

    expected_configurations: tuple[str, ...] | None = None
    name: str = "field-validation"
    family: str = "validation"

    def configure(
        self, planned_cases: Sequence[PlannedCase]
    ) -> "FieldValidationAnalyzer":
        if not planned_cases:
            raise AnalysisError("cannot configure analyzer from an empty plan")
        configurations: set[str] = set()
        for case in planned_cases:
            if case.spec.family != self.family:
                raise AnalysisError(
                    f"analyzer {self.name} cannot analyze family {case.spec.family}"
                )
            encoded = _canonical(case.spec.to_dict())
            if encoded in configurations:
                raise AnalysisError("frozen validation plan contains a duplicate case")
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
            raise AnalysisError("field validation received no results")

        seen_ids: set[str] = set()
        seen_configurations: set[str] = set()
        physics_groups: dict[str, list[_ValidationPoint]] = {}
        physics_documents: dict[str, Mapping[str, object]] = {}
        for index, document in enumerate(documents):
            if not isinstance(document, Mapping):
                raise AnalysisError(f"result {index} is not an object")
            point = _extract_point(document)
            if point.case_id in seen_ids:
                raise AnalysisError(f"duplicate case_id: {point.case_id}")
            seen_ids.add(point.case_id)
            encoded = _canonical(point.configuration)
            if encoded in seen_configurations:
                raise AnalysisError(
                    f"duplicate validation configuration at case {point.case_id}"
                )
            seen_configurations.add(encoded)
            physics = _physics_configuration(point.configuration)
            key = _canonical(physics)
            physics_groups.setdefault(key, []).append(point)
            physics_documents[key] = physics

        if self.expected_configurations is not None:
            expected = set(self.expected_configurations)
            if seen_configurations != expected:
                def labels(values: set[str]) -> list[str]:
                    return sorted(
                        hashlib.sha256(value.encode("utf-8")).hexdigest()[:16]
                        for value in values
                    )

                raise AnalysisError(
                    "results do not match frozen validation plan; "
                    f"missing={labels(expected - seen_configurations)}, "
                    f"unexpected={labels(seen_configurations - expected)}"
                )

        series_documents: list[dict[str, object]] = []
        diagnostics: list[str] = []
        failed_series = 0
        parallel_comparisons = 0
        invalid_timings = 0
        analytic_final_failures = 0
        scheme_pair_failures = 0

        for physics_key in sorted(physics_groups):
            points = physics_groups[physics_key]
            group_id = hashlib.sha256(physics_key.encode("utf-8")).hexdigest()[:16]
            by_step_scheme: dict[int, dict[str, list[_ValidationPoint]]] = {}
            for point in points:
                by_step_scheme.setdefault(point.steps, {}).setdefault(
                    point.scheme, []
                ).append(point)
            steps = sorted(by_step_scheme)
            scheme_sets = {tuple(sorted(groups)) for groups in by_step_scheme.values()}
            if len(scheme_sets) != 1:
                raise AnalysisError(
                    f"{group_id}: inconsistent scheme sets across time levels"
                )
            schemes = next(iter(scheme_sets))
            if len(schemes) > 1 and len(steps) < 2:
                raise AnalysisError(
                    f"{group_id}: cross-scheme validation requires at least two dt levels"
                )

            expected_layouts: set[tuple[int, Vector3, int]] | None = None
            level_documents: list[dict[str, object]] = []
            pair_history: dict[tuple[str, str], list[FieldComparison]] = {
                pair: [] for pair in itertools.combinations(schemes, 2)
            }
            series_failed = False

            for level_index, step_count in enumerate(steps):
                scheme_points = by_step_scheme[step_count]
                references: dict[str, _ValidationPoint] = {}
                scheme_documents: list[dict[str, object]] = []
                time_steps = {
                    point.time_step
                    for values in scheme_points.values()
                    for point in values
                }
                if len(time_steps) != 1:
                    raise AnalysisError(
                        f"{group_id}: N={step_count} has inconsistent time steps"
                    )

                for scheme in schemes:
                    candidates = scheme_points[scheme]
                    layouts = {_layout_signature(point) for point in candidates}
                    if len(layouts) != len(candidates):
                        raise AnalysisError(
                            f"{group_id}: {scheme}/N={step_count} has a duplicate layout"
                        )
                    serial = [point for point in candidates if _is_reference(point)]
                    if len(serial) != 1:
                        raise AnalysisError(
                            f"{group_id}: {scheme}/N={step_count} requires exactly one "
                            "MPI=1, OMP=1 reference"
                        )
                    if expected_layouts is None:
                        expected_layouts = layouts
                    elif layouts != expected_layouts:
                        raise AnalysisError(
                            f"{group_id}: {scheme}/N={step_count} has an incomplete layout matrix"
                        )
                    reference = serial[0]
                    references[scheme] = reference
                    try:
                        analytical = compare_fields(
                            reference.exact,
                            reference.numerical,
                            absolute_tolerance=ANALYTIC_ABSOLUTE_TOLERANCE,
                            relative_tolerance=ANALYTIC_RELATIVE_TOLERANCE,
                        )
                    except FieldComparisonError as error:
                        raise AnalysisError(
                            f"{group_id}: analytical field comparison failed: {error}"
                        ) from error
                    final_level = level_index == len(steps) - 1
                    analytical_accepted = analytical.passed if final_level else None
                    if final_level and not analytical.passed:
                        analytic_final_failures += 1
                        series_failed = True
                        diagnostics.append(
                            f"{group_id}: {reference.problem}/{scheme}/N={step_count} "
                            f"misses final analytical accuracy; {analytical.diagnostic()}"
                        )

                    variants: list[dict[str, object]] = []
                    for candidate in sorted(
                        (item for item in candidates if item is not reference),
                        key=_layout_signature,
                    ):
                        try:
                            comparison = compare_fields(
                                reference.numerical,
                                candidate.numerical,
                                absolute_tolerance=PARALLEL_ABSOLUTE_TOLERANCE,
                                relative_tolerance=PARALLEL_RELATIVE_TOLERANCE,
                            )
                        except FieldComparisonError as error:
                            raise AnalysisError(
                                f"{group_id}: parallel field comparison failed: {error}"
                            ) from error
                        parallel_comparisons += 1
                        timing_valid = comparison.passed and (
                            not final_level or analytical.passed
                        )
                        if not comparison.passed:
                            series_failed = True
                            diagnostics.append(
                                f"{group_id}: {candidate.problem}/{scheme} "
                                f"grid={candidate.mpi_grid} OMP={candidate.openmp_threads}; "
                                f"{comparison.diagnostic()}"
                            )
                        variants.append(
                            _case_document(candidate, timing_valid=timing_valid)
                            | {"comparison": _comparison_document(comparison)}
                        )

                    scheme_documents.append(
                        {
                            "scheme": scheme,
                            "reference": _case_document(
                                reference,
                                timing_valid=(not final_level or analytical.passed),
                            ),
                            "analytic_qualification": (
                                "final-accuracy" if final_level else "coarse-audit"
                            ),
                            "analytic_accepted": analytical_accepted,
                            "analytic_comparison": _comparison_document(analytical),
                            "variants": variants,
                        }
                    )

                pair_documents: list[dict[str, object]] = []
                for left, right in itertools.combinations(schemes, 2):
                    try:
                        comparison = compare_fields(
                            references[left].numerical,
                            references[right].numerical,
                            absolute_tolerance=SCHEME_ABSOLUTE_TOLERANCE,
                            relative_tolerance=SCHEME_RELATIVE_TOLERANCE,
                        )
                    except FieldComparisonError as error:
                        raise AnalysisError(
                            f"{group_id}: cross-scheme comparison failed: {error}"
                        ) from error
                    pair_history[left, right].append(comparison)
                    pair_documents.append(
                        {
                            "left_scheme": left,
                            "right_scheme": right,
                            "left_case_id": references[left].case_id,
                            "right_case_id": references[right].case_id,
                            "qualification": (
                                "final-accuracy"
                                if level_index == len(steps) - 1
                                else "coarse-audit"
                            ),
                            "comparison": _comparison_document(comparison),
                        }
                    )
                level_documents.append(
                    {
                        "steps": step_count,
                        "dt": float(next(iter(time_steps))),
                        "schemes": scheme_documents,
                        "scheme_pairs": pair_documents,
                    }
                )

            agreement_documents: list[dict[str, object]] = []
            schemes_failing_agreement: set[str] = set()
            for (left, right), comparisons in sorted(pair_history.items()):
                coarse = comparisons[0]
                final = comparisons[-1]
                l2_limit = coarse.l2_difference * (1.0 + SCHEME_TREND_RELAXATION)
                linf_limit = coarse.linf_difference * (1.0 + SCHEME_TREND_RELAXATION)
                non_growing = (
                    final.l2_difference <= l2_limit
                    and final.linf_difference <= linf_limit
                )
                identical_at_final = (
                    final.l2_difference == 0.0 and final.linf_difference == 0.0
                )
                meaningful_reduction = identical_at_final or (
                    final.l2_difference
                    <= coarse.l2_difference
                    * (1.0 - SCHEME_TREND_MINIMUM_REDUCTION)
                    or final.linf_difference
                    <= coarse.linf_difference
                    * (1.0 - SCHEME_TREND_MINIMUM_REDUCTION)
                )
                trend_passed = non_growing and meaningful_reduction
                status = "passed" if final.passed and trend_passed else "failed"
                if status == "failed":
                    schemes_failing_agreement.update((left, right))
                    scheme_pair_failures += 1
                    series_failed = True
                    diagnostics.append(
                        f"{group_id}: {left}/{right} scheme agreement failed; "
                        f"coarse L2/Linf={coarse.l2_difference:.8g}/"
                        f"{coarse.linf_difference:.8g}, final L2/Linf="
                        f"{final.l2_difference:.8g}/{final.linf_difference:.8g}; "
                        f"{final.diagnostic()}"
                    )
                agreement_documents.append(
                    {
                        "left_scheme": left,
                        "right_scheme": right,
                        "status": status,
                        "trend_passed": trend_passed,
                        "trend_relative_relaxation": SCHEME_TREND_RELAXATION,
                        "trend_minimum_relative_reduction": (
                            SCHEME_TREND_MINIMUM_REDUCTION
                        ),
                        "coarse_l2_difference": coarse.l2_difference,
                        "coarse_linf_difference": coarse.linf_difference,
                        "final_l2_difference": final.l2_difference,
                        "final_linf_difference": final.linf_difference,
                        "final_pointwise_accuracy_passed": final.passed,
                    }
                )

            # A finest-level time is eligible only after all numerical gates
            # relevant to that scheme have passed: analytical accuracy,
            # cross-scheme final accuracy/trend, and (for variants) equality
            # to the serial discrete reference.  Coarse analytical and
            # cross-scheme comparisons remain audits by design.
            final_level_document = level_documents[-1]
            final_scheme_documents = final_level_document["schemes"]
            assert isinstance(final_scheme_documents, list)
            for scheme_document in final_scheme_documents:
                assert isinstance(scheme_document, dict)
                scheme = scheme_document["scheme"]
                if scheme in schemes_failing_agreement:
                    reference_document = scheme_document["reference"]
                    assert isinstance(reference_document, dict)
                    reference_document["timing_valid"] = False
                    variant_documents = scheme_document["variants"]
                    assert isinstance(variant_documents, list)
                    for variant_document in variant_documents:
                        assert isinstance(variant_document, dict)
                        variant_document["timing_valid"] = False

            invalid_timings += sum(
                not bool(case_document["timing_valid"])
                for level_document in level_documents
                for scheme_document in level_document["schemes"]
                for case_document in (
                    scheme_document["reference"],
                    *scheme_document["variants"],
                )
            )

            status = "failed" if series_failed else "passed"
            if series_failed:
                failed_series += 1
            series_documents.append(
                {
                    "group_id": group_id,
                    "series_kind": "full-field-validation",
                    "problem": points[0].problem,
                    "final_time": float(points[0].final_time),
                    "controlled_configuration": dict(physics_documents[physics_key]),
                    "status": status,
                    "time_levels": level_documents,
                    "scheme_agreement": agreement_documents,
                }
            )

        series_documents.sort(
            key=lambda item: (str(item["problem"]), str(item["group_id"]))
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
                "parallel_comparison_count": parallel_comparisons,
                "invalid_timing_count": invalid_timings,
                "analytic_final_failure_count": analytic_final_failures,
                "scheme_pair_failure_count": scheme_pair_failures,
                "parallel_absolute_tolerance": PARALLEL_ABSOLUTE_TOLERANCE,
                "parallel_relative_tolerance": PARALLEL_RELATIVE_TOLERANCE,
                "analytic_absolute_tolerance": ANALYTIC_ABSOLUTE_TOLERANCE,
                "analytic_relative_tolerance": ANALYTIC_RELATIVE_TOLERANCE,
                "scheme_absolute_tolerance": SCHEME_ABSOLUTE_TOLERANCE,
                "scheme_relative_tolerance": SCHEME_RELATIVE_TOLERANCE,
                "scheme_trend_relative_relaxation": SCHEME_TREND_RELAXATION,
                "scheme_trend_minimum_relative_reduction": (
                    SCHEME_TREND_MINIMUM_REDUCTION
                ),
            },
            series=tuple(series_documents),
            diagnostics=tuple(diagnostics),
        )


__all__ = [
    "ANALYTIC_ABSOLUTE_TOLERANCE",
    "ANALYTIC_RELATIVE_TOLERANCE",
    "FieldValidationAnalyzer",
    "PARALLEL_ABSOLUTE_TOLERANCE",
    "PARALLEL_RELATIVE_TOLERANCE",
    "SCHEME_ABSOLUTE_TOLERANCE",
    "SCHEME_RELATIVE_TOLERANCE",
    "SCHEME_TREND_MINIMUM_REDUCTION",
    "SCHEME_TREND_RELAXATION",
    "VALIDATION_EXACT_CASE",
]
