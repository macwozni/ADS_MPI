"""Auditable spatial ``h`` and degree ``p`` convergence analysis.

Every spatial point is represented by two otherwise identical runs whose time
steps differ by exactly a factor of two.  The fine run supplies the spatial
error.  Time contamination is measured directly from the sampled field
difference ``u_dt - u_dt/2``; subtracting scalar error norms would only give a
lower bound and could hide cancellation.  A contaminated point is never used
in a local estimate or regression.
"""

from __future__ import annotations

from array import array
from collections.abc import Iterable, Mapping, Sequence
import csv
from dataclasses import dataclass, field, replace
from fractions import Fraction
import hashlib
import io
import json
import math
from ..framework.model import PlannedCase
from .base import AnalysisError, AnalysisReport


TEMPORAL_ERROR_FRACTION_LIMIT = 0.10
SPATIAL_EXACT_CASE = "spatial-cosine"
MINIMUM_H_LEVELS = 3
MINIMUM_P_LEVELS = 3
MINIMUM_P_PRE_PLATEAU_DECREASES = 2
MINIMUM_H_REGRESSION_R_SQUARED = 0.90
PLATEAU_REDUCTION_LIMIT = 1.25
PLATEAU_H_ORDER_LIMIT = 0.25
ROUNDOFF_RELATIVE_SCALE = 1.0e-10


Vector3 = tuple[int, int, int]


@dataclass(frozen=True)
class _FieldSamples:
    points_per_axis: int
    errors: array


@dataclass(frozen=True)
class _ResultPoint:
    case_id: str
    configuration: Mapping[str, object]
    family: str
    problem: str
    scheme: str
    exact_case: str
    final_time: Fraction
    time_step: Fraction
    steps: int
    mesh: Vector3
    test_degree: Vector3
    trial_degree: Vector3
    l2_error: float
    linf_error: float
    solution_l2_norm: float
    field_samples: _FieldSamples


@dataclass(frozen=True)
class _TemporalPair:
    coarse: _ResultPoint
    fine: _ResultPoint
    level_key: tuple[object, ...]
    level_document: Mapping[str, object]


def _mapping(value: object, field_name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise AnalysisError(f"{field_name} must be an object")
    return value


def _string(value: object, field_name: str) -> str:
    if not isinstance(value, str) or not value:
        raise AnalysisError(f"{field_name} must be a nonempty string")
    return value


def _integer(value: object, field_name: str) -> int:
    if type(value) is not int:
        raise AnalysisError(f"{field_name} must be an integer")
    return value


def _vector3(value: object, field_name: str) -> Vector3:
    if not isinstance(value, list) or len(value) != 3:
        raise AnalysisError(f"{field_name} must contain exactly three integers")
    result = tuple(
        _integer(component, f"{field_name}[{index}]")
        for index, component in enumerate(value)
    )
    if any(component <= 0 for component in result):
        raise AnalysisError(f"{field_name} entries must be positive")
    return result  # type: ignore[return-value]


def _finite(value: object, field_name: str, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise AnalysisError(f"{field_name} must be a finite number")
    result = float(value)
    if not math.isfinite(result):
        raise AnalysisError(f"{field_name} must be finite")
    if positive and result <= 0.0:
        raise AnalysisError(f"{field_name} must be positive")
    return result


def _sample_float(value: str, field_name: str) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise AnalysisError(f"{field_name} must be a finite number") from error
    if not math.isfinite(result):
        raise AnalysisError(f"{field_name} must be finite")
    return result


def _fraction(value: object, field_name: str) -> Fraction:
    if not isinstance(value, str):
        raise AnalysisError(f"{field_name} must be an exact-number string")
    try:
        result = Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise AnalysisError(
            f"{field_name} is not an exact number: {value!r}"
        ) from error
    if result <= 0:
        raise AnalysisError(f"{field_name} must be positive")
    return result


def _float_matches_exact(actual: float, expected: Fraction) -> bool:
    expected_float = float(expected)
    tolerance = max(64.0 * math.ulp(expected_float), 1.0e-14 * expected_float)
    return abs(actual - expected_float) <= tolerance


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


def _configuration_without(
    configuration: Mapping[str, object], *keys: str
) -> dict[str, object]:
    result = dict(configuration)
    for key in keys:
        result.pop(key, None)
    return result


def _extract_field_samples(
    document: Mapping[str, object],
    configuration: Mapping[str, object],
    domain: Mapping[str, object],
    case_id: str,
    final_time: Fraction,
    linf_error: float,
) -> _FieldSamples:
    sampling = _mapping(configuration.get("sampling"), f"{case_id}.sampling")
    points_per_axis = _integer(
        sampling.get("points_per_axis"), f"{case_id}.sampling.points_per_axis"
    )
    if points_per_axis < 2:
        raise AnalysisError(f"{case_id}: sampling requires at least two points")
    if sampling.get("write_samples") is not True:
        raise AnalysisError(
            f"{case_id}: spatial separation requires write_samples=true"
        )
    if (
        _integer(
            domain.get("sample_points_per_axis"),
            f"{case_id}.domain.sample_points_per_axis",
        )
        != points_per_axis
    ):
        raise AnalysisError(f"{case_id}: domain sampling size differs from plan")
    if domain.get("field_samples_written") is not True:
        raise AnalysisError(f"{case_id}: domain did not write field samples")

    artifacts = _mapping(
        document.get("analysis_artifacts"), f"{case_id}.analysis_artifacts"
    )
    artifact = artifacts.get("field_samples_csv")
    if isinstance(artifact, str):
        text = artifact
    else:
        read_text = getattr(artifact, "read_text", None)
        if not callable(read_text):
            raise AnalysisError(f"{case_id}: field_samples.csv is unavailable")
        text = read_text()
        if not isinstance(text, str):
            raise AnalysisError(f"{case_id}: field_samples.csv is unavailable")

    reader = csv.reader(io.StringIO(text, newline=""), strict=True)
    try:
        header = next(reader)
    except (StopIteration, csv.Error) as error:
        raise AnalysisError(f"{case_id}: field_samples.csv has no header") from error
    expected_header = ["x", "y", "z", "numerical", "exact", "error"]
    if header != expected_header:
        raise AnalysisError(
            f"{case_id}: field_samples.csv header must be "
            + ",".join(expected_header)
        )

    count = points_per_axis**3
    spacing = 1.0 / (points_per_axis - 1)
    physical_time = float(final_time)
    errors = array("d")
    observed_linf = 0.0
    sample_sum = 0.0
    sample_compensation = 0.0
    weighted_sum = 0.0
    weighted_compensation = 0.0
    try:
        for index, row in enumerate(reader):
            if index >= count:
                raise AnalysisError(
                    f"{case_id}: field_samples.csv has more than {count} rows"
                )
            if len(row) != 6:
                raise AnalysisError(
                    f"{case_id}: field sample row {index + 2} must have six fields"
                )
            values = tuple(
                _sample_float(value, f"{case_id}.field_samples[{index}][{column}]")
                for column, value in enumerate(row)
            )
            x, y, z, numerical, exact, error_value = values
            plane = points_per_axis * points_per_axis
            iz, remainder = divmod(index, plane)
            iy, ix = divmod(remainder, points_per_axis)
            expected_coordinates = (ix * spacing, iy * spacing, iz * spacing)
            if any(
                not math.isclose(
                    actual, expected, rel_tol=0.0, abs_tol=2.0e-14
                )
                for actual, expected in zip(
                    (x, y, z), expected_coordinates, strict=True
                )
            ):
                raise AnalysisError(
                    f"{case_id}: field sample row {index + 2} is not in the "
                    "declared regular-grid order"
                )
            expected_exact = (
                math.exp(-physical_time)
                * math.cos(math.pi * x)
                * math.cos(math.pi * y)
                * math.cos(math.pi * z)
            )
            if not math.isclose(
                exact, expected_exact, rel_tol=2.0e-13, abs_tol=2.0e-14
            ):
                raise AnalysisError(
                    f"{case_id}: field sample row {index + 2} has the wrong exact value"
                )
            if not math.isclose(
                error_value,
                numerical - exact,
                rel_tol=2.0e-11,
                abs_tol=5.0e-14,
            ):
                raise AnalysisError(
                    f"{case_id}: field sample row {index + 2} has an inconsistent error"
                )
            errors.append(error_value)
            observed_linf = max(observed_linf, abs(error_value))
            corrected = numerical - sample_compensation
            update = sample_sum + corrected
            sample_compensation = (update - sample_sum) - corrected
            sample_sum = update
            corrected = (index + 1) * numerical - weighted_compensation
            update = weighted_sum + corrected
            weighted_compensation = (update - weighted_sum) - corrected
            weighted_sum = update
    except csv.Error as error:
        raise AnalysisError(
            f"{case_id}: malformed field_samples.csv: {error}"
        ) from error

    if len(errors) != count:
        raise AnalysisError(
            f"{case_id}: field_samples.csv has {len(errors)} rows, expected {count}"
        )
    if not math.isclose(
        observed_linf, linf_error, rel_tol=2.0e-11, abs_tol=5.0e-14
    ):
        raise AnalysisError(
            f"{case_id}: sampled Linf error differs from the domain result"
        )
    observed_checksum = sample_sum + weighted_sum / (count + 1)
    expected_checksum = _finite(
        domain.get("field_checksum"), f"{case_id}.field_checksum"
    )
    if not math.isclose(
        observed_checksum,
        expected_checksum,
        rel_tol=2.0e-11,
        abs_tol=2.0e-12,
    ):
        raise AnalysisError(
            f"{case_id}: field sample checksum differs from the domain result"
        )
    return _FieldSamples(points_per_axis=points_per_axis, errors=errors)


def _extract_point(
    document: Mapping[str, object], expected_family: str
) -> _ResultPoint:
    case_id = _string(document.get("case_id"), "case_id")
    if document.get("kind") != "ads-benchmark-case-result":
        raise AnalysisError(f"{case_id}: unsupported result kind")
    if document.get("status") != "passed":
        raise AnalysisError(f"{case_id}: result status is not passed")

    configuration = _mapping(document.get("configuration"), f"{case_id}.configuration")
    family = _string(configuration.get("family"), f"{case_id}.family")
    if family != expected_family:
        raise AnalysisError(
            f"{case_id}: expected {expected_family} experiment family, got {family}"
        )
    problem = _string(configuration.get("problem"), f"{case_id}.problem")
    scheme = _string(configuration.get("scheme"), f"{case_id}.scheme").lower()
    exact_case = _string(configuration.get("exact_case"), f"{case_id}.exact_case")
    if exact_case != SPATIAL_EXACT_CASE:
        raise AnalysisError(
            f"{case_id}: spatial convergence requires exact case "
            f"{SPATIAL_EXACT_CASE!r}"
        )

    time = _mapping(configuration.get("time"), f"{case_id}.configuration.time")
    final_time = _fraction(time.get("final_time"), f"{case_id}.final_time")
    time_step = _fraction(time.get("time_step"), f"{case_id}.time_step")
    steps = _integer(time.get("steps"), f"{case_id}.steps")
    if steps <= 0 or time_step * steps != final_time:
        raise AnalysisError(f"{case_id}: inconsistent T, dt, and step count")

    mesh = _mapping(configuration.get("mesh"), f"{case_id}.mesh")
    elements = _vector3(mesh.get("elements"), f"{case_id}.mesh.elements")
    spaces = _mapping(configuration.get("spaces"), f"{case_id}.spaces")
    test_degree = _vector3(
        spaces.get("test_degree"), f"{case_id}.spaces.test_degree"
    )
    trial_degree = _vector3(
        spaces.get("trial_degree"), f"{case_id}.spaces.trial_degree"
    )

    domain = _mapping(document.get("domain_result"), f"{case_id}.domain_result")
    if domain.get("problem") != problem or domain.get("scheme") != scheme:
        raise AnalysisError(f"{case_id}: domain problem/scheme differs from plan")
    if domain.get("exact_case") != exact_case:
        raise AnalysisError(f"{case_id}: domain exact case differs from plan")
    if _integer(domain.get("solver_status"), f"{case_id}.solver_status") != 0:
        raise AnalysisError(f"{case_id}: solver_status is nonzero")
    if _integer(domain.get("steps"), f"{case_id}.domain.steps") != steps:
        raise AnalysisError(f"{case_id}: domain step count differs from plan")

    actual_final_time = _finite(
        domain.get("actual_final_time"), f"{case_id}.actual_final_time"
    )
    requested_final_time = _finite(
        domain.get("requested_final_time"), f"{case_id}.requested_final_time"
    )
    actual_time_step = _finite(
        domain.get("time_step"), f"{case_id}.domain.time_step"
    )
    if not _float_matches_exact(actual_final_time, final_time):
        raise AnalysisError(f"{case_id}: actual final time does not match T")
    if not _float_matches_exact(requested_final_time, final_time):
        raise AnalysisError(f"{case_id}: requested final time differs from plan")
    if not _float_matches_exact(actual_time_step, time_step):
        raise AnalysisError(f"{case_id}: domain time step differs from plan")

    l2_error = _finite(
        domain.get("l2_error"), f"{case_id}.l2_error", positive=True
    )
    linf_error = _finite(
        domain.get("linf_error"), f"{case_id}.linf_error", positive=True
    )
    solution_l2_norm = _finite(
        domain.get("solution_l2_norm"),
        f"{case_id}.solution_l2_norm",
        positive=True,
    )
    field_samples = _extract_field_samples(
        document,
        configuration,
        domain,
        case_id,
        final_time,
        linf_error,
    )

    return _ResultPoint(
        case_id=case_id,
        configuration=configuration,
        family=family,
        problem=problem,
        scheme=scheme,
        exact_case=exact_case,
        final_time=final_time,
        time_step=time_step,
        steps=steps,
        mesh=elements,
        test_degree=test_degree,
        trial_degree=trial_degree,
        l2_error=l2_error,
        linf_error=linf_error,
        solution_l2_norm=solution_l2_norm,
        field_samples=field_samples,
    )


def _degree_cohort(
    test_degree: Vector3, trial_degree: Vector3
) -> tuple[str, tuple[int, ...], Vector3]:
    enrichment = tuple(
        test - trial for test, trial in zip(test_degree, trial_degree, strict=True)
    )
    if len(set(enrichment)) != 1 or enrichment[0] not in (1, 2):
        raise AnalysisError(
            "p convergence requires component-wise p_test=p_trial+1 or +2"
        )
    enrichment_vector: Vector3 = enrichment  # type: ignore[assignment]
    if len(set(trial_degree)) == 1:
        return "isotropic-degree", (), enrichment_vector
    return "anisotropic-rotations", tuple(sorted(trial_degree)), enrichment_vector


def _series_descriptor(
    configuration: Mapping[str, object], family: str
) -> tuple[str, Mapping[str, object], Mapping[str, object]]:
    if family == "h":
        controlled = _configuration_without(configuration, "time", "mesh")
        cohort: dict[str, object] = {"sequence_kind": "h-refinement"}
        return (
            _canonical({"controlled": controlled, "cohort": cohort}),
            controlled,
            cohort,
        )

    spaces = _mapping(configuration.get("spaces"), "configuration.spaces")
    test_degree = _vector3(spaces.get("test_degree"), "spaces.test_degree")
    trial_degree = _vector3(spaces.get("trial_degree"), "spaces.trial_degree")
    sequence_kind, rotation_signature, enrichment = _degree_cohort(
        test_degree, trial_degree
    )
    cohort = {
        "sequence_kind": sequence_kind,
        "enrichment": list(enrichment),
        "rotation_signature": list(rotation_signature),
    }
    controlled = _configuration_without(configuration, "time", "spaces")
    return _canonical({"controlled": controlled, "cohort": cohort}), controlled, cohort


def _level_descriptor(
    point: _ResultPoint, family: str
) -> tuple[tuple[object, ...], Mapping[str, object]]:
    if family == "h":
        if len(set(point.mesh)) != 1:
            raise AnalysisError(
                f"{point.case_id}: h convergence requires isotropic element meshes"
            )
        return tuple(point.mesh), {"mesh": list(point.mesh), "h": 1.0 / point.mesh[0]}

    sequence_kind, _, enrichment = _degree_cohort(
        point.test_degree, point.trial_degree
    )
    return (
        (*point.trial_degree, *point.test_degree),
        {
            "sequence_kind": sequence_kind,
            "trial_degree": list(point.trial_degree),
            "test_degree": list(point.test_degree),
            "enrichment": list(enrichment),
        },
    )


def _pair_levels(
    points: Sequence[_ResultPoint], family: str, label: str
) -> tuple[_TemporalPair, ...]:
    levels: dict[tuple[object, ...], list[_ResultPoint]] = {}
    level_documents: dict[tuple[object, ...], Mapping[str, object]] = {}
    for point in points:
        level_key, level_document = _level_descriptor(point, family)
        levels.setdefault(level_key, []).append(point)
        level_documents[level_key] = level_document

    pairs: list[_TemporalPair] = []
    for level_key in sorted(levels):
        candidates = sorted(
            levels[level_key], key=lambda item: item.time_step, reverse=True
        )
        if len(candidates) != 2:
            raise AnalysisError(
                f"{label}: spatial level {list(level_key)} requires exactly two "
                f"temporal runs, got {len(candidates)}"
            )
        coarse, fine = candidates
        if coarse.final_time != fine.final_time:
            raise AnalysisError(
                f"{label}: temporal pair at {list(level_key)} has inconsistent T"
            )
        if coarse.time_step != 2 * fine.time_step:
            raise AnalysisError(
                f"{label}: temporal pair at {list(level_key)} must use dt and dt/2"
            )
        coarse_without_time = _configuration_without(coarse.configuration, "time")
        fine_without_time = _configuration_without(fine.configuration, "time")
        if _canonical(coarse_without_time) != _canonical(fine_without_time):
            raise AnalysisError(
                f"{label}: temporal pair differs in parameters other than time"
            )
        pairs.append(
            _TemporalPair(
                coarse=coarse,
                fine=fine,
                level_key=level_key,
                level_document=level_documents[level_key],
            )
        )
    return tuple(pairs)


def _temporal_field_difference(pair: _TemporalPair, metric_name: str) -> float:
    coarse_samples = pair.coarse.field_samples
    fine_samples = pair.fine.field_samples
    if coarse_samples.points_per_axis != fine_samples.points_per_axis:
        raise AnalysisError(
            f"{pair.coarse.case_id}/{pair.fine.case_id}: temporal pair uses "
            "different sample grids"
        )
    if len(coarse_samples.errors) != len(fine_samples.errors):
        raise AnalysisError(
            f"{pair.coarse.case_id}/{pair.fine.case_id}: temporal sample counts differ"
        )

    if metric_name == "linf":
        return max(
            abs(coarse - fine)
            for coarse, fine in zip(
                coarse_samples.errors, fine_samples.errors, strict=True
            )
        )
    if metric_name != "l2":
        raise AnalysisError(f"unsupported spatial metric: {metric_name}")

    points = coarse_samples.points_per_axis
    spacing = 1.0 / (points - 1)
    plane = points * points
    scale = 0.0
    scaled_square_sum = 1.0
    for index, (coarse, fine) in enumerate(
        zip(coarse_samples.errors, fine_samples.errors, strict=True)
    ):
        iz, remainder = divmod(index, plane)
        iy, ix = divmod(remainder, points)
        weight = 1.0
        if ix in (0, points - 1):
            weight *= 0.5
        if iy in (0, points - 1):
            weight *= 0.5
        if iz in (0, points - 1):
            weight *= 0.5
        component = math.sqrt(weight) * abs(coarse - fine)
        if not math.isfinite(component):
            return math.inf
        if component == 0.0:
            continue
        if scale < component:
            ratio = scale / component
            scaled_square_sum = 1.0 + scaled_square_sum * ratio * ratio
            scale = component
        else:
            ratio = component / scale
            scaled_square_sum += ratio * ratio

    if scale == 0.0:
        return 0.0
    return scale * math.sqrt(scaled_square_sum) * spacing**1.5


def _temporal_metric_audit(
    pair: _TemporalPair, metric_name: str, limit: float
) -> dict[str, object]:
    coarse_error = getattr(pair.coarse, f"{metric_name}_error")
    fine_error = getattr(pair.fine, f"{metric_name}_error")
    try:
        temporal_difference = _temporal_field_difference(pair, metric_name)
        temporal_fraction = temporal_difference / fine_error
    except OverflowError:
        temporal_difference = math.inf
        temporal_fraction = math.inf
    indicator_is_finite = math.isfinite(temporal_difference) and math.isfinite(
        temporal_fraction
    )
    reliable = indicator_is_finite and temporal_fraction <= limit
    return {
        "coarse_error": coarse_error,
        "fine_error": fine_error,
        "temporal_field_difference": (
            temporal_difference if math.isfinite(temporal_difference) else None
        ),
        "temporal_error_fraction": (
            temporal_fraction if math.isfinite(temporal_fraction) else None
        ),
        "temporal_indicator_status": (
            "finite" if indicator_is_finite else "nonfinite"
        ),
        "maximum_temporal_error_fraction": limit,
        "reliable_for_spatial_analysis": reliable,
        "sample_points_per_axis": pair.fine.field_samples.points_per_axis,
        "difference_method": (
            "same-grid-composite-trapezoid-l2"
            if metric_name == "l2"
            else "same-grid-pointwise-linf"
        ),
    }


def _linear_regression(xs: Sequence[float], ys: Sequence[float]) -> tuple[float, float]:
    if len(xs) != len(ys) or len(xs) < 2:
        raise AnalysisError("spatial regression requires at least two paired values")
    mean_x = math.fsum(xs) / len(xs)
    mean_y = math.fsum(ys) / len(ys)
    centered_x = [value - mean_x for value in xs]
    centered_y = [value - mean_y for value in ys]
    denominator = math.fsum(value * value for value in centered_x)
    if denominator <= 0.0:
        raise AnalysisError("spatial regression abscissas have zero variance")
    slope = math.fsum(
        left * right for left, right in zip(centered_x, centered_y, strict=True)
    ) / denominator
    total = math.fsum(value * value for value in centered_y)
    residual = math.fsum(
        (right - slope * left) ** 2
        for left, right in zip(centered_x, centered_y, strict=True)
    )
    r_squared = 1.0 if total == 0.0 else max(0.0, 1.0 - residual / total)
    if not math.isfinite(slope) or not math.isfinite(r_squared):
        raise AnalysisError("spatial regression produced a non-finite value")
    return slope, r_squared


def _plateau_start(
    pairs: Sequence[_TemporalPair],
    metric_audits: Sequence[Mapping[str, object]],
    family: str,
) -> int | None:
    if len(pairs) < 2:
        return None
    errors = [float(audit["fine_error"]) for audit in metric_audits]
    reliable = [bool(audit["reliable_for_spatial_analysis"]) for audit in metric_audits]
    low_improvement: list[bool] = []
    for index in range(len(pairs) - 1):
        if not (reliable[index] and reliable[index + 1]):
            low_improvement.append(False)
            continue
        if family == "h":
            coarse_h = float(pairs[index].level_document["h"])
            fine_h = float(pairs[index + 1].level_document["h"])
            order = math.log(errors[index] / errors[index + 1]) / math.log(
                coarse_h / fine_h
            )
            low_improvement.append(order <= PLATEAU_H_ORDER_LIMIT)
        else:
            low_improvement.append(
                errors[index] / errors[index + 1] <= PLATEAU_REDUCTION_LIMIT
            )

    suffix_count = 0
    for is_low in reversed(low_improvement):
        if not is_low:
            break
        suffix_count += 1
    if suffix_count == 0:
        return None

    transition = len(low_improvement) - suffix_count
    plateau_values = errors[transition:]
    band_ratio = max(plateau_values) / min(plateau_values)
    solution_scale = max(pair.fine.solution_l2_norm for pair in pairs)
    very_small = min(plateau_values) <= solution_scale * ROUNDOFF_RELATIVE_SCALE
    if (suffix_count >= 2 or very_small) and band_ratio <= PLATEAU_REDUCTION_LIMIT:
        return transition + 1
    return None


def _h_metric(
    pairs: Sequence[_TemporalPair], metric_name: str, limit: float
) -> tuple[dict[str, object], tuple[str, ...]]:
    audits = [_temporal_metric_audit(pair, metric_name, limit) for pair in pairs]
    errors = [float(audit["fine_error"]) for audit in audits]
    reliable = [bool(audit["reliable_for_spatial_analysis"]) for audit in audits]
    plateau_start = _plateau_start(pairs, audits, "h")

    local_orders: list[dict[str, object]] = []
    fatal_nonmonotonic: list[dict[str, object]] = []
    for index in range(len(pairs) - 1):
        coarse = pairs[index]
        fine = pairs[index + 1]
        endpoints_reliable = reliable[index] and reliable[index + 1]
        order = None
        if endpoints_reliable:
            order = math.log(errors[index] / errors[index + 1]) / math.log(
                float(coarse.level_document["h"])
                / float(fine.level_document["h"])
            )
            if not math.isfinite(order):
                raise AnalysisError(f"{metric_name}: non-finite h local order")
        in_plateau = plateau_start is not None and index + 1 >= plateau_start
        fatal = order is not None and order <= 0.0 and not in_plateau
        estimate = {
            "coarse_mesh": list(coarse.fine.mesh),
            "fine_mesh": list(fine.fine.mesh),
            "coarse_h": float(coarse.level_document["h"]),
            "fine_h": float(fine.level_document["h"]),
            "coarse_case_id": coarse.fine.case_id,
            "fine_case_id": fine.fine.case_id,
            "order": order,
            "available": endpoints_reliable,
            "unavailable_reason": (
                None if endpoints_reliable else "time-dominated endpoint"
            ),
            "in_plateau": in_plateau,
            "fatal_nonmonotonic": fatal,
            "used_in_regression": False,
        }
        local_orders.append(estimate)
        if fatal:
            fatal_nonmonotonic.append(dict(estimate))

    regression_indices = [
        index
        for index, is_reliable in enumerate(reliable)
        if is_reliable and (plateau_start is None or index < plateau_start)
    ]
    regression: dict[str, object] | None = None
    if len(regression_indices) >= MINIMUM_H_LEVELS:
        slope, r_squared = _linear_regression(
            [
                math.log(float(pairs[index].level_document["h"]))
                for index in regression_indices
            ],
            [math.log(errors[index]) for index in regression_indices],
        )
        regression = {
            "method": "centered-ordinary-least-squares-log-error-vs-log-h",
            "point_count": len(regression_indices),
            "point_case_ids": [
                pairs[index].fine.case_id for index in regression_indices
            ],
            "order": slope,
            "r_squared": r_squared,
            "minimum_r_squared": MINIMUM_H_REGRESSION_R_SQUARED,
        }
        selected = set(regression_indices)
        for index, estimate in enumerate(local_orders):
            estimate["used_in_regression"] = index in selected and index + 1 in selected

    time_dominated = [
        {
            "mesh": list(pair.fine.mesh),
            "fine_case_id": pair.fine.case_id,
            "temporal_error_fraction": audit["temporal_error_fraction"],
        }
        for pair, audit in zip(pairs, audits, strict=True)
        if not audit["reliable_for_spatial_analysis"]
    ]
    accepted = (
        not time_dominated
        and not fatal_nonmonotonic
        and regression is not None
        and float(regression["order"]) > 0.0
        and float(regression["r_squared"]) >= MINIMUM_H_REGRESSION_R_SQUARED
    )
    status = "passed" if accepted else ("unreliable" if time_dominated else "failed")
    diagnostics: list[str] = []
    if time_dominated:
        diagnostics.append(
            f"{metric_name}: {len(time_dominated)} time-dominated h point(s) excluded"
        )
    if plateau_start is not None:
        diagnostics.append(
            f"{metric_name}: plateau begins at mesh "
            f"{list(pairs[plateau_start].fine.mesh)}"
        )
    if fatal_nonmonotonic:
        diagnostics.append(
            f"{metric_name}: nonmonotonic h error outside a detected plateau"
        )
    if regression is None:
        diagnostics.append(
            f"{metric_name}: fewer than {MINIMUM_H_LEVELS} reliable pre-plateau "
            "points remain for regression"
        )
    elif float(regression["r_squared"]) < MINIMUM_H_REGRESSION_R_SQUARED:
        diagnostics.append(
            f"{metric_name}: h regression R^2 is below "
            f"{MINIMUM_H_REGRESSION_R_SQUARED:g}"
        )

    return (
        {
            "status": status,
            "acceptance_mode": "positive-h-regression",
            "acceptance_reason": (
                "requires a positive qualified log-error-vs-log-h regression"
            ),
            "temporal_fraction_limit": limit,
            "time_dominated_points": time_dominated,
            "plateau_detected": plateau_start is not None,
            "plateau_start_mesh": (
                None if plateau_start is None else list(pairs[plateau_start].fine.mesh)
            ),
            "regression": regression,
            "fatal_nonmonotonic_anomalies": fatal_nonmonotonic,
            "local_orders": local_orders,
            "rotation_error_spread": None,
        },
        tuple(diagnostics),
    )


def _p_metric(
    pairs: Sequence[_TemporalPair],
    metric_name: str,
    limit: float,
    sequence_kind: str,
) -> tuple[dict[str, object], tuple[str, ...]]:
    audits = [_temporal_metric_audit(pair, metric_name, limit) for pair in pairs]
    errors = [float(audit["fine_error"]) for audit in audits]
    reliable = [bool(audit["reliable_for_spatial_analysis"]) for audit in audits]
    plateau_start = (
        _plateau_start(pairs, audits, "p")
        if sequence_kind == "isotropic-degree"
        else None
    )
    ratios: list[dict[str, object]] = []
    fatal_nonmonotonic: list[dict[str, object]] = []
    if sequence_kind == "isotropic-degree":
        for index in range(len(pairs) - 1):
            coarse = pairs[index]
            fine = pairs[index + 1]
            endpoints_reliable = reliable[index] and reliable[index + 1]
            ratio = errors[index] / errors[index + 1] if endpoints_reliable else None
            in_plateau = plateau_start is not None and index + 1 >= plateau_start
            fatal = ratio is not None and ratio <= 1.0 and not in_plateau
            estimate = {
                "coarse_trial_degree": list(coarse.fine.trial_degree),
                "fine_trial_degree": list(fine.fine.trial_degree),
                "coarse_test_degree": list(coarse.fine.test_degree),
                "fine_test_degree": list(fine.fine.test_degree),
                "coarse_case_id": coarse.fine.case_id,
                "fine_case_id": fine.fine.case_id,
                "error_ratio": ratio,
                "error_decreased": None if ratio is None else ratio > 1.0,
                "available": endpoints_reliable,
                "unavailable_reason": (
                    None if endpoints_reliable else "time-dominated endpoint"
                ),
                "in_plateau": in_plateau,
                "fatal_nonmonotonic": fatal,
            }
            ratios.append(estimate)
            if fatal:
                fatal_nonmonotonic.append(dict(estimate))

    pre_plateau_decreases = sum(
        1
        for estimate in ratios
        if estimate["available"]
        and not estimate["in_plateau"]
        and float(estimate["error_ratio"]) > 1.0
    )

    time_dominated = [
        {
            "trial_degree": list(pair.fine.trial_degree),
            "test_degree": list(pair.fine.test_degree),
            "fine_case_id": pair.fine.case_id,
            "temporal_error_fraction": audit["temporal_error_fraction"],
        }
        for pair, audit in zip(pairs, audits, strict=True)
        if not audit["reliable_for_spatial_analysis"]
    ]
    rotation_spread = None
    if sequence_kind == "anisotropic-rotations" and all(reliable):
        rotation_spread = max(errors) / min(errors)

    sufficient_decrease = (
        sequence_kind != "isotropic-degree"
        or pre_plateau_decreases >= MINIMUM_P_PRE_PLATEAU_DECREASES
    )
    accepted = not time_dominated and not fatal_nonmonotonic and sufficient_decrease
    status = "passed" if accepted else ("unreliable" if time_dominated else "failed")
    diagnostics: list[str] = []
    if time_dominated:
        diagnostics.append(
            f"{metric_name}: {len(time_dominated)} time-dominated p point(s) excluded"
        )
    if plateau_start is not None:
        diagnostics.append(
            f"{metric_name}: p plateau begins at trial degree "
            f"{list(pairs[plateau_start].fine.trial_degree)}"
        )
    if fatal_nonmonotonic:
        diagnostics.append(
            f"{metric_name}: p error does not decrease outside a detected plateau"
        )
    if not sufficient_decrease:
        diagnostics.append(
            f"{metric_name}: fewer than {MINIMUM_P_PRE_PLATEAU_DECREASES} "
            "strict pre-plateau error decreases"
        )

    return (
        {
            "status": status,
            "acceptance_mode": (
                "ordered-error-decrease"
                if sequence_kind == "isotropic-degree"
                else "audit-only"
            ),
            "acceptance_reason": (
                "requires qualified pre-plateau error decreases"
                if sequence_kind == "isotropic-degree"
                else "fixed-axis ADI rotations have no prescribed equality tolerance"
            ),
            "temporal_fraction_limit": limit,
            "time_dominated_points": time_dominated,
            "plateau_detected": plateau_start is not None,
            "plateau_start_trial_degree": (
                None
                if plateau_start is None
                else list(pairs[plateau_start].fine.trial_degree)
            ),
            # This is deliberately a sequence of raw reductions, not a
            # polynomial order inferred from degree increments.
            "error_ratios": ratios,
            "strict_pre_plateau_decrease_count": pre_plateau_decreases,
            "minimum_strict_pre_plateau_decreases": (
                MINIMUM_P_PRE_PLATEAU_DECREASES
                if sequence_kind == "isotropic-degree"
                else None
            ),
            "fatal_nonmonotonic_anomalies": fatal_nonmonotonic,
            "rotation_error_spread": rotation_spread,
            "regression": None,
        },
        tuple(diagnostics),
    )


def _point_document(
    pair: _TemporalPair, limit: float
) -> dict[str, object]:
    l2 = _temporal_metric_audit(pair, "l2", limit)
    linf = _temporal_metric_audit(pair, "linf", limit)
    return {
        **dict(pair.level_document),
        "coarse_case_id": pair.coarse.case_id,
        "fine_case_id": pair.fine.case_id,
        "coarse_steps": pair.coarse.steps,
        "fine_steps": pair.fine.steps,
        "coarse_dt": float(pair.coarse.time_step),
        "fine_dt": float(pair.fine.time_step),
        "l2_error": pair.fine.l2_error,
        "linf_error": pair.fine.linf_error,
        "temporal_pair": {"l2": l2, "linf": linf},
    }


@dataclass(frozen=True)
class SpatialConvergenceAnalyzer:
    """Analyze either an ``h`` or a ``p`` family without mixing controls."""

    family: str = "h"
    name: str = "h-convergence"
    temporal_fraction_limit: float = TEMPORAL_ERROR_FRACTION_LIMIT
    expected_configurations: tuple[str, ...] | None = None

    def __post_init__(self) -> None:
        if self.family not in {"h", "p"}:
            raise ValueError("spatial analyzer family must be 'h' or 'p'")
        if not (0.0 < self.temporal_fraction_limit < 1.0):
            raise ValueError(
                "temporal fraction limit must lie strictly between 0 and 1"
            )

    def configure(
        self, planned_cases: Sequence[PlannedCase]
    ) -> "SpatialConvergenceAnalyzer":
        if not planned_cases:
            raise AnalysisError("cannot configure spatial analyzer from an empty plan")
        configurations: list[str] = []
        for case in planned_cases:
            if case.spec.family != self.family:
                raise AnalysisError(
                    f"analyzer {self.name} cannot analyze family {case.spec.family}"
                )
            configurations.append(_canonical(case.spec.to_dict()))
        if len(set(configurations)) != len(configurations):
            raise AnalysisError("frozen spatial plan contains duplicate configurations")
        return replace(self, expected_configurations=tuple(sorted(configurations)))

    def analyze(
        self,
        results: Iterable[Mapping[str, object]],
        *,
        source_run: str | None = None,
    ) -> AnalysisReport:
        documents = tuple(results)
        if not documents:
            raise AnalysisError(f"{self.family} analysis received no results")

        seen_ids: set[str] = set()
        seen_configurations: set[str] = set()
        groups: dict[str, list[_ResultPoint]] = {}
        controlled_configurations: dict[str, Mapping[str, object]] = {}
        cohorts: dict[str, Mapping[str, object]] = {}
        for index, document in enumerate(documents):
            if not isinstance(document, Mapping):
                raise AnalysisError(f"result {index} is not an object")
            point = _extract_point(document, self.family)
            if point.case_id in seen_ids:
                raise AnalysisError(f"duplicate case_id: {point.case_id}")
            seen_ids.add(point.case_id)
            configuration_key = _canonical(point.configuration)
            if configuration_key in seen_configurations:
                raise AnalysisError(
                    f"duplicate spatial configuration at case {point.case_id}"
                )
            seen_configurations.add(configuration_key)
            group_key, controlled, cohort = _series_descriptor(
                point.configuration, self.family
            )
            groups.setdefault(group_key, []).append(point)
            controlled_configurations[group_key] = controlled
            cohorts[group_key] = cohort

        if self.expected_configurations is not None:
            actual = tuple(sorted(seen_configurations))
            if actual != self.expected_configurations:
                expected_set = set(self.expected_configurations)
                actual_set = set(actual)

                def labels(values: set[str]) -> list[str]:
                    return sorted(
                        hashlib.sha256(value.encode("utf-8")).hexdigest()[:16]
                        for value in values
                    )

                raise AnalysisError(
                    "results do not match frozen spatial plan; "
                    f"missing={labels(expected_set - actual_set)}, "
                    f"unexpected={labels(actual_set - expected_set)}"
                )

        series_reports: list[dict[str, object]] = []
        diagnostics: list[str] = []
        failed_series = 0
        unreliable_metric_points = 0
        plateau_metrics = 0
        for group_key in sorted(groups):
            points = groups[group_key]
            label = f"{points[0].problem}/{points[0].scheme}"
            pairs = list(_pair_levels(points, self.family, label))
            if len({pair.coarse.final_time for pair in pairs}) != 1:
                raise AnalysisError(
                    f"{label}: inconsistent final time across spatial levels"
                )

            sequence_kind = str(cohorts[group_key]["sequence_kind"])
            if self.family == "h":
                pairs.sort(
                    key=lambda pair: float(pair.level_document["h"]), reverse=True
                )
                if len(pairs) < MINIMUM_H_LEVELS:
                    raise AnalysisError(
                        f"{label}: h convergence requires at least "
                        f"{MINIMUM_H_LEVELS} levels"
                    )
                if any(
                    float(pairs[index].level_document["h"])
                    <= float(pairs[index + 1].level_document["h"])
                    for index in range(len(pairs) - 1)
                ):
                    raise AnalysisError(
                        f"{label}: h levels are not strictly decreasing"
                    )
                l2, l2_diagnostics = _h_metric(
                    pairs, "l2", self.temporal_fraction_limit
                )
                linf, linf_diagnostics = _h_metric(
                    pairs, "linf", self.temporal_fraction_limit
                )
            else:
                if sequence_kind == "isotropic-degree":
                    pairs.sort(key=lambda pair: pair.fine.trial_degree[0])
                    if len(pairs) < MINIMUM_P_LEVELS:
                        raise AnalysisError(
                            f"{label}: p convergence requires at least "
                            f"{MINIMUM_P_LEVELS} degrees"
                        )
                    scalar_degrees = [pair.fine.trial_degree[0] for pair in pairs]
                    if len(set(scalar_degrees)) != len(scalar_degrees):
                        raise AnalysisError(
                            f"{label}: duplicate isotropic trial degree"
                        )
                else:
                    pairs.sort(key=lambda pair: pair.fine.trial_degree)
                    if len(pairs) < 2:
                        raise AnalysisError(
                            f"{label}: anisotropic audit requires at least two "
                            "rotations"
                        )
                l2, l2_diagnostics = _p_metric(
                    pairs,
                    "l2",
                    self.temporal_fraction_limit,
                    sequence_kind,
                )
                linf, linf_diagnostics = _p_metric(
                    pairs,
                    "linf",
                    self.temporal_fraction_limit,
                    sequence_kind,
                )

            status = (
                "passed"
                if l2["status"] == "passed" and linf["status"] == "passed"
                else "failed"
            )
            if status == "failed":
                failed_series += 1
            unreliable_metric_points += len(l2["time_dominated_points"]) + len(
                linf["time_dominated_points"]
            )
            plateau_metrics += int(bool(l2["plateau_detected"])) + int(
                bool(linf["plateau_detected"])
            )
            group_id = hashlib.sha256(group_key.encode("utf-8")).hexdigest()[:16]
            for message in (*l2_diagnostics, *linf_diagnostics):
                diagnostics.append(f"{group_id}: {message}")
            series_reports.append(
                {
                    "group_id": group_id,
                    "family": self.family,
                    "sequence_kind": sequence_kind,
                    "acceptance_mode": (
                        "ordered-error-decrease"
                        if sequence_kind == "isotropic-degree"
                        else (
                            "audit-only"
                            if sequence_kind == "anisotropic-rotations"
                            else "positive-h-regression"
                        )
                    ),
                    "problem": points[0].problem,
                    "scheme": points[0].scheme,
                    "final_time": float(points[0].final_time),
                    "controlled_configuration": dict(
                        controlled_configurations[group_key]
                    ),
                    "degree_cohort": dict(cohorts[group_key]),
                    "status": status,
                    "points": [
                        _point_document(pair, self.temporal_fraction_limit)
                        for pair in pairs
                    ],
                    "metrics": {"l2": l2, "linf": linf},
                }
            )

        series_reports.sort(
            key=lambda item: (
                str(item["problem"]),
                str(item["scheme"]),
                str(item["sequence_kind"]),
                str(item["group_id"]),
            )
        )
        return AnalysisReport(
            analyzer=self.name,
            family=self.family,
            status="passed" if failed_series == 0 else "failed",
            source_run=source_run,
            summary={
                "series_count": len(series_reports),
                "passed_series": len(series_reports) - failed_series,
                "failed_series": failed_series,
                "spatial_point_count": sum(
                    len(series["points"]) for series in series_reports
                ),
                "unreliable_metric_points": unreliable_metric_points,
                "plateau_metric_count": plateau_metrics,
                "temporal_error_fraction_limit": self.temporal_fraction_limit,
                "sequence_kinds": sorted(
                    {str(series["sequence_kind"]) for series in series_reports}
                ),
            },
            series=tuple(series_reports),
            diagnostics=tuple(diagnostics),
        )


@dataclass(frozen=True)
class HConvergenceAnalyzer(SpatialConvergenceAnalyzer):
    family: str = field(default="h", init=False)
    name: str = field(default="h-convergence", init=False)


@dataclass(frozen=True)
class PConvergenceAnalyzer(SpatialConvergenceAnalyzer):
    family: str = field(default="p", init=False)
    name: str = field(default="p-convergence", init=False)


__all__ = [
    "HConvergenceAnalyzer",
    "PConvergenceAnalyzer",
    "SPATIAL_EXACT_CASE",
    "SpatialConvergenceAnalyzer",
    "TEMPORAL_ERROR_FRACTION_LIMIT",
]
