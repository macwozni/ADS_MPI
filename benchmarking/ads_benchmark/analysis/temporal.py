"""Strict temporal-order analysis against manufactured analytical errors."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from fractions import Fraction
import hashlib
import json
import math
from typing import Any

from ..framework.model import PlannedCase
from .base import AnalysisError, AnalysisReport


EXPECTED_TEMPORAL_STEPS = (4, 8, 16, 32, 64, 128, 256, 512)
MINIMUM_ANALYTICAL_REDUCTION = 4.0
MINIMUM_REGRESSION_R_SQUARED = 0.98
REGRESSION_POINT_COUNT = 4
# A refinement jump exceeding four times the formal order is not plausible
# asymptotic behavior.  In particular, halving dt cannot credibly reduce a
# first-order error by more than 2**4=16 while the same series later follows a
# clean first-order tail.  Such a jump is retained as evidence of a defective
# coarse solve; it is never discarded as pre-asymptotic data.
MAXIMUM_LOCAL_ORDER_MULTIPLIER = 4.0


@dataclass(frozen=True)
class SchemeExpectation:
    theoretical_order: float
    accepted_minimum: float
    accepted_maximum: float


# DG is the second-order midpoint update used by this repository.  PR and BE
# are first-order for the currently implemented split/backward updates.  The
# intervals allow modest pre-asymptotic contamination without accepting a
# different formal order.
SCHEME_EXPECTATIONS: Mapping[str, SchemeExpectation] = {
    "dg": SchemeExpectation(2.0, 1.70, 2.30),
    "pr": SchemeExpectation(1.0, 0.80, 1.20),
    "be": SchemeExpectation(1.0, 0.80, 1.20),
}


@dataclass(frozen=True)
class _Point:
    case_id: str
    problem: str
    scheme: str
    final_time: Fraction
    time_step: Fraction
    steps: int
    l2_error: float
    linf_error: float
    solution_l2_norm: float
    configuration_without_time: Mapping[str, object]


def _mapping(value: object, field: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise AnalysisError(f"{field} must be an object")
    return value


def _string(value: object, field: str) -> str:
    if not isinstance(value, str) or not value:
        raise AnalysisError(f"{field} must be a nonempty string")
    return value


def _integer(value: object, field: str) -> int:
    if type(value) is not int:
        raise AnalysisError(f"{field} must be an integer")
    return value


def _finite(value: object, field: str, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise AnalysisError(f"{field} must be a finite number")
    normalized = float(value)
    if not math.isfinite(normalized):
        raise AnalysisError(f"{field} must be finite")
    if positive and normalized <= 0.0:
        raise AnalysisError(f"{field} must be positive")
    return normalized


def _fraction(value: object, field: str) -> Fraction:
    if not isinstance(value, str):
        raise AnalysisError(f"{field} must be an exact-number string")
    try:
        normalized = Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise AnalysisError(f"{field} is not an exact number: {value!r}") from error
    if normalized <= 0:
        raise AnalysisError(f"{field} must be positive")
    return normalized


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


def _extract_point(document: Mapping[str, object]) -> _Point:
    case_id = _string(document.get("case_id"), "case_id")
    if document.get("kind") != "ads-benchmark-case-result":
        raise AnalysisError(f"{case_id}: unsupported result kind")
    if document.get("status") != "passed":
        raise AnalysisError(f"{case_id}: result status is not passed")

    configuration = _mapping(document.get("configuration"), f"{case_id}.configuration")
    if configuration.get("family") != "temporal":
        raise AnalysisError(f"{case_id}: expected temporal experiment family")
    problem = _string(configuration.get("problem"), f"{case_id}.problem")
    scheme = _string(configuration.get("scheme"), f"{case_id}.scheme").lower()
    if scheme not in SCHEME_EXPECTATIONS:
        raise AnalysisError(f"{case_id}: unsupported temporal scheme {scheme!r}")
    exact_case = _string(configuration.get("exact_case"), f"{case_id}.exact_case")

    time = _mapping(configuration.get("time"), f"{case_id}.configuration.time")
    final_time = _fraction(time.get("final_time"), f"{case_id}.final_time")
    time_step = _fraction(time.get("time_step"), f"{case_id}.time_step")
    steps = _integer(time.get("steps"), f"{case_id}.steps")
    if steps <= 0 or time_step * steps != final_time:
        raise AnalysisError(f"{case_id}: inconsistent T, dt, and step count")

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
    if not _float_matches_exact(actual_final_time, final_time):
        raise AnalysisError(
            f"{case_id}: actual final time {actual_final_time:.17g} "
            f"does not match T={float(final_time):.17g}"
        )
    requested_final_time = _finite(
        domain.get("requested_final_time"), f"{case_id}.requested_final_time"
    )
    if not _float_matches_exact(requested_final_time, final_time):
        raise AnalysisError(f"{case_id}: requested final time differs from plan")
    actual_time_step = _finite(domain.get("time_step"), f"{case_id}.domain.time_step")
    if not _float_matches_exact(actual_time_step, time_step):
        raise AnalysisError(f"{case_id}: domain time step differs from plan")

    configuration_without_time = dict(configuration)
    configuration_without_time.pop("time", None)
    return _Point(
        case_id=case_id,
        problem=problem,
        scheme=scheme,
        final_time=final_time,
        time_step=time_step,
        steps=steps,
        # These are errors against the registered analytical solution, not
        # differences between two numerical refinements.
        l2_error=_finite(domain.get("l2_error"), f"{case_id}.l2_error", positive=True),
        linf_error=_finite(
            domain.get("linf_error"), f"{case_id}.linf_error", positive=True
        ),
        solution_l2_norm=_finite(
            domain.get("solution_l2_norm"),
            f"{case_id}.solution_l2_norm",
            positive=True,
        ),
        configuration_without_time=configuration_without_time,
    )


def _linear_regression(xs: Sequence[float], ys: Sequence[float]) -> tuple[float, float]:
    if len(xs) != len(ys) or len(xs) < 2:
        raise AnalysisError("regression requires at least two paired values")
    mean_x = math.fsum(xs) / len(xs)
    mean_y = math.fsum(ys) / len(ys)
    centered_x = [value - mean_x for value in xs]
    centered_y = [value - mean_y for value in ys]
    denominator = math.fsum(value * value for value in centered_x)
    if denominator <= 0.0:
        raise AnalysisError("regression time steps have zero variance")
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
        raise AnalysisError("regression produced a non-finite value")
    return slope, r_squared


def _plateau_transition(
    local_orders: Sequence[float],
    errors: Sequence[float],
    solution_norms: Sequence[float],
    expectation: SchemeExpectation,
) -> int | None:
    """Return the first transition in a credible finest-level plateau."""

    low_order_limit = max(0.20, 0.25 * expectation.theoretical_order)
    suffix_count = 0
    for order in reversed(local_orders):
        if order <= low_order_limit:
            suffix_count += 1
        else:
            break
    if suffix_count == 0:
        return None

    transition = len(local_orders) - suffix_count
    plateau_values = errors[transition:]
    band_ratio = max(plateau_values) / min(plateau_values)
    very_small = min(plateau_values) <= max(solution_norms) * 1.0e-10
    credible_length = suffix_count >= 2
    # One final flat transition is accepted only at an actual roundoff scale;
    # two or more transitions also model a solver-tolerance plateau.  A large
    # late error spike is not mislabeled as a plateau.
    if (credible_length or very_small) and band_ratio <= 1.25:
        return transition
    return None


def _analyze_metric(
    points: Sequence[_Point],
    metric_name: str,
    expectation: SchemeExpectation,
) -> tuple[dict[str, object], tuple[str, ...]]:
    errors = [getattr(point, f"{metric_name}_error") for point in points]
    dts = [float(point.time_step) for point in points]
    local_orders = [
        math.log(errors[index] / errors[index + 1])
        / math.log(dts[index] / dts[index + 1])
        for index in range(len(points) - 1)
    ]
    if any(not math.isfinite(value) for value in local_orders):
        raise AnalysisError(f"{metric_name}: local order is non-finite")

    plateau_transition = _plateau_transition(
        local_orders,
        errors,
        [point.solution_l2_norm for point in points],
        expectation,
    )
    usable_count = (
        len(points) if plateau_transition is None else plateau_transition + 1
    )
    if usable_count < REGRESSION_POINT_COUNT:
        raise AnalysisError(
            f"{metric_name}: plateau leaves only {usable_count} usable points"
        )

    monotonic_violations = [
        index for index, order in enumerate(local_orders) if order <= 0.0
    ]
    fatal_violations = [
        index
        for index in monotonic_violations
        if plateau_transition is None or index < plateau_transition
    ]
    regression_start = usable_count - REGRESSION_POINT_COUNT
    regression_points = points[regression_start:usable_count]
    regression_errors = errors[regression_start:usable_count]
    slope, r_squared = _linear_regression(
        [math.log(float(point.time_step)) for point in regression_points],
        [math.log(value) for value in regression_errors],
    )
    analytical_reduction = errors[0] / errors[usable_count - 1]
    maximum_local_order = (
        MAXIMUM_LOCAL_ORDER_MULTIPLIER * expectation.theoretical_order
    )
    implausible_orders = [
        index
        for index, order in enumerate(local_orders)
        if order > maximum_local_order
    ]
    accepted = (
        expectation.accepted_minimum
        <= slope
        <= expectation.accepted_maximum
        and r_squared >= MINIMUM_REGRESSION_R_SQUARED
        and analytical_reduction >= MINIMUM_ANALYTICAL_REDUCTION
        and not implausible_orders
        and not fatal_violations
    )

    estimates: list[dict[str, object]] = []
    for index, order in enumerate(local_orders):
        estimates.append(
            {
                "coarse_steps": points[index].steps,
                "fine_steps": points[index + 1].steps,
                "coarse_dt": dts[index],
                "fine_dt": dts[index + 1],
                "order": order,
                "used_in_regression": (
                    regression_start <= index < usable_count - 1
                ),
                "in_plateau": (
                    plateau_transition is not None and index >= plateau_transition
                ),
                "nonmonotonic": index in monotonic_violations,
                "fatal_nonmonotonic": index in fatal_violations,
                "implausible": index in implausible_orders,
            }
        )

    diagnostics: list[str] = []
    if plateau_transition is not None:
        diagnostics.append(
            f"{metric_name}: excluded plateau beginning at "
            f"N={points[plateau_transition + 1].steps}"
        )
    if fatal_violations:
        labels = ", ".join(
            f"N={points[index].steps}->{points[index + 1].steps}"
            for index in fatal_violations
        )
        diagnostics.append(
            f"{metric_name}: fatal nonmonotonic analytical error at {labels}"
        )
    elif monotonic_violations:
        diagnostics.append(
            f"{metric_name}: nonmonotonic transitions confined to detected plateau"
        )
    if implausible_orders:
        labels = ", ".join(
            f"N={points[index].steps}->{points[index + 1].steps} "
            f"(p={local_orders[index]:.6g})"
            for index in implausible_orders
        )
        diagnostics.append(
            f"{metric_name}: implausible local order above "
            f"{maximum_local_order:g} at {labels}"
        )
    if not accepted:
        diagnostics.append(
            f"{metric_name}: convergence checks failed with regression order "
            f"{slope:.6g} and analytical reduction "
            f"{analytical_reduction:.6g}"
        )

    return (
        {
            "status": "passed" if accepted else "failed",
            "expected_order": expectation.theoretical_order,
            "accepted_interval": [
                expectation.accepted_minimum,
                expectation.accepted_maximum,
            ],
            "coarse_error": errors[0],
            "finest_error": errors[-1],
            "finest_usable_error": errors[usable_count - 1],
            "analytical_reduction": analytical_reduction,
            "maximum_local_order": maximum_local_order,
            "implausible_local_orders": [
                {
                    "coarse_steps": points[index].steps,
                    "fine_steps": points[index + 1].steps,
                    "order": local_orders[index],
                }
                for index in implausible_orders
            ],
            "plateau_detected": plateau_transition is not None,
            "plateau_start_steps": (
                None
                if plateau_transition is None
                else points[plateau_transition + 1].steps
            ),
            "monotonic_violations": [
                {
                    "coarse_steps": points[index].steps,
                    "fine_steps": points[index + 1].steps,
                    "order": local_orders[index],
                }
                for index in monotonic_violations
            ],
            "fatal_nonmonotonic_anomalies": [
                {
                    "coarse_steps": points[index].steps,
                    "fine_steps": points[index + 1].steps,
                    "order": local_orders[index],
                }
                for index in fatal_violations
            ],
            "regression": {
                "method": "centered-ordinary-least-squares-log-log",
                "point_count": REGRESSION_POINT_COUNT,
                "first_steps": regression_points[0].steps,
                "last_steps": regression_points[-1].steps,
                "order": slope,
                "r_squared": r_squared,
                "minimum_r_squared": MINIMUM_REGRESSION_R_SQUARED,
            },
            "local_orders": estimates,
        },
        tuple(diagnostics),
    )


@dataclass(frozen=True)
class TemporalConvergenceAnalyzer:
    """Analyze complete fixed-space temporal series of at least four levels."""

    expected_steps: tuple[int, ...] = EXPECTED_TEMPORAL_STEPS
    expected_groups: tuple[str, ...] | None = None
    name: str = "temporal-convergence"
    family: str = "temporal"

    def __post_init__(self) -> None:
        if (
            len(self.expected_steps) < REGRESSION_POINT_COUNT
            or any(type(item) is not int or item <= 0 for item in self.expected_steps)
            or tuple(sorted(set(self.expected_steps))) != self.expected_steps
        ):
            raise ValueError(
                "temporal analyzer requires at least four unique increasing steps"
            )

    def configure(
        self, planned_cases: Sequence[PlannedCase]
    ) -> "TemporalConvergenceAnalyzer":
        """Derive intentional refinement levels from a validated frozen plan.

        A bare analyzer retains the official eight-level full-profile oracle.
        Explicitly filtered runs may select four or more complete levels.  All
        controlled series must contain the same planned levels; otherwise the
        frozen plan is not a comparable convergence experiment.
        """

        if not planned_cases:
            raise AnalysisError("cannot configure analyzer from an empty plan")
        planned_groups: dict[str, set[int]] = {}
        for case in planned_cases:
            if case.spec.family != self.family:
                raise AnalysisError(
                    f"analyzer {self.name} cannot analyze family {case.spec.family}"
                )
            configuration = case.spec.to_dict()
            configuration.pop("time", None)
            planned_groups.setdefault(_canonical(configuration), set()).add(
                case.spec.time.steps
            )
        level_sets = {
            tuple(sorted(levels)) for levels in planned_groups.values()
        }
        if len(level_sets) != 1:
            raise AnalysisError(
                "frozen plan has incomplete or inconsistent temporal series"
            )
        expected_steps = next(iter(level_sets))
        if len(expected_steps) < REGRESSION_POINT_COUNT:
            raise AnalysisError(
                "temporal convergence requires at least four planned levels"
            )
        return replace(
            self,
            expected_steps=expected_steps,
            expected_groups=tuple(sorted(planned_groups)),
        )

    def analyze(
        self,
        results: Iterable[Mapping[str, object]],
        *,
        source_run: str | None = None,
    ) -> AnalysisReport:
        documents = tuple(results)
        if not documents:
            raise AnalysisError("temporal analysis received no results")

        seen_ids: set[str] = set()
        groups: dict[str, list[_Point]] = {}
        group_configurations: dict[str, Mapping[str, object]] = {}
        for index, document in enumerate(documents):
            if not isinstance(document, Mapping):
                raise AnalysisError(f"result {index} is not an object")
            point = _extract_point(document)
            if point.case_id in seen_ids:
                raise AnalysisError(f"duplicate case_id: {point.case_id}")
            seen_ids.add(point.case_id)
            group_key = _canonical(point.configuration_without_time)
            groups.setdefault(group_key, []).append(point)
            group_configurations[group_key] = point.configuration_without_time

        if self.expected_groups is not None:
            expected_groups = set(self.expected_groups)
            actual_groups = set(groups)
            if actual_groups != expected_groups:
                missing = sorted(expected_groups - actual_groups)
                unexpected = sorted(actual_groups - expected_groups)

                def labels(values: Sequence[str]) -> list[str]:
                    return [
                        hashlib.sha256(value.encode("utf-8")).hexdigest()[:16]
                        for value in values
                    ]

                raise AnalysisError(
                    "results do not match frozen planned groups; "
                    f"missing={labels(missing)}, unexpected={labels(unexpected)}"
                )

        series_reports: list[dict[str, object]] = []
        diagnostics: list[str] = []
        failed_series = 0
        for group_key in sorted(groups):
            points = sorted(groups[group_key], key=lambda item: item.steps)
            label = f"{points[0].problem}/{points[0].scheme}"
            actual_steps = tuple(point.steps for point in points)
            if actual_steps != self.expected_steps:
                missing = sorted(set(self.expected_steps) - set(actual_steps))
                unexpected = sorted(set(actual_steps) - set(self.expected_steps))
                raise AnalysisError(
                    f"{label}: incomplete temporal levels; "
                    f"expected={list(self.expected_steps)}, actual={list(actual_steps)}, "
                    f"missing={missing}, unexpected={unexpected}"
                )
            if len({point.case_id for point in points}) != len(points):
                raise AnalysisError(f"{label}: duplicate case_id within series")
            if len({point.final_time for point in points}) != 1:
                raise AnalysisError(f"{label}: inconsistent final time T")
            if len({point.time_step for point in points}) != len(points):
                raise AnalysisError(f"{label}: duplicate temporal resolution")
            if any(
                points[index].time_step <= points[index + 1].time_step
                for index in range(len(points) - 1)
            ):
                raise AnalysisError(f"{label}: dt is not strictly decreasing")

            expectation = SCHEME_EXPECTATIONS[points[0].scheme]
            l2, l2_diagnostics = _analyze_metric(points, "l2", expectation)
            linf, linf_diagnostics = _analyze_metric(points, "linf", expectation)
            status = (
                "passed"
                if l2["status"] == "passed" and linf["status"] == "passed"
                else "failed"
            )
            if status == "failed":
                failed_series += 1
            group_id = hashlib.sha256(group_key.encode("utf-8")).hexdigest()[:16]
            for message in (*l2_diagnostics, *linf_diagnostics):
                diagnostics.append(f"{group_id}: {message}")
            series_reports.append(
                {
                    "group_id": group_id,
                    "family": "temporal",
                    "problem": points[0].problem,
                    "scheme": points[0].scheme,
                    "final_time": float(points[0].final_time),
                    "controlled_configuration": dict(group_configurations[group_key]),
                    "status": status,
                    "points": [
                        {
                            "case_id": point.case_id,
                            "steps": point.steps,
                            "dt": float(point.time_step),
                            "l2_error": point.l2_error,
                            "linf_error": point.linf_error,
                        }
                        for point in points
                    ],
                    "metrics": {"l2": l2, "linf": linf},
                }
            )

        series_reports.sort(
            key=lambda item: (item["problem"], item["scheme"], item["group_id"])
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
                "levels_per_series": len(self.expected_steps),
                "local_orders_per_metric": len(self.expected_steps) - 1,
                "expected_steps": list(self.expected_steps),
            },
            series=tuple(series_reports),
            diagnostics=tuple(diagnostics),
        )
