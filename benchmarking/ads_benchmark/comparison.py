"""Numerics-first A/B comparison for already verified frozen benchmark runs.

This module deliberately does not open result directories.  A caller must use
the normal result store and completion verifier to construct :class:`RunSnapshot`
instances.  Keeping filesystem ownership and historical-command verification
outside the statistical layer avoids a second, weaker result loader.

``execution_provenance`` is a normalized *compatibility identity*: it should
contain the machine, toolchain, linked-library, launcher, binding, and relevant
environment facts that must be equal for a meaningful performance comparison.
Source commit, executable digest, build-root path, run ID, and timestamps are
intentionally represented elsewhere and must not be included in that identity.
Legacy runs may pass ``None``; their numerical and descriptive timing data can
still be inspected, but they can never receive a green or regression verdict.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
import json
import math
from typing import Any

from .analysis.base import AnalysisError
from .analysis.statistics import SampleStatistics, summarize_samples
from .analysis.validation import extract_field_result
from .framework.model import PlannedCase
from .framework.planner import FrozenPlan
from .validation import FieldComparisonError, compare_fields


COMPARISON_SCHEMA_VERSION = 1
COMPARISON_KIND = "ads-benchmark-comparison"
COMPARISON_STATUSES = frozenset(
    {
        "no-regression-detected",
        "regression",
        "numerical-mismatch",
        "incompatible",
        "inconclusive",
    }
)


class ComparisonError(AnalysisError):
    """An A/B policy or caller-provided snapshot is malformed."""


def _canonical(value: object, field: str) -> str:
    try:
        return json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
    except (TypeError, ValueError) as error:
        raise ComparisonError(f"{field} is not strict JSON: {error}") from error


def _finite_nonnegative(value: object, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ComparisonError(f"{field} must be a finite nonnegative number")
    result = float(value)
    if not math.isfinite(result) or result < 0.0:
        raise ComparisonError(f"{field} must be a finite nonnegative number")
    return result


@dataclass(frozen=True)
class ComparisonPolicy:
    """Frozen decision policy for one candidate-over-baseline comparison."""

    regression_threshold: float = 0.05
    minimum_samples: int = 5
    absolute_field_tolerance: float = 1.0e-11
    relative_field_tolerance: float = 1.0e-10

    def __post_init__(self) -> None:
        _finite_nonnegative(self.regression_threshold, "regression_threshold")
        if type(self.minimum_samples) is not int or self.minimum_samples < 3:
            raise ComparisonError("minimum_samples must be an integer >= 3")
        _finite_nonnegative(
            self.absolute_field_tolerance, "absolute_field_tolerance"
        )
        _finite_nonnegative(
            self.relative_field_tolerance, "relative_field_tolerance"
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "ratio_definition": "candidate_median/baseline_median",
            "regression_rule": "median_ratio > 1 + regression_threshold",
            "regression_threshold": float(self.regression_threshold),
            "minimum_samples": self.minimum_samples,
            "absolute_field_tolerance": float(
                self.absolute_field_tolerance
            ),
            "relative_field_tolerance": float(
                self.relative_field_tolerance
            ),
        }


@dataclass(frozen=True)
class RunSnapshot:
    """One safely loaded frozen plan and its fully verified result documents.

    Use :meth:`from_loaded` with results returned by the shared run loader.  The
    comparator validates the exact case/configuration correspondence again,
    but it does not replace status/log/command verification at the storage
    boundary.
    """

    plan: FrozenPlan
    results: tuple[Mapping[str, object], ...]
    execution_provenance: Mapping[str, object] | None

    @classmethod
    def from_loaded(
        cls,
        plan: FrozenPlan,
        results: Iterable[Mapping[str, object]],
        *,
        execution_provenance: Mapping[str, object] | None,
    ) -> "RunSnapshot":
        if not isinstance(plan, FrozenPlan):
            raise ComparisonError("snapshot plan must be a FrozenPlan")
        if isinstance(results, Mapping):
            raise ComparisonError("snapshot results must be an iterable of documents")
        documents = tuple(results)
        if any(not isinstance(document, Mapping) for document in documents):
            raise ComparisonError("every snapshot result must be an object")
        if execution_provenance is not None:
            if not isinstance(execution_provenance, Mapping):
                raise ComparisonError("execution_provenance must be an object")
            _canonical(execution_provenance, "execution_provenance")
        return cls(
            plan=plan,
            results=documents,
            execution_provenance=(
                dict(execution_provenance)
                if execution_provenance is not None
                else None
            ),
        )


@dataclass(frozen=True)
class ComparisonReport:
    """Versioned, strict-JSON-compatible A/B report."""

    status: str
    baseline: Mapping[str, object]
    candidate: Mapping[str, object]
    compatibility: Mapping[str, object]
    policy: Mapping[str, object]
    summary: Mapping[str, object]
    cases: tuple[Mapping[str, object], ...]
    diagnostics: tuple[str, ...]

    def __post_init__(self) -> None:
        if self.status not in COMPARISON_STATUSES:
            raise ComparisonError(f"unsupported comparison status: {self.status}")

    def to_dict(self) -> dict[str, object]:
        document = {
            "schema_version": COMPARISON_SCHEMA_VERSION,
            "kind": COMPARISON_KIND,
            "status": self.status,
            "baseline": dict(self.baseline),
            "candidate": dict(self.candidate),
            "compatibility": dict(self.compatibility),
            "policy": dict(self.policy),
            "summary": dict(self.summary),
            "cases": [dict(case) for case in self.cases],
            "diagnostics": list(self.diagnostics),
        }
        _canonical(document, "comparison report")
        return document


def _hash_label(config_hash: str) -> str:
    return config_hash if config_hash.startswith("sha256:") else f"sha256:{config_hash}"


def _snapshot_document(snapshot: RunSnapshot) -> dict[str, object]:
    return {
        "run_id": snapshot.plan.run_id,
        "profile": snapshot.plan.profile_name,
        "config_hash": _hash_label(snapshot.plan.config_hash),
        "repository": snapshot.plan.repository.to_dict(),
        "case_count": len(snapshot.plan.cases),
        "execution_provenance_present": snapshot.execution_provenance is not None,
    }


def _empty_summary(case_count: int) -> dict[str, object]:
    return {
        "case_count": case_count,
        "numerically_checked_cases": 0,
        "numerical_failure_count": 0,
        "timing_case_count": 0,
        "compared_timing_count": 0,
        "regression_count": 0,
        "no_regression_count": 0,
        "inconclusive_timing_count": 0,
        "timing_not_applicable_count": 0,
    }


def _early_report(
    baseline: RunSnapshot,
    candidate: RunSnapshot,
    policy: ComparisonPolicy,
    *,
    status: str,
    compatibility: Mapping[str, object],
    diagnostics: Sequence[str],
) -> ComparisonReport:
    return ComparisonReport(
        status=status,
        baseline=_snapshot_document(baseline),
        candidate=_snapshot_document(candidate),
        compatibility=dict(compatibility),
        policy=policy.to_dict(),
        summary=_empty_summary(len(baseline.plan.cases)),
        cases=(),
        diagnostics=tuple(diagnostics),
    )


def _case_map(plan: FrozenPlan, label: str) -> tuple[dict[str, PlannedCase], list[str]]:
    cases: dict[str, PlannedCase] = {}
    diagnostics: list[str] = []
    for case in plan.cases:
        if case.case_id in cases:
            diagnostics.append(f"{label} frozen plan duplicates case {case.case_id}")
        else:
            cases[case.case_id] = case
    return cases, diagnostics


def _result_map(
    snapshot: RunSnapshot,
    expected: Mapping[str, PlannedCase],
    label: str,
) -> tuple[dict[str, Mapping[str, object]], list[str]]:
    documents: dict[str, Mapping[str, object]] = {}
    diagnostics: list[str] = []
    for position, document in enumerate(snapshot.results):
        case_id = document.get("case_id")
        if not isinstance(case_id, str) or not case_id:
            diagnostics.append(f"{label} result {position} has an invalid case_id")
            continue
        if case_id in documents:
            diagnostics.append(f"{label} results duplicate case {case_id}")
            continue
        documents[case_id] = document

    expected_ids = set(expected)
    actual_ids = set(documents)
    missing = sorted(expected_ids - actual_ids)
    extra = sorted(actual_ids - expected_ids)
    if missing:
        diagnostics.append(f"{label} results miss cases: {', '.join(missing)}")
    if extra:
        diagnostics.append(f"{label} results contain extra cases: {', '.join(extra)}")
    for case_id in sorted(expected_ids & actual_ids):
        try:
            actual = _canonical(
                documents[case_id].get("configuration"),
                f"{label} result {case_id} configuration",
            )
            planned = _canonical(
                expected[case_id].spec.to_dict(),
                f"{label} planned case {case_id}",
            )
        except ComparisonError as error:
            diagnostics.append(str(error))
            continue
        if actual != planned:
            diagnostics.append(
                f"{label} result {case_id} configuration differs from its frozen plan"
            )
    return documents, diagnostics


def _compatibility(
    baseline: RunSnapshot, candidate: RunSnapshot
) -> tuple[dict[str, object], list[str], bool]:
    """Return compatibility record, fatal diagnostics, and provenance completeness."""

    diagnostics: list[str] = []
    baseline_cases, baseline_case_errors = _case_map(
        baseline.plan, "baseline"
    )
    candidate_cases, candidate_case_errors = _case_map(
        candidate.plan, "candidate"
    )
    diagnostics.extend(baseline_case_errors)
    diagnostics.extend(candidate_case_errors)

    if baseline.plan.config_hash != candidate.plan.config_hash:
        diagnostics.append("baseline and candidate config_hash values differ")
    if set(baseline_cases) != set(candidate_cases):
        missing = sorted(set(baseline_cases) - set(candidate_cases))
        extra = sorted(set(candidate_cases) - set(baseline_cases))
        diagnostics.append(
            "baseline and candidate case_id sets differ: "
            f"missing_in_candidate={missing}, extra_in_candidate={extra}"
        )
    for case_id in sorted(set(baseline_cases) & set(candidate_cases)):
        try:
            left = _canonical(
                baseline_cases[case_id].spec.to_dict(),
                f"baseline case {case_id}",
            )
            right = _canonical(
                candidate_cases[case_id].spec.to_dict(),
                f"candidate case {case_id}",
            )
        except ComparisonError as error:
            diagnostics.append(str(error))
            continue
        if left != right:
            diagnostics.append(
                f"case {case_id} has different baseline and candidate configurations"
            )

    _, baseline_result_errors = _result_map(
        baseline, baseline_cases, "baseline"
    )
    _, candidate_result_errors = _result_map(
        candidate, candidate_cases, "candidate"
    )
    diagnostics.extend(baseline_result_errors)
    diagnostics.extend(candidate_result_errors)
    configuration_match = not diagnostics

    baseline_provenance = baseline.execution_provenance
    candidate_provenance = candidate.execution_provenance
    provenance_complete = (
        baseline_provenance is not None and candidate_provenance is not None
    )
    if provenance_complete:
        assert baseline_provenance is not None
        assert candidate_provenance is not None
        if _canonical(
            baseline_provenance, "baseline execution_provenance"
        ) != _canonical(candidate_provenance, "candidate execution_provenance"):
            diagnostics.append(
                "baseline and candidate execution provenance are incompatible"
            )

    compatibility_reasons = list(diagnostics)
    if not provenance_complete:
        compatibility_reasons.append("execution provenance is missing")
    compatibility = {
        "passed": not diagnostics and provenance_complete,
        "configuration_match": configuration_match,
        "execution_provenance": (
            "compatible"
            if provenance_complete and not any(
                "execution provenance are incompatible" in item
                for item in diagnostics
            )
            else "incompatible"
            if provenance_complete
            else "missing"
        ),
        "reasons": compatibility_reasons,
    }
    return compatibility, diagnostics, provenance_complete


def _compact_field_comparison(comparison: Any) -> dict[str, object]:
    document = comparison.to_dict()
    return {
        "passed": document["passed"],
        "point_count": document["point_count"],
        "mismatch_count": document["mismatch_count"],
        "l2_difference": document["l2_difference"],
        "linf_difference": document["linf_difference"],
        "absolute_tolerance": document["absolute_tolerance"],
        "relative_tolerance": document["relative_tolerance"],
        "baseline_checksum": document["reference"]["checksum"],
        "candidate_checksum": document["candidate"]["checksum"],
        "worst_point": document["worst_point"],
        "diagnostic": document["diagnostic"],
    }


def _statistics_document(statistics: SampleStatistics) -> dict[str, object]:
    return statistics.to_dict() | {
        "range_seconds": statistics.maximum - statistics.minimum,
    }


@dataclass(frozen=True)
class _TimingInput:
    samples: tuple[float, ...]
    reliable: bool


def _timing_input(
    document: Mapping[str, object], case: PlannedCase, label: str
) -> _TimingInput:
    timing = document.get("timing")
    if not isinstance(timing, Mapping):
        raise ComparisonError(f"{label} {case.case_id} timing must be an object")
    if timing.get("metric") != "physical_step_wall_seconds":
        raise ComparisonError(
            f"{label} {case.case_id} timing metric must be physical_step_wall_seconds"
        )
    raw_samples = timing.get("measured_samples")
    if not isinstance(raw_samples, list):
        raise ComparisonError(
            f"{label} {case.case_id} measured_samples must be an array"
        )
    if len(raw_samples) != case.spec.measurement.samples:
        raise ComparisonError(
            f"{label} {case.case_id} measured sample count differs from the frozen plan"
        )
    samples: list[float] = []
    for index, raw_value in enumerate(raw_samples):
        if isinstance(raw_value, bool) or not isinstance(raw_value, (int, float)):
            raise ComparisonError(
                f"{label} {case.case_id} measured_samples[{index}] must be a number"
            )
        value = float(raw_value)
        if not math.isfinite(value) or value <= 0.0:
            raise ComparisonError(
                f"{label} {case.case_id} measured_samples[{index}] "
                "must be finite and positive"
            )
        samples.append(value)
    reliable = timing.get("reliable")
    if type(reliable) is not bool:
        raise ComparisonError(
            f"{label} {case.case_id} timing.reliable must be a boolean"
        )
    return _TimingInput(samples=tuple(samples), reliable=reliable)


def _timing_applicable(case: PlannedCase) -> bool:
    weak = case.spec.weak_scaling
    return weak is None or weak.role == "measurement"


def compare_snapshots(
    baseline: RunSnapshot,
    candidate: RunSnapshot,
    policy: ComparisonPolicy | None = None,
) -> ComparisonReport:
    """Compare two complete frozen runs with a global numerics-first gate.

    The candidate/baseline median ratio is classified only after *every* case
    has passed full-field comparison.  A missing execution-provenance identity,
    too few samples, unreliable timing, or malformed timing makes the affected
    performance comparison inconclusive and can never create a regression.
    """

    if not isinstance(baseline, RunSnapshot) or not isinstance(candidate, RunSnapshot):
        raise ComparisonError("baseline and candidate must be RunSnapshot values")
    selected_policy = policy or ComparisonPolicy()
    if not isinstance(selected_policy, ComparisonPolicy):
        raise ComparisonError("policy must be a ComparisonPolicy")

    compatibility, compatibility_errors, provenance_complete = _compatibility(
        baseline, candidate
    )
    if compatibility_errors:
        return _early_report(
            baseline,
            candidate,
            selected_policy,
            status="incompatible",
            compatibility=compatibility,
            diagnostics=compatibility_errors,
        )

    baseline_cases = {case.case_id: case for case in baseline.plan.cases}
    baseline_results, _ = _result_map(
        baseline, baseline_cases, "baseline"
    )
    candidate_results, _ = _result_map(
        candidate, baseline_cases, "candidate"
    )

    case_reports: list[dict[str, object]] = []
    diagnostics: list[str] = []
    numerical_failures = 0
    for case_id in sorted(baseline_cases):
        case = baseline_cases[case_id]
        numerical: dict[str, object]
        try:
            baseline_point = extract_field_result(
                baseline_results[case_id], expected_family=case.spec.family
            )
            candidate_point = extract_field_result(
                candidate_results[case_id], expected_family=case.spec.family
            )
            comparison = compare_fields(
                baseline_point.numerical,
                candidate_point.numerical,
                absolute_tolerance=selected_policy.absolute_field_tolerance,
                relative_tolerance=selected_policy.relative_field_tolerance,
            )
            numerical = _compact_field_comparison(comparison)
            if not comparison.passed:
                numerical_failures += 1
                diagnostics.append(
                    f"{case_id}: candidate field differs from baseline; "
                    f"{comparison.diagnostic()}"
                )
        except (AnalysisError, FieldComparisonError, ComparisonError) as error:
            numerical_failures += 1
            numerical = {"passed": False, "error": str(error)}
            diagnostics.append(f"{case_id}: numerical gate failed: {error}")
        case_reports.append(
            {
                "case_id": case_id,
                "numerical": numerical,
            }
        )

    if numerical_failures:
        for case_report in case_reports:
            case_report["timing"] = {"status": "blocked-by-numerical-gate"}
        summary = _empty_summary(len(case_reports))
        summary.update(
            {
                "numerically_checked_cases": len(case_reports),
                "numerical_failure_count": numerical_failures,
            }
        )
        return ComparisonReport(
            status="numerical-mismatch",
            baseline=_snapshot_document(baseline),
            candidate=_snapshot_document(candidate),
            compatibility=compatibility,
            policy=selected_policy.to_dict(),
            summary=summary,
            cases=tuple(case_reports),
            diagnostics=tuple(diagnostics),
        )

    timing_case_count = 0
    compared_timing_count = 0
    regressions = 0
    no_regressions = 0
    inconclusive = 0
    not_applicable = 0

    for case_report in case_reports:
        case_id = str(case_report["case_id"])
        case = baseline_cases[case_id]
        if not _timing_applicable(case):
            case_report["timing"] = {"status": "not-applicable"}
            not_applicable += 1
            continue
        timing_case_count += 1
        try:
            baseline_timing = _timing_input(
                baseline_results[case_id], case, "baseline"
            )
            candidate_timing = _timing_input(
                candidate_results[case_id], case, "candidate"
            )
            baseline_statistics = summarize_samples(baseline_timing.samples)
            candidate_statistics = summarize_samples(candidate_timing.samples)
            ratio = candidate_statistics.median / baseline_statistics.median
            timing: dict[str, object] = {
                "metric": "physical_step_wall_seconds",
                "baseline": _statistics_document(baseline_statistics),
                "candidate": _statistics_document(candidate_statistics),
                "median_ratio_candidate_over_baseline": ratio,
            }
            reasons: list[str] = []
            if len(baseline_timing.samples) < selected_policy.minimum_samples:
                reasons.append(
                    "baseline has fewer than the configured minimum samples"
                )
            if len(candidate_timing.samples) < selected_policy.minimum_samples:
                reasons.append(
                    "candidate has fewer than the configured minimum samples"
                )
            if not baseline_timing.reliable:
                reasons.append("baseline timing is marked unreliable")
            if not candidate_timing.reliable:
                reasons.append("candidate timing is marked unreliable")
            if not provenance_complete:
                reasons.append("execution provenance is missing")

            if reasons:
                timing["status"] = "inconclusive"
                timing["reasons"] = reasons
                inconclusive += 1
                diagnostics.append(f"{case_id}: " + "; ".join(reasons))
            elif ratio > 1.0 + selected_policy.regression_threshold:
                timing["status"] = "regression"
                regressions += 1
                compared_timing_count += 1
            else:
                timing["status"] = "no-regression-detected"
                no_regressions += 1
                compared_timing_count += 1
            case_report["timing"] = timing
        except (AnalysisError, ComparisonError, ZeroDivisionError) as error:
            case_report["timing"] = {
                "status": "inconclusive",
                "reasons": [str(error)],
            }
            inconclusive += 1
            diagnostics.append(f"{case_id}: timing comparison failed: {error}")

    if regressions:
        status = "regression"
    elif inconclusive or timing_case_count == 0:
        status = "inconclusive"
        if timing_case_count == 0:
            diagnostics.append("comparison contains no timing-applicable cases")
    else:
        status = "no-regression-detected"

    summary = _empty_summary(len(case_reports))
    summary.update(
        {
            "numerically_checked_cases": len(case_reports),
            "numerical_failure_count": 0,
            "timing_case_count": timing_case_count,
            "compared_timing_count": compared_timing_count,
            "regression_count": regressions,
            "no_regression_count": no_regressions,
            "inconclusive_timing_count": inconclusive,
            "timing_not_applicable_count": not_applicable,
        }
    )
    return ComparisonReport(
        status=status,
        baseline=_snapshot_document(baseline),
        candidate=_snapshot_document(candidate),
        compatibility=compatibility,
        policy=selected_policy.to_dict(),
        summary=summary,
        cases=tuple(case_reports),
        diagnostics=tuple(diagnostics),
    )


__all__ = [
    "COMPARISON_KIND",
    "COMPARISON_SCHEMA_VERSION",
    "COMPARISON_STATUSES",
    "ComparisonError",
    "ComparisonPolicy",
    "ComparisonReport",
    "RunSnapshot",
    "compare_snapshots",
]
