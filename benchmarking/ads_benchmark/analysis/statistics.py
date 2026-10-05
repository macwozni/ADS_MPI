"""Small deterministic statistics used by scaling analyses."""

from __future__ import annotations

from dataclasses import dataclass
import math
from statistics import median
from typing import Iterable

from .base import AnalysisError


@dataclass(frozen=True)
class SampleStatistics:
    count: int
    minimum: float
    maximum: float
    median: float
    median_absolute_deviation: float
    relative_median_absolute_deviation: float | None

    def to_dict(self) -> dict[str, object]:
        return {
            "count": self.count,
            "minimum_seconds": self.minimum,
            "maximum_seconds": self.maximum,
            "median_seconds": self.median,
            "median_absolute_deviation_seconds": self.median_absolute_deviation,
            "relative_median_absolute_deviation": (
                self.relative_median_absolute_deviation
            ),
        }


@dataclass(frozen=True)
class ScalingMetrics:
    speedup: float
    efficiency: float


def weak_scaling_efficiency(
    *, baseline_seconds: float, measured_seconds: float
) -> float:
    """Return weak-scaling efficiency without assigning speedup semantics.

    A weak-scaling series keeps its declared local-work policy constant, so
    ideal elapsed time is constant as resources and global work grow.  Its
    efficiency is consequently the baseline time divided by the measured
    time; no resource-count factor belongs in this metric.
    """

    if any(
        isinstance(value, bool) or not isinstance(value, (int, float))
        for value in (baseline_seconds, measured_seconds)
    ):
        raise AnalysisError("weak-scaling times must be finite and positive")
    if (
        not math.isfinite(baseline_seconds)
        or not math.isfinite(measured_seconds)
        or baseline_seconds <= 0.0
        or measured_seconds <= 0.0
    ):
        raise AnalysisError("weak-scaling times must be finite and positive")
    return baseline_seconds / measured_seconds


def summarize_samples(values: Iterable[float]) -> SampleStatistics:
    """Summarize finite, nonnegative raw seconds without discarding samples."""

    raw_samples = tuple(values)
    if not raw_samples:
        raise AnalysisError("timing sample set must not be empty")
    if any(
        isinstance(value, bool) or not isinstance(value, (int, float))
        for value in raw_samples
    ):
        raise AnalysisError("timing samples must be numbers")
    samples = tuple(float(value) for value in raw_samples)
    if any(not math.isfinite(value) or value < 0.0 for value in samples):
        raise AnalysisError("timing samples must be finite and nonnegative")
    center = float(median(samples))
    absolute_deviations = tuple(abs(value - center) for value in samples)
    deviation = float(median(absolute_deviations))
    relative_deviation = deviation / center if center > 0.0 else None
    return SampleStatistics(
        count=len(samples),
        minimum=min(samples),
        maximum=max(samples),
        median=center,
        median_absolute_deviation=deviation,
        relative_median_absolute_deviation=relative_deviation,
    )


def scaling_metrics(
    *,
    baseline_seconds: float,
    measured_seconds: float,
    baseline_resources: int,
    resources: int,
) -> ScalingMetrics:
    """Compute strong-scaling speedup and efficiency for one valid point."""

    if any(
        isinstance(value, bool) or not isinstance(value, (int, float))
        for value in (baseline_seconds, measured_seconds)
    ):
        raise AnalysisError("scaling times must be finite and positive")
    if (
        not math.isfinite(baseline_seconds)
        or not math.isfinite(measured_seconds)
        or baseline_seconds <= 0.0
        or measured_seconds <= 0.0
    ):
        raise AnalysisError("scaling times must be finite and positive")
    if (
        type(baseline_resources) is not int
        or type(resources) is not int
        or baseline_resources <= 0
        or resources < baseline_resources
    ):
        raise AnalysisError(
            "scaling resources must be positive and not below the baseline"
        )
    speedup = baseline_seconds / measured_seconds
    resource_ratio = resources / baseline_resources
    return ScalingMetrics(speedup=speedup, efficiency=speedup / resource_ratio)


__all__ = [
    "SampleStatistics",
    "ScalingMetrics",
    "scaling_metrics",
    "summarize_samples",
    "weak_scaling_efficiency",
]
