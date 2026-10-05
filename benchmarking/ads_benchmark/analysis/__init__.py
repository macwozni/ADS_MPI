"""Reusable benchmark analysis pipeline and built-in analyzers."""

from .base import AnalysisError, AnalysisPipeline, AnalysisReport, Analyzer
from .io import WrittenAnalysis, load_run_results, write_report, write_run_report
from .spatial import HConvergenceAnalyzer, PConvergenceAnalyzer
from .statistics import (
    SampleStatistics,
    ScalingMetrics,
    scaling_metrics,
    summarize_samples,
)
from .strong import StrongScalingAnalyzer
from .temporal import TemporalConvergenceAnalyzer
from .validation import FieldValidationAnalyzer

__all__ = [
    "AnalysisError",
    "AnalysisPipeline",
    "AnalysisReport",
    "Analyzer",
    "HConvergenceAnalyzer",
    "PConvergenceAnalyzer",
    "SampleStatistics",
    "ScalingMetrics",
    "StrongScalingAnalyzer",
    "FieldValidationAnalyzer",
    "TemporalConvergenceAnalyzer",
    "WrittenAnalysis",
    "load_run_results",
    "write_report",
    "write_run_report",
    "scaling_metrics",
    "summarize_samples",
]
