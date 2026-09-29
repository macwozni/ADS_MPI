"""Reusable benchmark analysis pipeline and built-in analyzers."""

from .base import AnalysisError, AnalysisPipeline, AnalysisReport, Analyzer
from .io import WrittenAnalysis, load_run_results, write_report, write_run_report
from .spatial import HConvergenceAnalyzer, PConvergenceAnalyzer
from .temporal import TemporalConvergenceAnalyzer

__all__ = [
    "AnalysisError",
    "AnalysisPipeline",
    "AnalysisReport",
    "Analyzer",
    "HConvergenceAnalyzer",
    "PConvergenceAnalyzer",
    "TemporalConvergenceAnalyzer",
    "WrittenAnalysis",
    "load_run_results",
    "write_report",
    "write_run_report",
]
