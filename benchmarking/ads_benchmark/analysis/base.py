"""Problem-neutral analysis extension contract and result schema."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Protocol

from ..framework.errors import BenchmarkError
from ..framework.model import PlannedCase
from ..framework.registry import Registry


ANALYSIS_SCHEMA_VERSION = 1


class AnalysisError(BenchmarkError):
    """Input data cannot support a trustworthy benchmark analysis."""


@dataclass(frozen=True)
class AnalysisReport:
    """Common, JSON-serializable result returned by every analyzer."""

    analyzer: str
    family: str
    status: str
    source_run: str | None
    summary: Mapping[str, Any]
    series: tuple[Mapping[str, Any], ...]
    diagnostics: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": ANALYSIS_SCHEMA_VERSION,
            "kind": "ads-benchmark-analysis",
            "analyzer": self.analyzer,
            "family": self.family,
            "status": self.status,
            "source_run": self.source_run,
            "summary": dict(self.summary),
            "series": [dict(item) for item in self.series],
            "diagnostics": list(self.diagnostics),
        }


class Analyzer(Protocol):
    """One registered analysis strategy over common case-result documents."""

    name: str
    family: str

    def configure(self, planned_cases: Sequence[PlannedCase]) -> "Analyzer":
        """Bind analysis expectations to the validated frozen plan."""

    def analyze(
        self,
        results: Iterable[Mapping[str, object]],
        *,
        source_run: str | None = None,
    ) -> AnalysisReport:
        """Validate and analyze result documents, or raise ``AnalysisError``."""


class AnalysisPipeline:
    """Registry-backed dispatcher independent of experiment-family details."""

    def __init__(self, analyzers: Registry[Analyzer] | None = None) -> None:
        self.analyzers = analyzers if analyzers is not None else Registry("analyzer")

    def register(self, analyzer: Analyzer) -> None:
        self.analyzers.register(analyzer.name, analyzer)

    def analyze(
        self,
        analyzer_name: str,
        results: Iterable[Mapping[str, object]],
        *,
        planned_cases: Sequence[PlannedCase] | None = None,
        source_run: str | None = None,
    ) -> AnalysisReport:
        analyzer = self.analyzers.get(analyzer_name)
        if planned_cases is not None:
            analyzer = analyzer.configure(planned_cases)
        return analyzer.analyze(results, source_run=source_run)
