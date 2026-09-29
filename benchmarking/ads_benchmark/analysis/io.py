"""Strict run loading and dependency-free analysis report writers."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import csv
import io
import json
import os
from pathlib import Path
import tempfile
from typing import Any

from ..framework.executor import Executor
from ..framework.model import PlannedCase
from ..framework.storage import ResultStore
from .base import AnalysisError, AnalysisReport


@dataclass(frozen=True)
class WrittenAnalysis:
    json_path: Path
    csv_path: Path
    plot_path: Path | None
    plot_message: str | None


def load_run_results(
    executor: Executor,
    run_id: str,
    planned_cases: Sequence[PlannedCase],
) -> tuple[dict[str, Any], ...]:
    """Load only cases accepted by the executor's shared completion oracle.

    ``planned_cases`` must come from the already validated frozen manifest.
    This deliberately delegates every status/result/log check, adapter reparse,
    domain validation, strict persisted-result comparison, and sample-artifact
    check to the same verifier used by resume.
    """

    if not planned_cases:
        raise AnalysisError("frozen plan contains no cases to analyze")
    seen: set[str] = set()
    documents: list[dict[str, Any]] = []
    for case in sorted(planned_cases, key=lambda item: item.case_id):
        if case.case_id in seen:
            raise AnalysisError(f"duplicate planned case_id: {case.case_id}")
        seen.add(case.case_id)
        case_directory = (
            executor.store.results_root / run_id / "cases" / case.case_id
        )
        result = executor.verified_result(case_directory, case)
        if result is None:
            raise AnalysisError(
                f"case {case.case_id} is not a fully verified completed result"
            )
        documents.append(dict(result))
    return tuple(documents)


def _atomic_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except OSError:
            pass
        raise


def _atomic_bytes(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except OSError:
            pass
        raise


def _json_text(report: AnalysisReport) -> str:
    try:
        return json.dumps(
            report.to_dict(),
            indent=2,
            sort_keys=True,
            ensure_ascii=True,
            allow_nan=False,
        ) + "\n"
    except (TypeError, ValueError) as error:
        raise AnalysisError(f"analysis report is not strict JSON: {error}") from error


def _csv_text(report: AnalysisReport) -> str:
    fields = (
        "group_id",
        "problem",
        "scheme",
        "metric",
        "status",
        "expected_order",
        "accepted_min",
        "accepted_max",
        "regression_order",
        "regression_first_steps",
        "regression_last_steps",
        "plateau_start_steps",
        "coarse_error",
        "finest_usable_error",
        "analytical_reduction",
        "coarse_steps",
        "fine_steps",
        "coarse_dt",
        "fine_dt",
        "local_order",
        "maximum_local_order",
        "implausible",
        "nonmonotonic",
        "fatal_nonmonotonic",
        "used_in_regression",
    )
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for series in report.series:
        for metric_name in ("l2", "linf"):
            metric = series["metrics"][metric_name]  # type: ignore[index]
            common = {
                "group_id": series["group_id"],
                "problem": series["problem"],
                "scheme": series["scheme"],
                "metric": metric_name,
                "status": metric["status"],
                "expected_order": metric["expected_order"],
                "accepted_min": metric["accepted_interval"][0],
                "accepted_max": metric["accepted_interval"][1],
                "regression_order": metric["regression"]["order"],
                "regression_first_steps": metric["regression"]["first_steps"],
                "regression_last_steps": metric["regression"]["last_steps"],
                "plateau_start_steps": metric["plateau_start_steps"],
                "coarse_error": metric["coarse_error"],
                "finest_usable_error": metric["finest_usable_error"],
                "analytical_reduction": metric["analytical_reduction"],
            }
            for estimate in metric["local_orders"]:
                writer.writerow(
                    common
                    | {
                        "coarse_steps": estimate["coarse_steps"],
                        "fine_steps": estimate["fine_steps"],
                        "coarse_dt": estimate["coarse_dt"],
                        "fine_dt": estimate["fine_dt"],
                        "local_order": estimate["order"],
                        "maximum_local_order": metric["maximum_local_order"],
                        "implausible": estimate["implausible"],
                        "nonmonotonic": estimate["nonmonotonic"],
                        "fatal_nonmonotonic": estimate["fatal_nonmonotonic"],
                        "used_in_regression": estimate["used_in_regression"],
                    }
                )
    return stream.getvalue()


def _render_plot(report: AnalysisReport) -> tuple[bytes | None, str | None]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        from matplotlib import pyplot as plt
    except (ImportError, ModuleNotFoundError) as error:
        return None, f"plot omitted: matplotlib is unavailable ({error})"

    figure, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    try:
        for axis, metric_name, title in zip(
            axes, ("l2", "linf"), ("L2 error", "Linf error"), strict=True
        ):
            for series in report.series:
                points = series["points"]
                label = f"{series['problem']}/{series['scheme']}/{series['group_id']}"
                axis.loglog(
                    [point["dt"] for point in points],
                    [point[f"{metric_name}_error"] for point in points],
                    marker="o",
                    linewidth=1,
                    label=label,
                )
            axis.set_xlabel("dt")
            axis.set_ylabel(title)
            axis.grid(True, which="both", alpha=0.25)
        if len(report.series) <= 12:
            axes[1].legend(fontsize="x-small")
        buffer = io.BytesIO()
        figure.savefig(buffer, format="png", dpi=150)
    finally:
        plt.close(figure)
    return buffer.getvalue(), None


def write_report(
    report: AnalysisReport,
    output_directory: Path,
    *,
    plot: bool = False,
) -> WrittenAnalysis:
    """Atomically write JSON/CSV, optionally omitting only an unavailable plot."""

    output_directory = output_directory.resolve()
    json_path = output_directory / "analysis.json"
    csv_path = output_directory / "analysis.csv"
    stale_plot = output_directory / "convergence.png"
    try:
        stale_plot.unlink()
    except FileNotFoundError:
        pass
    _atomic_text(json_path, _json_text(report))
    _atomic_text(csv_path, _csv_text(report))

    plot_path: Path | None = None
    plot_message: str | None = None
    if plot:
        plot_bytes, plot_message = _render_plot(report)
        if plot_bytes is not None:
            plot_path = output_directory / "convergence.png"
            _atomic_bytes(plot_path, plot_bytes)
    return WrittenAnalysis(json_path, csv_path, plot_path, plot_message)


def write_run_report(
    report: AnalysisReport,
    store: ResultStore,
    run_id: str,
    *,
    plot: bool = False,
) -> WrittenAnalysis:
    """Write only below ``benchmarks/<run-id>/analysis`` via safe dirfds."""

    # A plot is derived from this exact report.  Remove a previous one first,
    # including when plotting was not requested or matplotlib is unavailable;
    # ResultStore unlinks the directory entry without following a symlink.
    store.remove_analysis_artifact(run_id, "convergence.png")
    json_path = store.write_analysis_artifact(
        run_id, "analysis.json", _json_text(report)
    )
    csv_path = store.write_analysis_artifact(
        run_id, "analysis.csv", _csv_text(report)
    )
    plot_path: Path | None = None
    plot_message: str | None = None
    if plot:
        plot_bytes, plot_message = _render_plot(report)
        if plot_bytes is not None:
            plot_path = store.write_analysis_artifact(
                run_id, "convergence.png", plot_bytes
            )
    return WrittenAnalysis(json_path, csv_path, plot_path, plot_message)
