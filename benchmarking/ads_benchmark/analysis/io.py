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
from ..framework.errors import StorageError
from ..framework.model import PlannedCase
from ..framework.storage import ResultStore
from .base import AnalysisError, AnalysisReport


@dataclass(frozen=True)
class WrittenAnalysis:
    json_path: Path
    csv_path: Path
    plot_path: Path | None
    plot_message: str | None


@dataclass(frozen=True)
class GeneratedArtifactReference:
    """Lazy, safety-checked reference to one generated case artifact."""

    store: ResultStore
    case_directory: Path
    case_id: str
    name: str

    def read_text(self) -> str:
        try:
            return self.store.read_generated_artifact(
                self.case_directory, self.name
            )
        except StorageError as error:
            raise AnalysisError(
                f"case {self.case_id} {self.name} cannot be loaded: {error}"
            ) from error


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
        document = dict(result)
        if case.spec.sampling.write_samples:
            document["analysis_artifacts"] = {
                "field_samples_csv": GeneratedArtifactReference(
                    store=executor.store,
                    case_directory=case_directory,
                    case_id=case.case_id,
                    name="field_samples.csv",
                ),
            }
        documents.append(document)
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


def _temporal_csv_text(report: AnalysisReport) -> str:
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


def _compact_json(value: object) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )


def _spatial_csv_text(report: AnalysisReport) -> str:
    """Render point qualification and derived estimates without losing vectors."""

    fields = (
        "group_id",
        "family",
        "sequence_kind",
        "acceptance_mode",
        "acceptance_reason",
        "problem",
        "scheme",
        "series_status",
        "metric",
        "metric_status",
        "row_kind",
        "coarse_case_id",
        "fine_case_id",
        "coarse_steps",
        "fine_steps",
        "coarse_dt",
        "fine_dt",
        "mesh",
        "h",
        "trial_degree",
        "test_degree",
        "enrichment",
        "coarse_error",
        "fine_error",
        "temporal_field_difference",
        "temporal_error_fraction",
        "temporal_indicator_status",
        "maximum_temporal_error_fraction",
        "reliable_for_spatial_analysis",
        "sample_points_per_axis",
        "difference_method",
        "coarse_mesh",
        "fine_mesh",
        "coarse_h",
        "fine_h",
        "coarse_trial_degree",
        "fine_trial_degree",
        "coarse_test_degree",
        "fine_test_degree",
        "local_order",
        "error_ratio",
        "error_decreased",
        "available",
        "unavailable_reason",
        "in_plateau",
        "fatal_nonmonotonic",
        "used_in_regression",
        "regression_order",
        "regression_r_squared",
        "regression_point_count",
        "plateau_detected",
        "rotation_error_spread",
    )
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()

    vector_fields = {
        "mesh",
        "trial_degree",
        "test_degree",
        "enrichment",
        "coarse_mesh",
        "fine_mesh",
        "coarse_trial_degree",
        "fine_trial_degree",
        "coarse_test_degree",
        "fine_test_degree",
    }

    def write(values: Mapping[str, object]) -> None:
        normalized = dict(values)
        for field in vector_fields:
            value = normalized.get(field)
            if value is not None:
                normalized[field] = _compact_json(value)
        writer.writerow(normalized)

    for series in report.series:
        common = {
            "group_id": series["group_id"],
            "family": series["family"],
            "sequence_kind": series["sequence_kind"],
            "acceptance_mode": series["acceptance_mode"],
            "problem": series["problem"],
            "scheme": series["scheme"],
            "series_status": series["status"],
        }
        for metric_name in ("l2", "linf"):
            metric = series["metrics"][metric_name]  # type: ignore[index]
            metric_common = common | {
                "metric": metric_name,
                "metric_status": metric["status"],
                "acceptance_reason": metric.get("acceptance_reason"),
                "plateau_detected": metric["plateau_detected"],
                "rotation_error_spread": metric["rotation_error_spread"],
            }
            for point in series["points"]:  # type: ignore[assignment]
                audit = point["temporal_pair"][metric_name]
                write(
                    metric_common
                    | {
                        "row_kind": "point",
                        "coarse_case_id": point["coarse_case_id"],
                        "fine_case_id": point["fine_case_id"],
                        "coarse_steps": point["coarse_steps"],
                        "fine_steps": point["fine_steps"],
                        "coarse_dt": point["coarse_dt"],
                        "fine_dt": point["fine_dt"],
                        "mesh": point.get("mesh"),
                        "h": point.get("h"),
                        "trial_degree": point.get("trial_degree"),
                        "test_degree": point.get("test_degree"),
                        "enrichment": point.get("enrichment"),
                        "coarse_error": audit["coarse_error"],
                        "fine_error": audit["fine_error"],
                        "temporal_field_difference": audit[
                            "temporal_field_difference"
                        ],
                        "temporal_error_fraction": audit[
                            "temporal_error_fraction"
                        ],
                        "temporal_indicator_status": audit[
                            "temporal_indicator_status"
                        ],
                        "maximum_temporal_error_fraction": audit[
                            "maximum_temporal_error_fraction"
                        ],
                        "reliable_for_spatial_analysis": audit[
                            "reliable_for_spatial_analysis"
                        ],
                        "sample_points_per_axis": audit[
                            "sample_points_per_axis"
                        ],
                        "difference_method": audit["difference_method"],
                    }
                )

            estimates = (
                metric["local_orders"]
                if report.family == "h"
                else metric["error_ratios"]
            )
            for estimate in estimates:
                estimate_values = dict(estimate)
                local_order = estimate_values.pop("order", None)
                write(
                    metric_common
                    | estimate_values
                    | {
                        "row_kind": (
                            "local-order" if report.family == "h" else "error-ratio"
                        ),
                        "local_order": local_order,
                    }
                )

            regression = metric["regression"]
            if regression is not None:
                write(
                    metric_common
                    | {
                        "row_kind": "regression",
                        "regression_order": regression["order"],
                        "regression_r_squared": regression["r_squared"],
                        "regression_point_count": regression["point_count"],
                    }
                )
            elif series["sequence_kind"] == "anisotropic-rotations":
                write(metric_common | {"row_kind": "rotation-summary"})
    return stream.getvalue()


def _csv_text(report: AnalysisReport) -> str:
    if report.family == "temporal":
        return _temporal_csv_text(report)
    if report.family in {"h", "p"}:
        return _spatial_csv_text(report)
    raise AnalysisError(f"no CSV renderer registered for family {report.family}")


def _plot_coordinates(
    report: AnalysisReport, series: Mapping[str, Any]
) -> tuple[list[float], str, str]:
    points = series["points"]
    if report.family == "temporal":
        return [float(point["dt"]) for point in points], "dt", "log"
    if report.family == "h":
        return [float(point["h"]) for point in points], "h", "log"
    if series["sequence_kind"] == "isotropic-degree":
        return (
            [float(point["trial_degree"][0]) for point in points],
            "trial degree",
            "linear",
        )
    return [float(index) for index in range(1, len(points) + 1)], "rotation", "linear"


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
                xs, xlabel, xscale = _plot_coordinates(report, series)
                axis.plot(
                    xs,
                    [point[f"{metric_name}_error"] for point in points],
                    marker="o",
                    linewidth=1,
                    label=label,
                )
                axis.set_xscale(xscale)
            axis.set_yscale("log")
            axis.set_xlabel(xlabel)
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
