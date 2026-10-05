"""Deterministic JSON/CSV serialization for A/B comparison reports."""

from __future__ import annotations

import csv
import io
import json
import os
from pathlib import Path
import time

from .comparison import ComparisonReport
from .framework.errors import StorageError


def comparison_json_text(report: ComparisonReport) -> str:
    return json.dumps(
        report.to_dict(),
        indent=2,
        sort_keys=True,
        ensure_ascii=True,
        allow_nan=False,
    ) + "\n"


def comparison_csv_text(report: ComparisonReport) -> str:
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(
        (
            "case_id",
            "numerical_passed",
            "timing_status",
            "baseline_count",
            "baseline_median_seconds",
            "baseline_mad_seconds",
            "baseline_range_seconds",
            "candidate_count",
            "candidate_median_seconds",
            "candidate_mad_seconds",
            "candidate_range_seconds",
            "candidate_over_baseline_median_ratio",
        )
    )
    for case in report.cases:
        numerical = case.get("numerical")
        timing = case.get("timing")
        numerical_document = numerical if isinstance(numerical, dict) else {}
        timing_document = timing if isinstance(timing, dict) else {}
        baseline = timing_document.get("baseline")
        candidate = timing_document.get("candidate")
        baseline_document = baseline if isinstance(baseline, dict) else {}
        candidate_document = candidate if isinstance(candidate, dict) else {}
        writer.writerow(
            (
                case.get("case_id", ""),
                numerical_document.get("passed", ""),
                timing_document.get("status", ""),
                baseline_document.get("count", ""),
                baseline_document.get("median_seconds", ""),
                baseline_document.get(
                    "median_absolute_deviation_seconds", ""
                ),
                baseline_document.get("range_seconds", ""),
                candidate_document.get("count", ""),
                candidate_document.get("median_seconds", ""),
                candidate_document.get(
                    "median_absolute_deviation_seconds", ""
                ),
                candidate_document.get("range_seconds", ""),
                timing_document.get(
                    "median_ratio_candidate_over_baseline", ""
                ),
            )
        )
    return output.getvalue()


def write_comparison_file(path: Path, content: str) -> Path:
    """Atomically replace one explicitly requested regular output file."""

    destination = path.absolute()
    parent = destination.parent
    if parent.is_symlink() or not parent.is_dir():
        raise StorageError(
            f"comparison output parent is missing or unsafe: {parent}"
        )
    if destination.is_symlink() or (destination.exists() and not destination.is_file()):
        raise StorageError(
            f"comparison output is not a regular file: {destination}"
        )
    temporary = parent / (
        f".{destination.name}.{os.getpid()}.{time.time_ns()}.tmp"
    )
    descriptor: int | None = None
    try:
        descriptor = os.open(
            temporary,
            os.O_WRONLY
            | os.O_CREAT
            | os.O_EXCL
            | getattr(os, "O_NOFOLLOW", 0),
            0o644,
        )
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            descriptor = None
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, destination)
        parent_fd = os.open(
            parent,
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0),
        )
        try:
            os.fsync(parent_fd)
        finally:
            os.close(parent_fd)
    except (OSError, UnicodeError) as error:
        if descriptor is not None:
            os.close(descriptor)
        try:
            temporary.unlink()
        except OSError:
            pass
        raise StorageError(
            f"cannot write comparison output {destination}: {error}"
        ) from error
    return destination


__all__ = [
    "comparison_csv_text",
    "comparison_json_text",
    "write_comparison_file",
]
