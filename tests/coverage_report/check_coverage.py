#!/usr/bin/env python3
"""Validate LCOV coverage for every active ADS core source.

LCOV itself can gate line coverage, but it does not provide matching gates for
functions and branches.  This checker applies all three thresholds and also
guards the denominator: every source declared in ``src/sources.mk`` must be
present in the tracefile or be proven by fresh gcov metadata to contain no
executable lines. No undeclared source may slip into the report.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
import json
import os
from pathlib import Path
import subprocess
import sys


class CoverageError(ValueError):
    """Raised when a tracefile or source manifest is incomplete or invalid."""


@dataclass
class SourceCoverage:
    lines: dict[int, int] = field(default_factory=dict)
    functions: set[tuple[int, str]] = field(default_factory=set)
    function_counts: dict[str, int] = field(default_factory=dict)
    function_lines: dict[str, int] = field(default_factory=dict)
    function_id_counts: dict[tuple[str, str], int] = field(default_factory=dict)
    functions_found: int | None = None
    functions_hit: int | None = None
    branches: dict[tuple[int, str, str], int] = field(default_factory=dict)


@dataclass(frozen=True)
class Metric:
    hit: int
    found: int

    @property
    def percent(self) -> float:
        if self.found == 0:
            raise CoverageError("coverage metric has an empty denominator")
        return 100.0 * self.hit / self.found


def _count(value: str, context: str) -> int:
    try:
        return int(value)
    except ValueError as error:
        raise CoverageError(f"invalid execution count in {context}: {value!r}") from error


def _merge_count(mapping: dict, key, count: int) -> None:
    mapping[key] = mapping.get(key, 0) + count


def parse_tracefile(path: Path, base_directory: Path) -> dict[Path, SourceCoverage]:
    """Read the executable-line, function, and branch records from an LCOV file."""
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise CoverageError(f"cannot read LCOV tracefile {path}: {error}") from error

    coverage: dict[Path, SourceCoverage] = {}
    current: SourceCoverage | None = None

    for line_number, raw_line in enumerate(lines, start=1):
        if not raw_line:
            continue
        if raw_line.startswith("SF:"):
            source_text = raw_line[3:]
            if not source_text:
                raise CoverageError(f"{path}:{line_number}: empty SF record")
            source_path = Path(source_text)
            if not source_path.is_absolute():
                source_path = base_directory / source_path
            current = coverage.setdefault(source_path.resolve(), SourceCoverage())
            continue
        if raw_line == "end_of_record":
            current = None
            continue
        if current is None:
            # TN and other tracefile metadata legitimately occur outside SF records.
            continue

        context = f"{path}:{line_number}"
        if raw_line.startswith("DA:"):
            fields = raw_line[3:].split(",", 2)
            if len(fields) < 2:
                raise CoverageError(f"{context}: malformed DA record")
            source_line = _count(fields[0], context)
            executions = _count(fields[1], context)
            _merge_count(current.lines, source_line, executions)
        elif raw_line.startswith("FN:"):
            fields = raw_line[3:].split(",", 2)
            if len(fields) < 2:
                raise CoverageError(f"{context}: malformed FN record")
            source_line = _count(fields[0], context)
            # LCOV 2.x may emit FN:start,end,name; older versions emit FN:start,name.
            name = fields[2] if len(fields) == 3 and fields[1].isdigit() else ",".join(fields[1:])
            if not name:
                raise CoverageError(f"{context}: empty function name")
            current.functions.add((source_line, name))
        elif raw_line.startswith("FNDA:"):
            fields = raw_line[5:].split(",", 1)
            if len(fields) != 2 or not fields[1]:
                raise CoverageError(f"{context}: malformed FNDA record")
            _merge_count(current.function_counts, fields[1], _count(fields[0], context))
        elif raw_line.startswith("FNL:"):
            fields = raw_line[4:].split(",", 2)
            if len(fields) < 2 or not fields[0]:
                raise CoverageError(f"{context}: malformed FNL record")
            function_id = fields[0]
            source_line = _count(fields[1], context)
            previous = current.function_lines.setdefault(function_id, source_line)
            if previous != source_line:
                raise CoverageError(f"{context}: conflicting FNL record")
        elif raw_line.startswith("FNA:"):
            fields = raw_line[4:].split(",", 2)
            if len(fields) != 3 or not fields[0] or not fields[2]:
                raise CoverageError(f"{context}: malformed FNA record")
            _merge_count(
                current.function_id_counts,
                (fields[0], fields[2]),
                _count(fields[1], context),
            )
        elif raw_line.startswith("FNF:"):
            if current.functions_found is not None:
                raise CoverageError(f"{context}: duplicate FNF record")
            current.functions_found = _count(raw_line[4:], context)
        elif raw_line.startswith("FNH:"):
            if current.functions_hit is not None:
                raise CoverageError(f"{context}: duplicate FNH record")
            current.functions_hit = _count(raw_line[4:], context)
        elif raw_line.startswith("BRDA:"):
            fields = raw_line[5:].split(",", 3)
            if len(fields) != 4:
                raise CoverageError(f"{context}: malformed BRDA record")
            source_line = _count(fields[0], context)
            executions = 0 if fields[3] == "-" else _count(fields[3], context)
            _merge_count(current.branches, (source_line, fields[1], fields[2]), executions)

    if not coverage:
        raise CoverageError(f"LCOV tracefile contains no source records: {path}")
    return coverage


def parse_source_manifest(path: Path, variable: str = "SRC_FILES") -> list[str]:
    """Parse the simple backslash-continued source list owned by src/."""
    try:
        raw_lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise CoverageError(f"cannot read source manifest {path}: {error}") from error

    entries: list[str] = []
    collecting = False
    for raw_line in raw_lines:
        line = raw_line.split("#", 1)[0].strip()
        if not collecting:
            prefix = f"{variable} :="
            if not line.startswith(prefix):
                continue
            line = line[len(prefix) :].strip()
            collecting = True
        continued = line.endswith("\\")
        if continued:
            line = line[:-1].strip()
        entries.extend(line.split())
        if not continued:
            break

    if not collecting or not entries:
        raise CoverageError(f"{path} has no non-empty {variable} assignment")
    duplicates = sorted({entry for entry in entries if entries.count(entry) > 1})
    if duplicates:
        raise CoverageError(
            f"{path} lists duplicate {variable} entries: {', '.join(duplicates)}"
        )
    return entries


def expected_sources(manifest: Path, source_root: Path) -> list[Path]:
    root = source_root.resolve()
    paths: list[Path] = []
    for entry in parse_source_manifest(manifest):
        path = (root / entry).resolve()
        try:
            path.relative_to(root)
        except ValueError as error:
            raise CoverageError(f"manifest source escapes source root: {entry}") from error
        if not path.is_file():
            raise CoverageError(f"manifest source does not exist: {path}")
        paths.append(path)
    return paths


def source_metrics(source: SourceCoverage) -> dict[str, Metric]:
    if (source.functions_found is None) != (source.functions_hit is None):
        raise CoverageError("LCOV source record has only one of FNF/FNH")
    if source.function_id_counts:
        # LCOV 2.x can attach several generated Fortran procedures to one FNL
        # location ID.  Its own summary counts the distinct FNA records.
        functions = Metric(
            hit=sum(count > 0 for count in source.function_id_counts.values()),
            found=len(source.function_id_counts),
        )
    elif source.functions:
        functions = Metric(
            hit=sum(
                source.function_counts.get(name, 0) > 0
                for _, name in source.functions
            ),
            found=len(source.functions),
        )
    elif source.functions_found is not None:
        functions = Metric(
            hit=source.functions_hit,
            found=source.functions_found,
        )
        if functions.hit > functions.found:
            raise CoverageError("LCOV source record has FNH greater than FNF")
    else:
        functions = Metric(hit=0, found=0)
    return {
        "lines": Metric(
            hit=sum(count > 0 for count in source.lines.values()),
            found=len(source.lines),
        ),
        "functions": functions,
        "branches": Metric(
            hit=sum(count > 0 for count in source.branches.values()),
            found=len(source.branches),
        ),
    }


def total_metrics(coverage: dict[Path, SourceCoverage]) -> dict[str, Metric]:
    per_source = [source_metrics(item) for item in coverage.values()]
    return {
        name: Metric(
            hit=sum(metrics[name].hit for metrics in per_source),
            found=sum(metrics[name].found for metrics in per_source),
        )
        for name in ("lines", "functions", "branches")
    }


def verify_noninstrumentable_sources(
    missing: list[Path], object_root: Path, gcov: str
) -> list[Path]:
    """Require gcov metadata proving that every absent source has no code."""
    if not missing:
        return []
    if not object_root.is_dir():
        raise CoverageError(f"coverage object root does not exist: {object_root}")

    verified: list[Path] = []
    failures: list[str] = []
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    for source in missing:
        metadata = object_root / f"{source.stem}.gcno"
        if not metadata.is_file():
            failures.append(f"{source} (missing {metadata.name})")
            continue
        try:
            result = subprocess.run(
                [
                    gcov,
                    "--no-output",
                    "--branch-counts",
                    "--object-directory",
                    str(object_root),
                    str(source),
                ],
                check=False,
                capture_output=True,
                text=True,
                env=environment,
            )
        except OSError as error:
            raise CoverageError(f"cannot execute gcov tool {gcov}: {error}") from error
        diagnostic = result.stdout + result.stderr
        if result.returncode == 0 and "No executable lines" in diagnostic:
            verified.append(source)
        else:
            failures.append(str(source))

    if failures:
        raise CoverageError(
            "active sources missing from LCOV tracefile and not proven "
            "non-instrumentable by gcov: " + ", ".join(failures)
        )
    return verified


def _threshold(value: str) -> float:
    try:
        threshold = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(f"not a number: {value}") from error
    if not 0.0 <= threshold <= 100.0:
        raise argparse.ArgumentTypeError("must be between 0 and 100")
    return threshold


def _metric_payload(metric: Metric, threshold: float) -> dict[str, object]:
    return {
        "hit": metric.hit,
        "found": metric.found,
        "percent": round(metric.percent, 6),
        "minimum_percent": threshold,
        "passed": metric.percent >= threshold,
    }


def _relative(path: Path, base: Path) -> str:
    try:
        return path.relative_to(base).as_posix()
    except ValueError:
        return str(path)


def build_summary(
    coverage: dict[Path, SourceCoverage],
    expected: list[Path],
    noninstrumentable: list[Path],
    thresholds: dict[str, float],
    repository_root: Path,
) -> dict[str, object]:
    expected_set = set(expected)
    present_set = set(coverage)
    unexpected = sorted(present_set - expected_set)
    if unexpected:
        raise CoverageError(
            "undeclared sources present in LCOV tracefile: "
            + ", ".join(_relative(path, repository_root) for path in unexpected)
        )

    totals = total_metrics(coverage)
    for name, metric in totals.items():
        if metric.found == 0:
            raise CoverageError(f"LCOV tracefile contains no {name} records")

    files: dict[str, object] = {}
    for source_path in expected:
        if source_path not in coverage:
            files[_relative(source_path, repository_root)] = {
                "instrumentable": False,
                "lines": {"hit": 0, "found": 0, "percent": None},
                "functions": {"hit": 0, "found": 0, "percent": None},
                "branches": {"hit": 0, "found": 0, "percent": None},
            }
            continue
        metrics = source_metrics(coverage[source_path])
        files[_relative(source_path, repository_root)] = {
            "instrumentable": True,
            **{
                name: {
                    "hit": metric.hit,
                    "found": metric.found,
                    "percent": round(metric.percent, 6) if metric.found else None,
                }
                for name, metric in metrics.items()
            },
        }

    total_payload = {
        name: _metric_payload(metric, thresholds[name])
        for name, metric in totals.items()
    }
    return {
        "schema_version": 1,
        "passed": all(item["passed"] for item in total_payload.values()),
        "active_source_files": len(expected),
        "instrumentable_source_files": len(coverage),
        "noninstrumentable_sources": [
            _relative(path, repository_root) for path in noninstrumentable
        ],
        "totals": total_payload,
        "files": files,
    }


def write_json(path: Path, payload: dict[str, object]) -> None:
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + ".tmp")
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        temporary.replace(path)
    except OSError as error:
        raise CoverageError(f"cannot write coverage summary {path}: {error}") from error


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--tracefile", type=Path, required=True)
    result.add_argument("--source-manifest", type=Path, required=True)
    result.add_argument("--source-root", type=Path, required=True)
    result.add_argument("--repository-root", type=Path, required=True)
    result.add_argument("--object-root", type=Path, required=True)
    result.add_argument("--gcov", required=True)
    result.add_argument("--summary-json", type=Path, required=True)
    result.add_argument("--min-lines", type=_threshold, required=True)
    result.add_argument("--min-functions", type=_threshold, required=True)
    result.add_argument("--min-branches", type=_threshold, required=True)
    return result


def main(arguments: list[str] | None = None) -> int:
    options = parser().parse_args(arguments)
    repository_root = options.repository_root.resolve()
    thresholds = {
        "lines": options.min_lines,
        "functions": options.min_functions,
        "branches": options.min_branches,
    }
    try:
        coverage = parse_tracefile(options.tracefile, repository_root)
        expected = expected_sources(options.source_manifest, options.source_root)
        missing = sorted(set(expected) - set(coverage))
        noninstrumentable = verify_noninstrumentable_sources(
            missing, options.object_root.resolve(), options.gcov
        )
        summary = build_summary(
            coverage, expected, noninstrumentable, thresholds, repository_root
        )
        write_json(options.summary_json, summary)
    except CoverageError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2

    print(
        f"Coverage: {summary['active_source_files']} active src files "
        f"({summary['instrumentable_source_files']} instrumentable)"
    )
    if summary["noninstrumentable_sources"]:
        print(
            "  no executable lines: "
            + ", ".join(summary["noninstrumentable_sources"])
        )
    for name in ("lines", "functions", "branches"):
        item = summary["totals"][name]
        status = "PASS" if item["passed"] else "FAIL"
        print(
            f"  {name:9s} {item['percent']:8.3f}% "
            f"({item['hit']}/{item['found']}), minimum {item['minimum_percent']:.3f}% "
            f"[{status}]"
        )
    if not summary["passed"]:
        print("ERROR: one or more coverage thresholds were not met", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
