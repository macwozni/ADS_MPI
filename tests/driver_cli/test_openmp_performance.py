#!/usr/bin/env python3
"""Compare release iGRM-L2 runtime and output with one and four threads.

The MPI launcher and release executable are supplied through ``MPIEXEC`` and
``EXECUTABLE``.  Timing thresholds are intentionally relative: both OpenMP
variants are measured in alternating order during the same invocation.
"""

from __future__ import annotations

import argparse
import contextlib
from dataclasses import dataclass
import io
import json
import math
import os
from pathlib import Path
import platform
import shlex
import signal
import statistics
import subprocess
import sys
import tempfile
import time
import unittest
from unittest import mock
import xml.etree.ElementTree as ET


WORKLOAD_ARGUMENTS = (
    "32",
    "32",
    "32",
    "3",
    "3",
    "3",
    "2",
    "2",
    "2",
    "1",
    "1",
    "1",
    "dg",
)
FIELD_ABS_TOLERANCE = 1.0e-12
FIELD_REL_TOLERANCE = 1.0e-12
EXPECTED_POINT_COUNT = 31**3
RUNTIME_DIAGNOSTICS = (
    "ERROR: AddressSanitizer",
    "Fortran runtime error",
    "ERROR STOP",
    "Segmentation fault",
    "ThreadSanitizer: data race",
)


class PerformanceTestError(RuntimeError):
    """A configuration, execution, or numerical-validation failure."""


@dataclass(frozen=True)
class Configuration:
    mpi_command: tuple[str, ...]
    mpi_flags: tuple[str, ...]
    mpi_np_flag: str
    executable: Path
    warmups: int
    samples: int
    timeout: float
    minimum_speedup: float
    maximum_regression: float
    baseline: Path | None
    output: Path | None
    omp_proc_bind: str
    omp_places: str


@dataclass(frozen=True)
class Measurement:
    phase: str
    pair: int
    threads: int
    seconds: float
    values: tuple[float, ...]


@dataclass(frozen=True)
class GateResult:
    passed: bool
    failures: tuple[str, ...]
    omp1_median_seconds: float
    omp4_median_seconds: float
    paired_speedups: tuple[float, ...]
    median_speedup: float
    median_time_speedup: float
    faster_pairs: int
    required_faster_pairs: int
    baseline_ratios: dict[str, float]


def _environment_integer(name: str, default: int, minimum: int) -> int:
    raw = os.environ.get(name, str(default))
    try:
        value = int(raw)
    except ValueError as error:
        raise PerformanceTestError(f"{name} must be an integer, got {raw!r}") from error
    if value < minimum:
        raise PerformanceTestError(f"{name} must be >= {minimum}, got {value}")
    return value


def _warmup_count() -> int:
    """Read the canonical plural setting, accepting the old singular alias."""
    if "PERFORMANCE_WARMUPS" in os.environ:
        return _environment_integer("PERFORMANCE_WARMUPS", 1, 0)
    return _environment_integer("PERFORMANCE_WARMUP", 1, 0)


def _environment_float(name: str, default: float) -> float:
    raw = os.environ.get(name, str(default))
    try:
        value = float(raw)
    except ValueError as error:
        raise PerformanceTestError(f"{name} must be a number, got {raw!r}") from error
    if not math.isfinite(value) or value <= 0.0:
        raise PerformanceTestError(f"{name} must be finite and positive, got {raw!r}")
    return value


def configuration_from_environment() -> Configuration:
    mpi_raw = os.environ.get("MPIEXEC", "").strip()
    executable_raw = os.environ.get("EXECUTABLE", "").strip()
    if not mpi_raw:
        raise PerformanceTestError("MPIEXEC is required")
    if not executable_raw:
        raise PerformanceTestError("EXECUTABLE is required")

    mpi_command = tuple(shlex.split(mpi_raw))
    if not mpi_command:
        raise PerformanceTestError("MPIEXEC must name an MPI launcher")
    mpi_flags = tuple(shlex.split(os.environ.get("MPIEXEC_FLAGS", "")))
    mpi_np_flag = os.environ.get("MPI_NP_FLAG", "-n").strip()
    if not mpi_np_flag or any(character.isspace() for character in mpi_np_flag):
        raise PerformanceTestError("MPI_NP_FLAG must be one nonempty argument")

    omp_proc_bind = os.environ.get("OMP_PROC_BIND", "close").strip()
    omp_places = os.environ.get("OMP_PLACES", "cores").strip()
    if not omp_proc_bind:
        raise PerformanceTestError("OMP_PROC_BIND must not be empty")
    if not omp_places:
        raise PerformanceTestError("OMP_PLACES must not be empty")

    executable = Path(executable_raw).expanduser().resolve()
    if not executable.is_file():
        raise PerformanceTestError(f"EXECUTABLE is not a file: {executable}")
    if not os.access(executable, os.X_OK):
        raise PerformanceTestError(f"EXECUTABLE is not executable: {executable}")

    baseline_raw = os.environ.get("PERFORMANCE_BASELINE", "").strip()
    output_raw = os.environ.get("PERFORMANCE_OUTPUT", "").strip()
    baseline = Path(baseline_raw).expanduser().resolve() if baseline_raw else None
    output = Path(output_raw).expanduser().resolve() if output_raw else None

    return Configuration(
        mpi_command=mpi_command,
        mpi_flags=mpi_flags,
        mpi_np_flag=mpi_np_flag,
        executable=executable,
        warmups=_warmup_count(),
        samples=_environment_integer("PERFORMANCE_SAMPLES", 3, 1),
        timeout=_environment_float("PERFORMANCE_TIMEOUT", 300.0),
        minimum_speedup=_environment_float("PERFORMANCE_MIN_SPEEDUP", 1.10),
        maximum_regression=_environment_float(
            "PERFORMANCE_MAX_REGRESSION", 1.15
        ),
        baseline=baseline,
        output=output,
        omp_proc_bind=omp_proc_bind,
        omp_places=omp_places,
    )


def _terminate_process_group(process: subprocess.Popen[str]) -> str:
    if process.poll() is None:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    try:
        output, _ = process.communicate(timeout=5)
        return output
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        output, _ = process.communicate()
        return output


def _output_tail(output: str, lines: int = 40) -> str:
    selected = output.rstrip().splitlines()[-lines:]
    return "\n".join(f"  {line}" for line in selected)


def _runtime_diagnostics(output: str) -> list[str]:
    lowered = output.lower()
    return [item for item in RUNTIME_DIAGNOSTICS if item.lower() in lowered]


def read_result_values(path: Path) -> tuple[float, ...]:
    if not path.is_file():
        raise PerformanceTestError(f"missing VTK result: {path}")
    try:
        root = ET.parse(path).getroot()
    except (ET.ParseError, OSError) as error:
        raise PerformanceTestError(f"cannot parse {path}: {error}") from error

    result_arrays = [
        element
        for element in root.iter()
        if element.tag.rsplit("}", 1)[-1] == "DataArray"
        and element.attrib.get("Name") == "Result"
    ]
    if len(result_arrays) != 1:
        raise PerformanceTestError(
            f"{path} must contain exactly one Result DataArray, found "
            f"{len(result_arrays)}"
        )

    result = result_arrays[0]
    if result.attrib.get("type") != "Float64":
        raise PerformanceTestError(
            f"{path} Result must use Float64, got {result.attrib.get('type')!r}"
        )
    if result.attrib.get("format", "ascii").lower() != "ascii":
        raise PerformanceTestError(f"{path} Result must use ASCII VTK storage")

    tokens = "".join(result.itertext()).split()
    if not tokens:
        raise PerformanceTestError(f"{path} Result is empty")
    try:
        values = tuple(float(token.replace("D", "E").replace("d", "e")) for token in tokens)
    except ValueError as error:
        raise PerformanceTestError(f"{path} Result contains a non-number") from error
    if not all(math.isfinite(value) for value in values):
        raise PerformanceTestError(f"{path} Result contains a non-finite value")
    if len(values) != EXPECTED_POINT_COUNT:
        raise PerformanceTestError(
            f"{path} Result has {len(values)} values, expected {EXPECTED_POINT_COUNT}"
        )
    return values


def compare_fields(
    left: tuple[float, ...],
    right: tuple[float, ...],
    label: str,
) -> float:
    if len(left) != len(right):
        raise PerformanceTestError(
            f"{label}: Result sizes differ ({len(left)} != {len(right)})"
        )

    maximum_difference = 0.0
    worst_index = 0
    for index, (left_value, right_value) in enumerate(zip(left, right)):
        difference = abs(left_value - right_value)
        maximum_difference = max(maximum_difference, difference)
        allowed = FIELD_ABS_TOLERANCE + FIELD_REL_TOLERANCE * max(
            abs(left_value), abs(right_value)
        )
        if difference > allowed:
            worst_index = index
            raise PerformanceTestError(
                f"{label}: Result differs at index {worst_index}: "
                f"left={left_value:.17g}, right={right_value:.17g}, "
                f"abs diff={difference:.3e}, allowed={allowed:.3e}"
            )
    return maximum_difference


def run_once(
    configuration: Configuration,
    work_root: Path,
    phase: str,
    pair: int,
    threads: int,
) -> Measurement:
    run_directory = work_root / f"{phase}-{pair:02d}-omp{threads}"
    run_directory.mkdir()
    command = [
        *configuration.mpi_command,
        *configuration.mpi_flags,
        configuration.mpi_np_flag,
        "1",
        str(configuration.executable),
        *WORKLOAD_ARGUMENTS,
    ]
    environment = os.environ.copy()
    environment.update(
        {
            "OMP_DYNAMIC": "FALSE",
            "OMP_NUM_THREADS": str(threads),
            "OMP_PROC_BIND": configuration.omp_proc_bind,
            "OMP_PLACES": configuration.omp_places,
        }
    )

    started = time.perf_counter()
    try:
        process = subprocess.Popen(
            command,
            cwd=run_directory,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            start_new_session=True,
        )
    except OSError as error:
        raise PerformanceTestError(
            f"cannot start {phase} pair {pair} OMP{threads}: {error}"
        ) from error
    timed_out = False
    try:
        output, _ = process.communicate(timeout=configuration.timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        output = _terminate_process_group(process)
    except BaseException:
        _terminate_process_group(process)
        raise
    elapsed = time.perf_counter() - started
    (run_directory / "driver.log").write_text(output, encoding="utf-8")

    if timed_out:
        raise PerformanceTestError(
            f"{phase} pair {pair} OMP{threads} timed out after "
            f"{configuration.timeout:g}s\n{_output_tail(output)}"
        )
    if process.returncode != 0:
        raise PerformanceTestError(
            f"{phase} pair {pair} OMP{threads} returned {process.returncode}\n"
            f"{_output_tail(output)}"
        )
    diagnostics = _runtime_diagnostics(output)
    if diagnostics:
        raise PerformanceTestError(
            f"{phase} pair {pair} OMP{threads} emitted runtime diagnostics: "
            f"{', '.join(diagnostics)}\n{_output_tail(output)}"
        )

    values = read_result_values(run_directory / "step0.vti")
    if not math.isfinite(elapsed) or elapsed <= 0.0:
        raise PerformanceTestError(
            f"{phase} pair {pair} OMP{threads} produced invalid time {elapsed!r}"
        )
    print(f"RUN {phase} pair={pair} OMP{threads} seconds={elapsed:.9f}")
    return Measurement(phase, pair, threads, elapsed, values)


def load_baseline(
    path: Path | None,
    configuration: Configuration,
) -> dict[str, float] | None:
    if path is None:
        return None
    if not path.is_file():
        raise PerformanceTestError(f"PERFORMANCE_BASELINE does not exist: {path}")
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise PerformanceTestError(f"cannot read baseline {path}: {error}") from error
    if not isinstance(document, dict):
        raise PerformanceTestError(f"baseline {path} must contain a JSON object")
    if document.get("schema_version") != 1:
        raise PerformanceTestError(
            f"baseline {path} has unsupported schema_version "
            f"{document.get('schema_version')!r}"
        )

    workload = document.get("workload")
    expected_workload = {
        "driver": "igrm_l2",
        "arguments": list(WORKLOAD_ARGUMENTS),
        "mpi_ranks": 1,
        "openmp_threads": [1, 4],
    }
    if workload != expected_workload:
        raise PerformanceTestError(
            f"baseline {path} was produced for an incompatible workload"
        )

    baseline_environment = document.get("environment")
    if not isinstance(baseline_environment, dict):
        raise PerformanceTestError(f"baseline {path} has an invalid environment")
    for name, expected in (
        ("mpiexec", list(configuration.mpi_command)),
        ("mpiexec_flags", list(configuration.mpi_flags)),
        ("mpi_np_flag", configuration.mpi_np_flag),
        ("omp_proc_bind", configuration.omp_proc_bind),
        ("omp_places", configuration.omp_places),
    ):
        if baseline_environment.get(name) != expected:
            raise PerformanceTestError(
                f"baseline {path} uses incompatible {name}: "
                f"{baseline_environment.get(name)!r} != {expected!r}"
            )

    summary = document.get("summary")
    if not isinstance(summary, dict):
        raise PerformanceTestError(f"baseline {path} has an invalid summary")

    baseline: dict[str, float] = {}
    for name in ("omp1_median_seconds", "omp4_median_seconds"):
        try:
            value = float(summary[name])
        except (KeyError, TypeError, ValueError) as error:
            raise PerformanceTestError(
                f"baseline {path} is missing a numeric {name}"
            ) from error
        if not math.isfinite(value) or value <= 0.0:
            raise PerformanceTestError(
                f"baseline {path} has invalid {name}: {value!r}"
            )
        baseline[name] = value
    return baseline


def evaluate_performance(
    omp1_seconds: list[float],
    omp4_seconds: list[float],
    minimum_speedup: float,
    baseline: dict[str, float] | None = None,
    maximum_regression: float = 1.15,
) -> GateResult:
    if not omp1_seconds or len(omp1_seconds) != len(omp4_seconds):
        raise ValueError("OMP1 and OMP4 need the same nonzero number of samples")
    if minimum_speedup <= 0.0 or maximum_regression <= 0.0:
        raise ValueError("performance thresholds must be positive")
    if not all(
        math.isfinite(value) and value > 0.0
        for value in (*omp1_seconds, *omp4_seconds)
    ):
        raise ValueError("all measurements must be finite and positive")

    omp1_median = statistics.median(omp1_seconds)
    omp4_median = statistics.median(omp4_seconds)
    paired_speedups = tuple(
        omp1 / omp4 for omp1, omp4 in zip(omp1_seconds, omp4_seconds)
    )
    median_speedup = statistics.median(paired_speedups)
    median_time_speedup = omp1_median / omp4_median
    faster_pairs = sum(
        omp4 < omp1 for omp1, omp4 in zip(omp1_seconds, omp4_seconds)
    )
    required_faster_pairs = len(omp1_seconds) // 2 + 1
    failures: list[str] = []

    if median_speedup < minimum_speedup:
        failures.append(
            f"median paired speedup {median_speedup:.6f} is below "
            f"{minimum_speedup:.6f}"
        )
    if faster_pairs < required_faster_pairs:
        failures.append(
            f"OMP4 was faster in only {faster_pairs}/{len(omp1_seconds)} pairs; "
            f"at least {required_faster_pairs} are required"
        )

    baseline_ratios: dict[str, float] = {}
    if baseline is not None:
        for name, current in (
            ("omp1_median_seconds", omp1_median),
            ("omp4_median_seconds", omp4_median),
        ):
            if name not in baseline:
                raise ValueError(f"baseline is missing {name}")
            baseline_value = baseline[name]
            if not math.isfinite(baseline_value) or baseline_value <= 0.0:
                raise ValueError(f"baseline {name} must be finite and positive")
            ratio = current / baseline_value
            baseline_ratios[name] = ratio
            if ratio > maximum_regression:
                failures.append(
                    f"{name} regression factor {ratio:.6f} exceeds "
                    f"{maximum_regression:.6f}"
                )

    return GateResult(
        passed=not failures,
        failures=tuple(failures),
        omp1_median_seconds=omp1_median,
        omp4_median_seconds=omp4_median,
        paired_speedups=paired_speedups,
        median_speedup=median_speedup,
        median_time_speedup=median_time_speedup,
        faster_pairs=faster_pairs,
        required_faster_pairs=required_faster_pairs,
        baseline_ratios=baseline_ratios,
    )


def make_output_document(
    configuration: Configuration,
    omp1_seconds: list[float],
    omp4_seconds: list[float],
    gate: GateResult,
    maximum_field_difference: float,
) -> dict[str, object]:
    return {
        "schema_version": 1,
        "passed": gate.passed,
        "failures": list(gate.failures),
        "environment": {
            "executable": str(configuration.executable),
            "mpiexec": list(configuration.mpi_command),
            "mpiexec_flags": list(configuration.mpi_flags),
            "mpi_np_flag": configuration.mpi_np_flag,
            "omp_proc_bind": configuration.omp_proc_bind,
            "omp_places": configuration.omp_places,
            "platform": platform.platform(),
            "python": platform.python_version(),
        },
        "workload": {
            "driver": "igrm_l2",
            "arguments": list(WORKLOAD_ARGUMENTS),
            "mpi_ranks": 1,
            "openmp_threads": [1, 4],
        },
        "settings": {
            "warmups": configuration.warmups,
            "samples": configuration.samples,
            "timeout_seconds": configuration.timeout,
            "minimum_speedup": configuration.minimum_speedup,
            "maximum_regression": configuration.maximum_regression,
            "baseline": str(configuration.baseline) if configuration.baseline else None,
        },
        "samples": [
            {
                "pair": index,
                "omp1_seconds": omp1,
                "omp4_seconds": omp4,
                "speedup": omp1 / omp4,
            }
            for index, (omp1, omp4) in enumerate(
                zip(omp1_seconds, omp4_seconds), start=1
            )
        ],
        "summary": {
            "omp1_median_seconds": gate.omp1_median_seconds,
            "omp4_median_seconds": gate.omp4_median_seconds,
            "median_paired_speedup": gate.median_speedup,
            "median_time_speedup": gate.median_time_speedup,
            "faster_pairs": gate.faster_pairs,
            "required_faster_pairs": gate.required_faster_pairs,
            "baseline_regression_factors": gate.baseline_ratios,
            "maximum_field_difference": maximum_field_difference,
        },
    }


def write_output(path: Path | None, document: dict[str, object]) -> None:
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"WROTE {path}")


def remove_stale_output(path: Path | None) -> None:
    if path is None:
        return
    try:
        path.unlink(missing_ok=True)
    except OSError as error:
        raise PerformanceTestError(
            f"cannot remove stale performance report {path}: {error}"
        ) from error


def run_performance_test(configuration: Configuration) -> int:
    same_baseline_and_output = (
        configuration.baseline is not None
        and configuration.output is not None
        and configuration.baseline == configuration.output
    )
    if same_baseline_and_output:
        try:
            baseline = load_baseline(configuration.baseline, configuration)
        finally:
            remove_stale_output(configuration.output)
    else:
        remove_stale_output(configuration.output)
        baseline = load_baseline(configuration.baseline, configuration)
    omp1_seconds: list[float] = []
    omp4_seconds: list[float] = []
    maximum_field_difference = 0.0
    reference_values: tuple[float, ...] | None = None
    reference_label = ""

    with tempfile.TemporaryDirectory(prefix="mpi-ads-openmp-performance.") as temp:
        work_root = Path(temp)
        total_pairs = configuration.warmups + configuration.samples
        for pair_index in range(total_pairs):
            phase = "warmup" if pair_index < configuration.warmups else "sample"
            phase_pair = (
                pair_index + 1
                if phase == "warmup"
                else pair_index - configuration.warmups + 1
            )
            order = (1, 4) if pair_index % 2 == 0 else (4, 1)
            pair_measurements: dict[int, Measurement] = {}
            for threads in order:
                measurement = run_once(
                    configuration,
                    work_root,
                    phase,
                    phase_pair,
                    threads,
                )
                pair_measurements[threads] = measurement
                measurement_label = f"{phase} pair {phase_pair} OMP{threads}"
                if reference_values is None:
                    reference_values = measurement.values
                    reference_label = measurement_label
                else:
                    difference = compare_fields(
                        reference_values,
                        measurement.values,
                        f"{reference_label} vs {measurement_label}",
                    )
                    maximum_field_difference = max(
                        maximum_field_difference, difference
                    )

            difference = compare_fields(
                pair_measurements[1].values,
                pair_measurements[4].values,
                f"{phase} pair {phase_pair}",
            )
            maximum_field_difference = max(maximum_field_difference, difference)
            print(
                f"FIELD {phase} pair={phase_pair} values="
                f"{len(pair_measurements[1].values)} max_abs_diff={difference:.3e}"
            )
            if phase == "sample":
                omp1_seconds.append(pair_measurements[1].seconds)
                omp4_seconds.append(pair_measurements[4].seconds)

    gate = evaluate_performance(
        omp1_seconds,
        omp4_seconds,
        configuration.minimum_speedup,
        baseline,
        configuration.maximum_regression,
    )
    print(f"OMP1 median seconds: {gate.omp1_median_seconds:.9f}")
    print(f"OMP4 median seconds: {gate.omp4_median_seconds:.9f}")
    print(f"median paired speedup: {gate.median_speedup:.6f}")
    print(f"median-time speedup: {gate.median_time_speedup:.6f}")
    print(
        f"OMP4 faster pairs: {gate.faster_pairs}/{configuration.samples} "
        f"(required {gate.required_faster_pairs})"
    )
    for name, ratio in gate.baseline_ratios.items():
        print(f"baseline regression {name}: {ratio:.6f}")

    document = make_output_document(
        configuration,
        omp1_seconds,
        omp4_seconds,
        gate,
        maximum_field_difference,
    )
    write_output(configuration.output, document)
    if gate.passed:
        print("PASS OpenMP performance regression")
        return 0
    for failure in gate.failures:
        print(f"FAIL {failure}")
    return 1


class PerformanceSelfTests(unittest.TestCase):
    @staticmethod
    def _configuration(root: Path, **overrides: object) -> Configuration:
        executable = root / "igrm_l2"
        executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        executable.chmod(0o755)
        values: dict[str, object] = {
            "mpi_command": ("mpiexec",),
            "mpi_flags": (),
            "mpi_np_flag": "-n",
            "executable": executable.resolve(),
            "warmups": 1,
            "samples": 3,
            "timeout": 300.0,
            "minimum_speedup": 1.10,
            "maximum_regression": 1.15,
            "baseline": None,
            "output": None,
            "omp_proc_bind": "close",
            "omp_places": "cores",
        }
        values.update(overrides)
        return Configuration(**values)  # type: ignore[arg-type]

    @staticmethod
    def _baseline_document(configuration: Configuration) -> dict[str, object]:
        return {
            "schema_version": 1,
            "environment": {
                "mpiexec": list(configuration.mpi_command),
                "mpiexec_flags": list(configuration.mpi_flags),
                "mpi_np_flag": configuration.mpi_np_flag,
                "omp_proc_bind": configuration.omp_proc_bind,
                "omp_places": configuration.omp_places,
            },
            "workload": {
                "driver": "igrm_l2",
                "arguments": list(WORKLOAD_ARGUMENTS),
                "mpi_ranks": 1,
                "openmp_threads": [1, 4],
            },
            "summary": {
                "omp1_median_seconds": 4.0,
                "omp4_median_seconds": 1.0,
            },
        }

    @staticmethod
    def _write_vti(path: Path, values: list[float | str]) -> None:
        payload = " ".join(str(value) for value in values)
        path.write_text(
            '<?xml version="1.0"?>\n'
            '<VTKFile type="ImageData">\n'
            "<ImageData><Piece><PointData>\n"
            '<DataArray type="Float64" Name="Result" format="ascii">\n'
            f"{payload}\n"
            "</DataArray></PointData></Piece></ImageData></VTKFile>\n",
            encoding="utf-8",
        )

    def test_speedup_and_majority_pass(self) -> None:
        result = evaluate_performance([4.0, 4.2, 3.8], [2.0, 2.1, 1.9], 1.10)
        self.assertTrue(result.passed)
        self.assertEqual(result.faster_pairs, 3)
        self.assertAlmostEqual(result.median_speedup, 2.0)

    def test_low_median_speedup_fails(self) -> None:
        result = evaluate_performance([1.0, 1.0, 1.0], [0.95, 0.95, 0.95], 1.10)
        self.assertFalse(result.passed)
        self.assertTrue(any("median paired speedup" in item for item in result.failures))

    def test_missing_majority_fails(self) -> None:
        result = evaluate_performance([2.0, 1.0, 1.0], [1.0, 2.0, 2.0], 0.49)
        self.assertFalse(result.passed)
        self.assertTrue(any("faster in only" in item for item in result.failures))

    def test_baseline_regression_fails(self) -> None:
        baseline = {
            "omp1_median_seconds": 1.0,
            "omp4_median_seconds": 0.5,
        }
        result = evaluate_performance(
            [1.20, 1.20, 1.20],
            [0.60, 0.60, 0.60],
            1.10,
            baseline,
            1.15,
        )
        self.assertFalse(result.passed)
        self.assertEqual(result.baseline_ratios["omp1_median_seconds"], 1.20)
        self.assertEqual(result.baseline_ratios["omp4_median_seconds"], 1.20)
        self.assertEqual(sum("regression factor" in item for item in result.failures), 2)

    def test_baseline_limit_is_inclusive(self) -> None:
        baseline = {
            "omp1_median_seconds": 1.0,
            "omp4_median_seconds": 0.5,
        }
        result = evaluate_performance(
            [1.15, 1.15, 1.15],
            [0.575, 0.575, 0.575],
            1.10,
            baseline,
            1.15,
        )
        self.assertTrue(result.passed)

    def test_invalid_samples_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_performance([1.0], [], 1.10)
        with self.assertRaises(ValueError):
            evaluate_performance([math.nan], [1.0], 1.10)

    def test_configuration_parses_all_forwarded_settings(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            executable = self._configuration(root).executable
            environment = {
                "MPIEXEC": "launcher --launcher-option",
                "MPIEXEC_FLAGS": "--bind-to core",
                "MPI_NP_FLAG": "--ranks",
                "EXECUTABLE": str(executable),
                "PERFORMANCE_WARMUPS": "2",
                "PERFORMANCE_SAMPLES": "5",
                "PERFORMANCE_TIMEOUT": "17.5",
                "PERFORMANCE_MIN_SPEEDUP": "1.25",
                "PERFORMANCE_MAX_REGRESSION": "1.05",
                "PERFORMANCE_BASELINE": "baseline.json",
                "PERFORMANCE_OUTPUT": "result.json",
                "OMP_PROC_BIND": "spread",
                "OMP_PLACES": "threads",
            }
            with mock.patch.dict(os.environ, environment, clear=True):
                configuration = configuration_from_environment()
            self.assertEqual(configuration.mpi_command, ("launcher", "--launcher-option"))
            self.assertEqual(configuration.mpi_flags, ("--bind-to", "core"))
            self.assertEqual(configuration.mpi_np_flag, "--ranks")
            self.assertEqual(configuration.warmups, 2)
            self.assertEqual(configuration.samples, 5)
            self.assertEqual(configuration.timeout, 17.5)
            self.assertEqual(configuration.minimum_speedup, 1.25)
            self.assertEqual(configuration.maximum_regression, 1.05)
            self.assertEqual(configuration.omp_proc_bind, "spread")
            self.assertEqual(configuration.omp_places, "threads")
            self.assertEqual(configuration.baseline, (Path.cwd() / "baseline.json").resolve())
            self.assertEqual(configuration.output, (Path.cwd() / "result.json").resolve())

    def test_configuration_rejects_invalid_environment(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            executable = self._configuration(root).executable
            base = {"MPIEXEC": "mpiexec", "EXECUTABLE": str(executable)}
            invalid_cases = (
                ({"MPIEXEC": ""}, "MPIEXEC is required"),
                ({"EXECUTABLE": ""}, "EXECUTABLE is required"),
                ({"MPI_NP_FLAG": "two words"}, "MPI_NP_FLAG"),
                ({"PERFORMANCE_SAMPLES": "0"}, "PERFORMANCE_SAMPLES"),
                ({"PERFORMANCE_WARMUPS": "-1"}, "PERFORMANCE_WARMUPS"),
                ({"PERFORMANCE_TIMEOUT": "nan"}, "PERFORMANCE_TIMEOUT"),
                ({"OMP_PROC_BIND": ""}, "OMP_PROC_BIND"),
                ({"OMP_PLACES": ""}, "OMP_PLACES"),
            )
            for changes, message in invalid_cases:
                environment = dict(base)
                environment.update(changes)
                with self.subTest(changes=changes):
                    with mock.patch.dict(os.environ, environment, clear=True):
                        with self.assertRaisesRegex(PerformanceTestError, message):
                            configuration_from_environment()

    def test_baseline_requires_matching_schema_workload_and_environment(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            configuration = self._configuration(root)
            baseline_path = root / "baseline.json"
            document = self._baseline_document(configuration)
            baseline_path.write_text(json.dumps(document), encoding="utf-8")
            self.assertEqual(
                load_baseline(baseline_path, configuration),
                {"omp1_median_seconds": 4.0, "omp4_median_seconds": 1.0},
            )

            incompatible_documents = []
            for mutator in (
                lambda item: item.update(schema_version=2),
                lambda item: item["workload"].update(driver="heat"),
                lambda item: item["environment"].update(mpiexec=["other-mpiexec"]),
                lambda item: item["environment"].update(
                    mpiexec_flags=["--bind-to", "none"]
                ),
                lambda item: item["environment"].update(mpi_np_flag="--np"),
                lambda item: item["environment"].update(omp_places="sockets"),
            ):
                candidate = self._baseline_document(configuration)
                mutator(candidate)
                incompatible_documents.append(candidate)
            for index, candidate in enumerate(incompatible_documents):
                baseline_path.write_text(json.dumps(candidate), encoding="utf-8")
                with self.subTest(incompatible=index):
                    with self.assertRaises(PerformanceTestError):
                        load_baseline(baseline_path, configuration)

    def test_baseline_rejects_missing_malformed_and_invalid_metrics(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            configuration = self._configuration(root)
            baseline_path = root / "baseline.json"
            with self.assertRaisesRegex(PerformanceTestError, "does not exist"):
                load_baseline(baseline_path, configuration)
            baseline_path.write_text("{", encoding="utf-8")
            with self.assertRaisesRegex(PerformanceTestError, "cannot read baseline"):
                load_baseline(baseline_path, configuration)
            document = self._baseline_document(configuration)
            document["summary"]["omp4_median_seconds"] = math.inf
            baseline_path.write_text(json.dumps(document), encoding="utf-8")
            with self.assertRaisesRegex(PerformanceTestError, "invalid omp4"):
                load_baseline(baseline_path, configuration)

    def test_vti_parser_and_field_comparison_cover_full_result(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "result.vti"
            values = [float(index % 17) for index in range(EXPECTED_POINT_COUNT)]
            self._write_vti(path, values)
            parsed = read_result_values(path)
            self.assertEqual(len(parsed), EXPECTED_POINT_COUNT)
            self.assertEqual(compare_fields(parsed, tuple(values), "same"), 0.0)

            changed = list(parsed)
            changed[-1] += 1.0e-3
            with self.assertRaisesRegex(PerformanceTestError, "differs at index"):
                compare_fields(parsed, tuple(changed), "changed")

    def test_vti_parser_rejects_bad_count_and_nonfinite_values(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "result.vti"
            self._write_vti(path, [0.0])
            with self.assertRaisesRegex(PerformanceTestError, "expected"):
                read_result_values(path)
            self._write_vti(path, ["nan"] * EXPECTED_POINT_COUNT)
            with self.assertRaisesRegex(PerformanceTestError, "non-finite"):
                read_result_values(path)

    def test_output_document_has_stable_schema_and_complete_samples(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            configuration = self._configuration(root)
            gate = evaluate_performance([4.0, 4.2, 3.8], [1.0, 1.1, 0.9], 1.10)
            document = make_output_document(
                configuration,
                [4.0, 4.2, 3.8],
                [1.0, 1.1, 0.9],
                gate,
                0.0,
            )
            self.assertEqual(document["schema_version"], 1)
            self.assertEqual(document["workload"]["arguments"], list(WORKLOAD_ARGUMENTS))
            self.assertEqual(len(document["samples"]), 3)
            self.assertEqual(document["summary"]["maximum_field_difference"], 0.0)
            json.dumps(document, sort_keys=True)

    def test_all_runs_are_compared_with_one_common_field_reference(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            configuration = self._configuration(root, warmups=0, samples=2)
            reference = (0.0,) * EXPECTED_POINT_COUNT
            changed_values = list(reference)
            changed_values[-1] = 1.0e-3
            changed = tuple(changed_values)
            measurements = (
                Measurement("sample", 1, 1, 4.0, reference),
                Measurement("sample", 1, 4, 1.0, reference),
                Measurement("sample", 2, 4, 1.0, changed),
                Measurement("sample", 2, 1, 4.0, changed),
            )
            with mock.patch(f"{__name__}.run_once", side_effect=measurements):
                with contextlib.redirect_stdout(io.StringIO()):
                    with self.assertRaisesRegex(
                        PerformanceTestError,
                        "sample pair 1 OMP1 vs sample pair 2 OMP4",
                    ):
                        run_performance_test(configuration)

    def test_failed_run_removes_stale_output_report(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            output = root / "openmp-performance.json"
            output.write_text("stale\n", encoding="utf-8")
            configuration = self._configuration(
                root,
                warmups=0,
                samples=1,
                output=output,
            )
            with mock.patch(
                f"{__name__}.run_once",
                side_effect=PerformanceTestError("synthetic launch failure"),
            ):
                with self.assertRaisesRegex(PerformanceTestError, "synthetic"):
                    run_performance_test(configuration)
            self.assertFalse(output.exists())

    def test_invalid_external_baseline_removes_stale_output_report(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            baseline = root / "baseline.json"
            baseline.write_text("{", encoding="utf-8")
            output = root / "openmp-performance.json"
            output.write_text("stale\n", encoding="utf-8")
            configuration = self._configuration(
                root,
                baseline=baseline,
                output=output,
            )
            with self.assertRaisesRegex(PerformanceTestError, "cannot read baseline"):
                run_performance_test(configuration)
            self.assertTrue(baseline.exists())
            self.assertFalse(output.exists())

    def test_baseline_can_share_output_path_without_becoming_stale(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            output = root / "openmp-performance.json"
            configuration = self._configuration(
                root,
                warmups=0,
                samples=1,
                baseline=output,
                output=output,
            )
            output.write_text(
                json.dumps(self._baseline_document(configuration)),
                encoding="utf-8",
            )
            with mock.patch(
                f"{__name__}.run_once",
                side_effect=PerformanceTestError("synthetic launch failure"),
            ):
                with self.assertRaisesRegex(PerformanceTestError, "synthetic"):
                    run_performance_test(configuration)
            self.assertFalse(output.exists())

    def test_stale_output_removal_error_is_reported(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "openmp-performance.json"
            output.write_text("stale\n", encoding="utf-8")
            with mock.patch.object(
                Path,
                "unlink",
                side_effect=PermissionError("synthetic permission failure"),
            ):
                with self.assertRaisesRegex(
                    PerformanceTestError, "cannot remove stale performance report"
                ):
                    remove_stale_output(output)


def run_self_tests() -> int:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(PerformanceSelfTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="test performance-runner logic without launching MPI",
    )
    arguments = parser.parse_args(argv)
    if arguments.self_test:
        return run_self_tests()
    try:
        return run_performance_test(configuration_from_environment())
    except PerformanceTestError as error:
        print(f"FAIL OpenMP performance regression: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
