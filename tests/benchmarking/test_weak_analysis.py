from __future__ import annotations

import csv
from copy import deepcopy
from fractions import Fraction
import io
import math
from pathlib import Path
import tempfile
import unittest

from ads_benchmark.analysis.base import AnalysisError
from ads_benchmark.analysis.io import write_report
from ads_benchmark.analysis.weak import WeakScalingAnalyzer
from ads_benchmark.framework.model import (
    CaseSpec,
    MeasurementSpec,
    MpiSpec,
    PlannedCase,
    SamplingSpec,
    TimeSpec,
    WeakScalingSpec,
)


def _case(
    local: tuple[int, int, int],
    grid: tuple[int, int, int],
    threads: int,
    *,
    role: str = "measurement",
) -> PlannedCase:
    mesh = tuple(
        elements * processes
        for elements, processes in zip(local, grid, strict=True)
    )
    spec = CaseSpec(
        family="weak",
        problem="igrm_l2",
        scheme="dg",
        exact_case="spatial-cosine",
        time=TimeSpec(final_time="0.1", time_step="0.025", steps=4),
        mesh=mesh,  # type: ignore[arg-type]
        test_degree=(3, 3, 3),
        trial_degree=(2, 2, 2),
        mpi=MpiSpec(ranks=math.prod(grid), process_grid=grid),
        openmp_threads=threads,
        sampling=SamplingSpec(points_per_axis=3, write_samples=True),
        measurement=MeasurementSpec(
            warmups=2,
            samples=7,
            timeout_seconds="60",
            minimum_sample_seconds="0.001",
        ),
        build_profile="release",
        launcher="mpi",
        openmp_dynamic=False,
        openmp_proc_bind="close",
        openmp_places="cores",
        weak_scaling=WeakScalingSpec(
            local_elements=local,
            workload_basis="per-rank",
            role=role,
        ),
    )
    return PlannedCase(
        case_id=(
            f"igrm_l2-dg-local{'x'.join(map(str, local))}-"
            f"grid{'x'.join(map(str, grid))}-t{threads}-{role}"
        ),
        spec=spec,
    )


def _checksum(values: list[float]) -> float:
    ordinary = math.fsum(values)
    weighted = math.fsum(
        (index + 1) * value for index, value in enumerate(values)
    )
    return ordinary + weighted / (len(values) + 1)


def _result(
    case: PlannedCase,
    median_seconds: float,
    *,
    short: bool = False,
    perturb: bool = False,
) -> dict[str, object]:
    final_time = float(Fraction(case.spec.time.final_time))
    coordinates: list[tuple[float, float, float]] = []
    exact: list[float] = []
    for iz in range(3):
        z = iz / 2
        for iy in range(3):
            y = iy / 2
            for ix in range(3):
                x = ix / 2
                coordinates.append((x, y, z))
                exact.append(
                    math.exp(-final_time)
                    * math.cos(math.pi * x)
                    * math.cos(math.pi * y)
                    * math.cos(math.pi * z)
                )
    numerical = list(exact)
    if perturb:
        numerical[1] += 0.1
    errors = [
        actual - expected
        for actual, expected in zip(numerical, exact, strict=True)
    ]
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(("x", "y", "z", "numerical", "exact", "error"))
    for coordinate, actual, expected, error in zip(
        coordinates, numerical, exact, errors, strict=True
    ):
        writer.writerow(
            tuple(format(value, ".17g") for value in coordinate)
            + (
                format(actual, ".17g"),
                format(expected, ".17g"),
                format(error, ".17g"),
            )
        )

    samples = [
        median_seconds * factor
        for factor in (1.0, 1.02, 0.98, 1.0, 1.0, 1.04, 0.96)
    ]
    if short:
        samples[0] = 0.0005
    return {
        "schema_version": 1,
        "kind": "ads-benchmark-case-result",
        "case_id": case.case_id,
        "status": "passed",
        "configuration": case.spec.to_dict(),
        "domain_result": {
            "problem": case.spec.problem,
            "scheme": case.spec.scheme,
            "exact_case": case.spec.exact_case,
            "solver_status": 0,
            "steps": case.spec.time.steps,
            "time_step": float(Fraction(case.spec.time.time_step)),
            "requested_final_time": final_time,
            "actual_final_time": final_time,
            "sample_points_per_axis": 3,
            "field_samples_written": True,
            "l2_error": math.sqrt(math.fsum(value * value for value in errors)),
            "linf_error": max(abs(value) for value in errors),
            "solution_l2_norm": math.sqrt(
                math.fsum(value * value for value in numerical)
            ),
            "field_checksum": _checksum(numerical),
            "physical_step_wall_seconds": samples[-1],
        },
        "timing": {
            "wall_seconds": 1.0,
            "metric": "physical_step_wall_seconds",
            "warmup_samples": [median_seconds, median_seconds],
            "measured_samples": samples,
            "warmup_process_wall_seconds": [1.1, 1.1],
            "measured_process_wall_seconds": [1.1] * 7,
            "minimum_reliable_seconds": 0.001,
            "reliable": min(samples) >= 0.001,
            "openmp_environment": {
                "OMP_NUM_THREADS": str(case.spec.openmp_threads),
                "OMP_DYNAMIC": "FALSE",
                "OMP_PROC_BIND": "close",
                "OMP_PLACES": "cores",
            },
        },
        "analysis_artifacts": {"field_samples_csv": stream.getvalue()},
    }


def _matrix() -> tuple[tuple[PlannedCase, ...], list[dict[str, object]]]:
    measurements = tuple(
        _case((2, 2, 2), grid, threads)
        for threads in (1, 2)
        for grid in ((1, 1, 1), (2, 1, 1))
    )
    reference = _case(
        (4, 2, 2), (1, 1, 1), 1, role="field-reference"
    )
    cases = measurements + (reference,)
    seconds = (8.0, 8.8, 4.0, 4.4, 12.0)
    return cases, [
        _result(case, elapsed)
        for case, elapsed in zip(cases, seconds, strict=True)
    ]


class WeakScalingAnalysisTests(unittest.TestCase):
    def test_complete_matrix_reports_weak_efficiency_without_speedup(self) -> None:
        cases, results = _matrix()
        report = WeakScalingAnalyzer().configure(cases).analyze(results)

        self.assertEqual(report.status, "passed")
        self.assertEqual(report.summary["case_count"], 5)
        self.assertEqual(report.summary["measurement_case_count"], 4)
        self.assertEqual(report.summary["field_reference_case_count"], 1)
        self.assertEqual(len(report.series), 2)
        self.assertEqual({series["openmp_threads"] for series in report.series}, {1, 2})
        omp_one = next(series for series in report.series if series["openmp_threads"] == 1)
        self.assertEqual(
            [point["weak_scaling_efficiency"] for point in omp_one["points"]],
            [1.0, 8.0 / 8.8],
        )
        self.assertTrue(all("speedup" not in point for point in omp_one["points"]))
        self.assertEqual(
            omp_one["points"][1]["field_reference_case_id"], cases[-1].case_id
        )

        with tempfile.TemporaryDirectory() as temporary:
            written = write_report(report, Path(temporary), plot=True)
            with written.csv_path.open(encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 4)
            self.assertIn("weak_scaling_efficiency", rows[0])
            self.assertNotIn("speedup", rows[0])
            self.assertTrue(all(row["measured_samples"].startswith("[") for row in rows))
            if written.plot_path is not None:
                self.assertEqual(written.plot_path.name, "weak-scaling.png")
                self.assertGreater(written.plot_path.stat().st_size, 0)
            else:
                self.assertIn("matplotlib is unavailable", written.plot_message or "")

    def test_unreliable_measurement_fails_series_and_is_not_given_efficiency(self) -> None:
        cases, results = _matrix()
        results[1] = _result(cases[1], 8.8, short=True)
        report = WeakScalingAnalyzer().configure(cases).analyze(results)

        self.assertEqual(report.status, "failed")
        self.assertEqual(report.summary["unreliable_measurement_count"], 1)
        point = next(
            point
            for series in report.series
            for point in series["points"]
            if point["case_id"] == cases[1].case_id
        )
        self.assertFalse(point["timing_valid"])
        self.assertIsNone(point["weak_scaling_efficiency"])

    def test_invalid_mpi_one_baseline_is_never_replaced_by_a_later_level(self) -> None:
        cases, results = _matrix()
        results[0] = _result(cases[0], 8.0, short=True)
        report = WeakScalingAnalyzer().configure(cases).analyze(results)

        self.assertEqual(report.status, "failed")
        omp_one = next(
            series for series in report.series if series["openmp_threads"] == 1
        )
        self.assertEqual(omp_one["status"], "failed")
        self.assertIsNone(omp_one["timing_reference_case_id"])
        self.assertTrue(
            all(
                point["weak_scaling_efficiency"] is None
                for point in omp_one["points"]
            )
        )
        self.assertTrue(
            any("MPI=1 timing baseline" in item for item in report.diagnostics)
        )

    def test_field_reference_is_a_gate_but_not_a_timing_point(self) -> None:
        cases, results = _matrix()
        results[-1] = _result(cases[-1], 12.0, short=True)
        report = WeakScalingAnalyzer().configure(cases).analyze(results)
        self.assertEqual(report.status, "passed")
        self.assertEqual(report.summary["unreliable_measurement_count"], 0)
        case_ids = {
            point["case_id"] for series in report.series for point in series["points"]
        }
        self.assertNotIn(cases[-1].case_id, case_ids)

    def test_field_mismatch_invalidates_only_affected_measurement(self) -> None:
        cases, results = _matrix()
        results[1] = _result(cases[1], 8.8, perturb=True)
        report = WeakScalingAnalyzer().configure(cases).analyze(results)
        self.assertEqual(report.status, "failed")
        self.assertEqual(report.summary["field_failure_count"], 1)
        self.assertEqual(report.summary["invalid_timing_count"], 1)

    def test_missing_same_mesh_reference_and_wrong_global_mesh_are_rejected(self) -> None:
        cases, results = _matrix()
        with self.assertRaisesRegex(AnalysisError, "same-mesh MPI=1, OMP=1"):
            WeakScalingAnalyzer().analyze(results[:-1])

        wrong_mesh = deepcopy(results)
        wrong_mesh[1]["configuration"]["mesh"]["elements"] = [99, 2, 2]
        with self.assertRaisesRegex(AnalysisError, "does not equal per-rank"):
            WeakScalingAnalyzer().analyze(wrong_mesh)


if __name__ == "__main__":
    unittest.main(verbosity=2)
