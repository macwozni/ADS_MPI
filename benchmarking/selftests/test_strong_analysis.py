from __future__ import annotations

import csv
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction
import io
import math
from pathlib import Path
import tempfile
import unittest

from ads_benchmark.analysis.base import AnalysisError
from ads_benchmark.analysis.io import write_report
from ads_benchmark.analysis.strong import StrongScalingAnalyzer
from ads_benchmark.framework.model import (
    CaseSpec,
    MeasurementSpec,
    MpiSpec,
    PlannedCase,
    SamplingSpec,
    TimeSpec,
)


def _case(
    grid: tuple[int, int, int],
    threads: int,
    *,
    proc_bind: str = "close",
    places: str = "cores",
) -> PlannedCase:
    spec = CaseSpec(
        family="strong",
        problem="igrm_l2",
        scheme="dg",
        exact_case="spatial-cosine",
        time=TimeSpec(final_time="0.1", time_step="0.025", steps=4),
        mesh=(2, 2, 2),
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
        openmp_proc_bind=proc_bind,
        openmp_places=places,
    )
    return PlannedCase(
        case_id=(
            f"igrm_l2-dg-{grid[0]}{grid[1]}{grid[2]}-t{threads}-"
            f"{proc_bind}-{places}"
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
    reliable = min(samples) >= 0.001
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
            "reliable": reliable,
            "openmp_environment": {
                "OMP_NUM_THREADS": str(case.spec.openmp_threads),
                "OMP_DYNAMIC": "FALSE",
                "OMP_PROC_BIND": str(case.spec.openmp_proc_bind),
                "OMP_PLACES": str(case.spec.openmp_places),
            },
        },
        "analysis_artifacts": {"field_samples_csv": stream.getvalue()},
    }


def _matrix() -> tuple[tuple[PlannedCase, ...], list[dict[str, object]]]:
    cases = tuple(
        _case(grid, threads)
        for grid in ((1, 1, 1), (2, 1, 1))
        for threads in (1, 2)
    )
    medians = (8.0, 4.0, 4.0, 2.0)
    return cases, [
        _result(case, seconds)
        for case, seconds in zip(cases, medians, strict=True)
    ]


class StrongScalingAnalysisTests(unittest.TestCase):
    def test_complete_matrix_reports_samples_statistics_speedup_and_efficiency(self) -> None:
        cases, results = _matrix()
        report = StrongScalingAnalyzer().configure(cases).analyze(results)
        self.assertEqual(report.status, "passed")
        self.assertEqual(report.summary["case_count"], 4)
        self.assertEqual(report.summary["invalid_timing_count"], 0)
        self.assertEqual(len(report.series), 2)
        points = [point for series in report.series for point in series["points"]]
        self.assertEqual({point["resources"] for point in points}, {1, 2, 4})
        self.assertTrue(all(len(point["measured_samples"]) == 7 for point in points))
        hybrid = next(point for point in points if point["resources"] == 4)
        self.assertEqual(hybrid["statistics"]["median_seconds"], 2.0)
        self.assertEqual(hybrid["speedup"], 4.0)
        self.assertEqual(hybrid["efficiency"], 1.0)

        with tempfile.TemporaryDirectory() as temporary:
            written = write_report(report, Path(temporary), plot=True)
            with written.csv_path.open(encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 4)
            self.assertEqual(
                {row["mpi_grid"] for row in rows}, {"[1,1,1]", "[2,1,1]"}
            )
            self.assertTrue(all(row["measured_samples"].startswith("[") for row in rows))
            self.assertTrue(
                all(
                    row["measured_process_wall_seconds"].startswith("[")
                    for row in rows
                )
            )
            if written.plot_path is not None:
                self.assertEqual(written.plot_path.name, "strong-scaling.png")
                self.assertGreater(written.plot_path.stat().st_size, 0)
            else:
                self.assertIn("matplotlib is unavailable", written.plot_message or "")

    def test_short_region_is_reported_and_excluded_from_scaling(self) -> None:
        cases, results = _matrix()
        results[-1] = _result(cases[-1], 2.0, short=True)
        report = StrongScalingAnalyzer().configure(cases).analyze(results)
        self.assertEqual(report.status, "failed")
        self.assertEqual(report.summary["unreliable_measurement_count"], 1)
        self.assertEqual(report.summary["invalid_timing_count"], 1)
        point = next(
            point
            for series in report.series
            for point in series["points"]
            if point["case_id"] == cases[-1].case_id
        )
        self.assertFalse(point["timing_valid"])
        self.assertIsNone(point["speedup"])

    def test_smallest_valid_resource_becomes_timing_baseline(self) -> None:
        cases, results = _matrix()
        results[0] = _result(cases[0], 8.0, short=True)
        report = StrongScalingAnalyzer().configure(cases).analyze(results)

        self.assertEqual(report.status, "failed")
        self.assertTrue(
            all(
                series["timing_reference_case_id"] == cases[1].case_id
                for series in report.series
            )
        )
        fallback = next(
            point
            for series in report.series
            for point in series["points"]
            if point["case_id"] == cases[1].case_id
        )
        self.assertEqual(fallback["resources"], 2)
        self.assertEqual(fallback["speedup"], 1.0)
        self.assertEqual(fallback["efficiency"], 1.0)

    def test_openmp_binding_policies_form_distinct_physics_groups(self) -> None:
        cases = tuple(
            _case((1, 1, 1), threads, proc_bind=policy)
            for policy in ("close", "spread")
            for threads in (1, 2)
        )
        report = StrongScalingAnalyzer().configure(cases).analyze(
            [_result(case, 8.0 / case.spec.openmp_threads) for case in cases]
        )

        self.assertEqual(report.status, "passed")
        self.assertEqual(len(report.series), 2)
        policies = {
            series["controlled_configuration"]["openmp"]["proc_bind"]
            for series in report.series
        }
        self.assertEqual(policies, {"close", "spread"})
        self.assertTrue(
            all(
                "threads" not in series["controlled_configuration"]["openmp"]
                for series in report.series
            )
        )

    def test_corrupt_repetition_contract_and_result_schema_are_rejected(self) -> None:
        cases, results = _matrix()
        corruptions = []

        missing_sample = deepcopy(results)
        missing_sample[0]["timing"]["measured_samples"].pop()
        corruptions.append((missing_sample, "measured count"))

        wrong_threshold = deepcopy(results)
        wrong_threshold[0]["timing"]["minimum_reliable_seconds"] = 0.002
        corruptions.append((wrong_threshold, "threshold differs"))

        wrong_environment = deepcopy(results)
        wrong_environment[0]["timing"]["openmp_environment"][
            "OMP_PROC_BIND"
        ] = "spread"
        corruptions.append((wrong_environment, "environment differs"))

        wrong_final_sample = deepcopy(results)
        wrong_final_sample[0]["domain_result"][
            "physical_step_wall_seconds"
        ] = 123.0
        corruptions.append((wrong_final_sample, "last measured sample"))

        wrong_schema = deepcopy(results)
        wrong_schema[0]["schema_version"] = 2
        corruptions.append((wrong_schema, "schema_version"))

        extra_root_key = deepcopy(results)
        extra_root_key[0]["unexpected"] = True
        corruptions.append((extra_root_key, "invalid keys"))

        analyzer = StrongScalingAnalyzer().configure(cases)
        for documents, message in corruptions:
            with self.subTest(message=message), self.assertRaisesRegex(
                AnalysisError, message
            ):
                analyzer.analyze(documents)

    def test_plot_is_explicitly_omitted_for_ambiguous_series_count(self) -> None:
        cases, results = _matrix()
        report = StrongScalingAnalyzer().configure(cases).analyze(results)
        crowded = replace(report, series=report.series * 7)
        with tempfile.TemporaryDirectory() as temporary:
            written = write_report(crowded, Path(temporary), plot=True)
        self.assertIsNone(written.plot_path)
        self.assertIn("filter to at most 12", written.plot_message or "")

    def test_field_mismatch_invalidates_all_samples_for_that_configuration(self) -> None:
        cases, results = _matrix()
        results[-1] = _result(cases[-1], 2.0, perturb=True)
        report = StrongScalingAnalyzer().configure(cases).analyze(results)
        self.assertEqual(report.status, "failed")
        self.assertEqual(report.summary["field_failure_count"], 1)
        self.assertEqual(report.summary["invalid_timing_count"], 1)

    def test_missing_serial_field_reference_is_an_explicit_error(self) -> None:
        cases, results = _matrix()
        filtered = [
            result
            for case, result in zip(cases, results, strict=True)
            if not (
                case.spec.mpi.ranks == 1 and case.spec.openmp_threads == 1
            )
        ]
        with self.assertRaisesRegex(AnalysisError, "exactly one MPI=1, OMP=1"):
            StrongScalingAnalyzer().analyze(filtered)


if __name__ == "__main__":
    unittest.main(verbosity=2)
