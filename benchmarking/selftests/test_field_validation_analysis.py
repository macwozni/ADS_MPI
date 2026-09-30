from __future__ import annotations

import csv
from fractions import Fraction
import io
import math
from pathlib import Path
import tempfile
import unittest

from ads_benchmark.analysis.base import AnalysisError
from ads_benchmark.analysis.io import write_report
from ads_benchmark.analysis.validation import FieldValidationAnalyzer
from ads_benchmark.framework.model import (
    CaseSpec,
    MeasurementSpec,
    MpiSpec,
    PlannedCase,
    SamplingSpec,
    TimeSpec,
)


def _case(scheme: str, steps: int, grid: tuple[int, int, int], omp: int) -> PlannedCase:
    spec = CaseSpec(
        family="validation",
        problem="igrm_l2",
        scheme=scheme,
        exact_case="spatial-cosine",
        time=TimeSpec(
            final_time="0.01",
            time_step=str(Fraction(1, 100 * steps)),
            steps=steps,
        ),
        mesh=(2, 2, 2),
        test_degree=(3, 3, 3),
        trial_degree=(2, 2, 2),
        mpi=MpiSpec(ranks=math.prod(grid), process_grid=grid),
        openmp_threads=omp,
        sampling=SamplingSpec(points_per_axis=3, write_samples=True),
        measurement=MeasurementSpec(warmups=0, samples=1, timeout_seconds="60"),
        build_profile="release",
        launcher="mpi",
    )
    return PlannedCase(
        case_id=f"igrm_l2-{scheme}-n{steps}-{grid[0]}{grid[1]}{grid[2]}-t{omp}",
        spec=spec,
    )


def _checksum(values: list[float]) -> float:
    sample_sum = 0.0
    sample_compensation = 0.0
    weighted_sum = 0.0
    weighted_compensation = 0.0
    for index, value in enumerate(values):
        corrected = value - sample_compensation
        updated = sample_sum + corrected
        sample_compensation = (updated - sample_sum) - corrected
        sample_sum = updated
        corrected = (index + 1) * value - weighted_compensation
        updated = weighted_sum + corrected
        weighted_compensation = (updated - weighted_sum) - corrected
        weighted_sum = updated
    return sample_sum + weighted_sum / (len(values) + 1)


def _result(
    case: PlannedCase,
    *,
    perturb_parallel: bool = False,
    inaccurate_reference: bool = False,
    field_bias: float = 0.0,
    be_error_amplitude: float | None = None,
) -> dict[str, object]:
    final_time = float(Fraction(case.spec.time.final_time))
    coordinates: list[tuple[float, float, float]] = []
    exact: list[float] = []
    numerical: list[float] = []
    for iz in range(3):
        z = iz / 2
        for iy in range(3):
            y = iy / 2
            for ix in range(3):
                x = ix / 2
                coordinates.append((x, y, z))
                expected = (
                    math.exp(-final_time)
                    * math.cos(math.pi * x)
                    * math.cos(math.pi * y)
                    * math.cos(math.pi * z)
                )
                exact.append(expected)
                spatial_error = 0.005 * (x + y + z)
                scheme_amplitude = (
                    0.004 * (128 / case.spec.time.steps)
                    if be_error_amplitude is None
                    else be_error_amplitude
                )
                scheme_error = (
                    scheme_amplitude * (1.0 + x)
                    if case.spec.scheme == "be"
                    else 0.0
                )
                numerical.append(
                    expected + spatial_error + scheme_error + field_bias
                )

    if inaccurate_reference:
        numerical = [value + 0.2 for value in numerical]
    if perturb_parallel:
        numerical[1], numerical[2] = numerical[2], numerical[1]
    errors = [value - expected for value, expected in zip(numerical, exact, strict=True)]
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(("x", "y", "z", "numerical", "exact", "error"))
    for coordinate, value, expected, error in zip(
        coordinates, numerical, exact, errors, strict=True
    ):
        writer.writerow(
            tuple(format(component, ".17g") for component in coordinate)
            + (
                format(value, ".17g"),
                format(expected, ".17g"),
                format(error, ".17g"),
            )
        )

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
        },
        "timing": {"wall_seconds": 1.25},
        "analysis_artifacts": {"field_samples_csv": stream.getvalue()},
    }


def _matrix() -> tuple[tuple[PlannedCase, ...], list[dict[str, object]]]:
    cases = tuple(
        _case(scheme, steps, grid, 1)
        for scheme in ("dg", "be")
        for steps in (128, 256)
        for grid in ((1, 1, 1), (1, 2, 1))
    )
    return cases, [_result(case) for case in cases]


class FieldValidationAnalysisTests(unittest.TestCase):
    def test_complete_matrix_passes_and_writes_readable_rows(self) -> None:
        cases, results = _matrix()
        report = FieldValidationAnalyzer().configure(cases).analyze(results)
        self.assertEqual(report.status, "passed")
        self.assertEqual(report.summary["case_count"], 8)
        self.assertEqual(report.summary["parallel_comparison_count"], 4)
        self.assertEqual(report.summary["invalid_timing_count"], 0)
        self.assertEqual(len(report.series), 1)
        self.assertEqual(
            report.series[0]["scheme_agreement"][0]["status"], "passed"
        )

        with tempfile.TemporaryDirectory() as temporary_directory:
            written = write_report(
                report, Path(temporary_directory), plot=True
            )
            with written.csv_path.open(encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(
                {row["row_kind"] for row in rows},
                {"analytic-level", "parallel-variant", "scheme-pair"},
            )
            parallel = next(row for row in rows if row["row_kind"] == "parallel-variant")
            self.assertEqual(parallel["mpi_grid"], "[1,2,1]")
            self.assertTrue(parallel["worst_coordinates"].startswith("["))
            final_pair = next(
                row
                for row in rows
                if row["row_kind"] == "scheme-pair"
                and row["qualification"] == "final-accuracy"
            )
            self.assertEqual(final_pair["agreement_status"], "passed")
            self.assertEqual(final_pair["trend_passed"], "True")
            self.assertNotEqual(final_pair["coarse_l2_difference"], "")
            self.assertNotEqual(final_pair["final_linf_difference"], "")
            self.assertIsNone(written.plot_path)
            self.assertIn("no convergence plot", written.plot_message or "")

    def test_parallel_local_permutation_invalidates_only_eligible_timing(self) -> None:
        cases, results = _matrix()
        target = next(
            index
            for index, case in enumerate(cases)
            if case.spec.scheme == "dg"
            and case.spec.time.steps == 256
            and case.spec.mpi.ranks == 2
        )
        results[target] = _result(cases[target], perturb_parallel=True)
        report = FieldValidationAnalyzer().configure(cases).analyze(results)
        self.assertEqual(report.status, "failed")
        self.assertEqual(report.summary["invalid_timing_count"], 1)
        self.assertTrue(any("grid=(1, 2, 1)" in item for item in report.diagnostics))
        variants = [
            variant
            for level in report.series[0]["time_levels"]
            for scheme in level["schemes"]
            for variant in scheme["variants"]
        ]
        invalid = [item for item in variants if not item["timing_valid"]]
        self.assertEqual([item["case_id"] for item in invalid], [cases[target].case_id])

    def test_missing_serial_reference_is_an_explicit_analysis_error(self) -> None:
        cases, results = _matrix()
        filtered = [
            result
            for case, result in zip(cases, results, strict=True)
            if not (
                case.spec.scheme == "dg"
                and case.spec.time.steps == 128
                and case.spec.mpi.ranks == 1
            )
        ]
        with self.assertRaisesRegex(AnalysisError, "exactly one MPI=1, OMP=1"):
            FieldValidationAnalyzer().analyze(filtered)

    def test_final_analytical_accuracy_is_a_gate(self) -> None:
        cases, results = _matrix()
        target = next(
            index
            for index, case in enumerate(cases)
            if case.spec.scheme == "dg"
            and case.spec.time.steps == 256
            and case.spec.mpi.ranks == 1
        )
        results[target] = _result(cases[target], inaccurate_reference=True)
        report = FieldValidationAnalyzer().configure(cases).analyze(results)
        self.assertEqual(report.status, "failed")
        self.assertEqual(report.summary["analytic_final_failure_count"], 1)
        self.assertTrue(any("analytical accuracy" in item for item in report.diagnostics))

    def test_identically_wrong_final_layouts_are_not_timing_eligible(self) -> None:
        cases, _ = _matrix()
        cases = tuple(case for case in cases if case.spec.scheme == "dg")
        results = [
            _result(
                case,
                field_bias=0.06 if case.spec.time.steps == 256 else 0.0,
            )
            for case in cases
        ]
        report = FieldValidationAnalyzer().configure(cases).analyze(results)
        self.assertEqual(report.status, "failed")
        self.assertEqual(report.summary["analytic_final_failure_count"], 1)
        self.assertEqual(report.summary["invalid_timing_count"], 2)
        final_scheme = report.series[0]["time_levels"][-1]["schemes"][0]
        self.assertFalse(final_scheme["reference"]["timing_valid"])
        self.assertTrue(final_scheme["variants"][0]["comparison"]["passed"])
        self.assertFalse(final_scheme["variants"][0]["timing_valid"])

    def test_flat_cross_scheme_difference_fails_convergence_trend(self) -> None:
        cases, results = _matrix()
        for index, case in enumerate(cases):
            if case.spec.scheme == "be" and case.spec.time.steps == 256:
                results[index] = _result(case, be_error_amplitude=0.004)
        report = FieldValidationAnalyzer().configure(cases).analyze(results)
        agreement = report.series[0]["scheme_agreement"][0]
        self.assertEqual(report.status, "failed")
        self.assertTrue(agreement["final_pointwise_accuracy_passed"])
        self.assertFalse(agreement["trend_passed"])
        self.assertEqual(agreement["status"], "failed")

    def test_unit_amplitude_scheme_error_above_one_percent_is_rejected(self) -> None:
        cases, results = _matrix()
        for index, case in enumerate(cases):
            if case.spec.scheme == "be" and case.spec.time.steps == 256:
                # At x=1 the BE/DG difference is exactly 0.015 against a
                # unit-amplitude field, which must not pass a declared 1% gate.
                results[index] = _result(case, be_error_amplitude=0.0075)
        report = FieldValidationAnalyzer().configure(cases).analyze(results)
        agreement = report.series[0]["scheme_agreement"][0]
        self.assertEqual(report.status, "failed")
        self.assertFalse(agreement["final_pointwise_accuracy_passed"])
        self.assertEqual(report.summary["scheme_relative_tolerance"], 0.0)

    def test_frozen_plan_mismatch_is_rejected_before_partial_claim(self) -> None:
        cases, results = _matrix()
        with self.assertRaisesRegex(AnalysisError, "frozen validation plan"):
            FieldValidationAnalyzer().configure(cases).analyze(results[:-1])


if __name__ == "__main__":
    unittest.main()
