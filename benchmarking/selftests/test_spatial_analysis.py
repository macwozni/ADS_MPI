from __future__ import annotations

import csv
from copy import deepcopy
from fractions import Fraction
import io
import json
import math
from pathlib import Path
import tempfile
import unittest

from ads_benchmark.analysis.base import AnalysisError
from ads_benchmark.analysis.io import write_report
from ads_benchmark.analysis.spatial import (
    HConvergenceAnalyzer,
    PConvergenceAnalyzer,
)
from ads_benchmark.framework.model import (
    CaseSpec,
    MeasurementSpec,
    MpiSpec,
    PlannedCase,
    SamplingSpec,
    TimeSpec,
)


SYNTHETIC_SAMPLE_POINTS = 5


def _field_samples_csv(
    l2_error: float, linf_error: float, *, sign: float = 1.0
) -> str:
    """Build a regular field whose trapezoidal L2/Linf are prescribed."""

    if not (0.0 < l2_error <= linf_error):
        raise ValueError("synthetic errors must satisfy 0 < L2 <= Linf")
    points = SYNTHETIC_SAMPLE_POINTS
    spacing = 1.0 / (points - 1)
    peak_weight = spacing**3
    background_squared = (
        l2_error * l2_error - peak_weight * linf_error * linf_error
    ) / (1.0 - peak_weight)
    if background_squared < 0.0:
        raise ValueError("synthetic Linf/L2 ratio is too large for the grid")
    background = math.sqrt(background_squared)
    center = points // 2
    rows = [["x", "y", "z", "numerical", "exact", "error"]]
    for iz in range(points):
        z = iz * spacing
        for iy in range(points):
            y = iy * spacing
            for ix in range(points):
                x = ix * spacing
                exact = (
                    math.exp(-0.1)
                    * math.cos(math.pi * x)
                    * math.cos(math.pi * y)
                    * math.cos(math.pi * z)
                )
                magnitude = (
                    linf_error
                    if (ix, iy, iz) == (center, center, center)
                    else background
                )
                error = sign * magnitude
                rows.append([x, y, z, exact + error, exact, error])
    stream = io.StringIO(newline="")
    csv.writer(stream, lineterminator="\n").writerows(rows)
    return stream.getvalue()


def _attach_field_samples(
    result: dict[str, object], *, sign: float = 1.0
) -> None:
    domain = result["domain_result"]
    text = _field_samples_csv(
        domain["l2_error"], domain["linf_error"], sign=sign
    )
    reader = csv.reader(io.StringIO(text, newline=""))
    next(reader)
    sample_sum = 0.0
    sample_compensation = 0.0
    weighted_sum = 0.0
    weighted_compensation = 0.0
    count = 0
    for count, row in enumerate(reader, start=1):
        numerical = float(row[3])
        corrected = numerical - sample_compensation
        update = sample_sum + corrected
        sample_compensation = (update - sample_sum) - corrected
        sample_sum = update
        corrected = count * numerical - weighted_compensation
        update = weighted_sum + corrected
        weighted_compensation = (update - weighted_sum) - corrected
        weighted_sum = update
    domain["field_checksum"] = sample_sum + weighted_sum / (count + 1)
    result["analysis_artifacts"] = {
        "field_samples_csv": text
    }


def _result(
    *,
    family: str,
    level_label: str,
    mesh: tuple[int, int, int],
    trial: tuple[int, int, int],
    test: tuple[int, int, int],
    steps: int,
    l2_error: float,
    linf_error: float,
    problem: str = "synthetic_problem",
    scheme: str = "dg",
) -> dict[str, object]:
    final_time = Fraction(1, 10)
    time_step = final_time / steps
    configuration = {
        "family": family,
        "problem": problem,
        "scheme": scheme,
        "exact_case": "spatial-cosine",
        "time": {
            "final_time": "0.1",
            "time_step": str(time_step),
            "steps": steps,
        },
        "mesh": {"elements": list(mesh)},
        "spaces": {
            "test_degree": list(test),
            "trial_degree": list(trial),
        },
        "mpi": {"ranks": 1, "process_grid": [1, 1, 1]},
        "openmp": {"threads": 1},
        "sampling": {
            "points_per_axis": SYNTHETIC_SAMPLE_POINTS,
            "write_samples": True,
        },
        "measurement": {
            "warmups": 0,
            "samples": 1,
            "timeout_seconds": "60",
        },
        "build": {"profile": "debug"},
        "launcher": "direct",
    }
    case_id = f"{family}-{level_label}-{steps}-{scheme}-{trial}-{test}"
    result = {
        "schema_version": 1,
        "kind": "ads-benchmark-case-result",
        "case_id": case_id,
        "status": "passed",
        "configuration": configuration,
        "timing": {"wall_seconds": 0.01},
        "domain_result": {
            "schema_version": 1,
            "kind": "ads-manufactured-transient-result",
            "problem": problem,
            "scheme": scheme,
            "exact_case": "spatial-cosine",
            "requested_final_time": 0.1,
            "actual_final_time": 0.1,
            "time_step": float(time_step),
            "steps": steps,
            "l2_error": l2_error,
            "linf_error": linf_error,
            "solution_l2_norm": 0.3,
            "solver_status": 0,
            "sample_points_per_axis": SYNTHETIC_SAMPLE_POINTS,
            "field_samples_written": True,
        },
    }
    _attach_field_samples(result)
    return result


def _paired_results(
    *,
    family: str,
    level_label: str,
    mesh: tuple[int, int, int],
    trial: tuple[int, int, int],
    test: tuple[int, int, int],
    fine_l2: float,
    fine_linf: float,
    coarse_factor: float = 1.04,
    fine_steps: int = 20,
    problem: str = "synthetic_problem",
    scheme: str = "dg",
) -> list[dict[str, object]]:
    return [
        _result(
            family=family,
            level_label=level_label,
            mesh=mesh,
            trial=trial,
            test=test,
            steps=fine_steps // 2,
            l2_error=coarse_factor * fine_l2,
            linf_error=coarse_factor * fine_linf,
            problem=problem,
            scheme=scheme,
        ),
        _result(
            family=family,
            level_label=level_label,
            mesh=mesh,
            trial=trial,
            test=test,
            steps=fine_steps,
            l2_error=fine_l2,
            linf_error=fine_linf,
            problem=problem,
            scheme=scheme,
        ),
    ]


def h_results(
    meshes: tuple[int, ...] = (2, 4, 8, 16),
    *,
    order: float = 4.0,
    scheme: str = "dg",
) -> list[dict[str, object]]:
    results: list[dict[str, object]] = []
    for elements in meshes:
        error = (1.0 / elements) ** order
        results.extend(
            _paired_results(
                family="h",
                level_label=f"n{elements}",
                mesh=(elements, elements, elements),
                trial=(3, 3, 3),
                test=(4, 4, 4),
                fine_l2=error,
                fine_linf=2.5 * error,
                scheme=scheme,
            )
        )
    return results


def p_results(
    degrees: tuple[int, ...] = (1, 2, 3, 4),
    *,
    enrichment: int = 1,
    errors: tuple[float, ...] | None = None,
) -> list[dict[str, object]]:
    if errors is None:
        errors = tuple(0.2**degree for degree in degrees)
    results: list[dict[str, object]] = []
    for degree, error in zip(degrees, errors, strict=True):
        results.extend(
            _paired_results(
                family="p",
                level_label=f"p{degree}-e{enrichment}",
                mesh=(4, 4, 4),
                trial=(degree, degree, degree),
                test=(
                    degree + enrichment,
                    degree + enrichment,
                    degree + enrichment,
                ),
                fine_l2=error,
                fine_linf=3.0 * error,
            )
        )
    return results


def _planned_case(result: dict[str, object]) -> PlannedCase:
    configuration = result["configuration"]
    time = configuration["time"]
    mesh = tuple(configuration["mesh"]["elements"])
    test = tuple(configuration["spaces"]["test_degree"])
    trial = tuple(configuration["spaces"]["trial_degree"])
    spec = CaseSpec(
        family=configuration["family"],
        problem=configuration["problem"],
        scheme=configuration["scheme"],
        exact_case=configuration["exact_case"],
        time=TimeSpec(
            final_time=time["final_time"],
            time_step=time["time_step"],
            steps=time["steps"],
        ),
        mesh=mesh,
        test_degree=test,
        trial_degree=trial,
        mpi=MpiSpec(ranks=1, process_grid=(1, 1, 1)),
        openmp_threads=1,
        sampling=SamplingSpec(
            points_per_axis=SYNTHETIC_SAMPLE_POINTS, write_samples=True
        ),
        measurement=MeasurementSpec(
            warmups=0, samples=1, timeout_seconds="60"
        ),
        build_profile="debug",
        launcher="direct",
    )
    return PlannedCase(case_id=result["case_id"], spec=spec)


class HConvergenceAnalysisTests(unittest.TestCase):
    def setUp(self) -> None:
        self.analyzer = HConvergenceAnalyzer()

    def test_h_orders_use_fine_member_of_each_temporal_pair(self) -> None:
        report = self.analyzer.analyze(h_results(), source_run="synthetic-h")
        self.assertEqual(report.status, "passed")
        self.assertEqual(report.source_run, "synthetic-h")
        self.assertEqual(report.summary["spatial_point_count"], 4)
        series = report.series[0]
        self.assertEqual(series["acceptance_mode"], "positive-h-regression")
        self.assertEqual(series["sequence_kind"], "h-refinement")
        self.assertEqual([point["mesh"] for point in series["points"]], [
            [2, 2, 2],
            [4, 4, 4],
            [8, 8, 8],
            [16, 16, 16],
        ])
        for metric_name in ("l2", "linf"):
            metric = series["metrics"][metric_name]
            self.assertEqual(metric["status"], "passed")
            self.assertAlmostEqual(metric["regression"]["order"], 4.0, places=12)
            self.assertEqual(len(metric["local_orders"]), 3)
            for estimate in metric["local_orders"]:
                self.assertAlmostEqual(estimate["order"], 4.0, places=12)
                self.assertTrue(estimate["available"])
        first = series["points"][0]
        self.assertIn("coarse_case_id", first)
        self.assertIn("fine_case_id", first)
        self.assertAlmostEqual(
            first["temporal_pair"]["l2"]["temporal_error_fraction"],
            0.04,
            places=12,
        )

    def test_time_dominated_point_is_audited_and_not_used_for_orders(self) -> None:
        results = h_results()
        # mesh=8 coarse member: temporal change is now 20% of the fine error.
        target = next(
            result
            for result in results
            if result["configuration"]["mesh"]["elements"] == [8, 8, 8]
            and result["configuration"]["time"]["steps"] == 10
        )
        fine_error = next(
            result["domain_result"]["l2_error"]
            for result in results
            if result["configuration"]["mesh"]["elements"] == [8, 8, 8]
            and result["configuration"]["time"]["steps"] == 20
        )
        target["domain_result"]["l2_error"] = 1.2 * fine_error
        _attach_field_samples(target)

        report = self.analyzer.analyze(results)
        self.assertEqual(report.status, "failed")
        metric = report.series[0]["metrics"]["l2"]
        self.assertEqual(metric["status"], "unreliable")
        self.assertEqual(len(metric["time_dominated_points"]), 1)
        self.assertIsNone(metric["local_orders"][1]["order"])
        self.assertIsNone(metric["local_orders"][2]["order"])
        self.assertIn("time-dominated", metric["local_orders"][1]["unavailable_reason"])

    def test_nonfinite_derived_indicator_is_reported_as_unreliable(self) -> None:
        results = h_results()
        coarse = next(
            result
            for result in results
            if result["configuration"]["mesh"]["elements"] == [2, 2, 2]
            and result["configuration"]["time"]["steps"] == 10
        )
        fine = next(
            result
            for result in results
            if result["configuration"]["mesh"]["elements"] == [2, 2, 2]
            and result["configuration"]["time"]["steps"] == 20
        )

        coarse_domain = coarse["domain_result"]
        coarse_domain["l2_error"] = 1.0e150
        coarse_domain["linf_error"] = 1.0e150
        _attach_field_samples(coarse)
        fine_domain = fine["domain_result"]
        fine_domain["l2_error"] = 1.0e-200
        fine_domain["linf_error"] = 1.0e-200
        _attach_field_samples(fine)

        report = self.analyzer.analyze(results)
        self.assertEqual(report.status, "failed")
        point = report.series[0]["points"][0]
        for metric_name in ("l2", "linf"):
            audit = point["temporal_pair"][metric_name]
            self.assertEqual(audit["temporal_indicator_status"], "nonfinite")
            self.assertIsNone(audit["temporal_error_fraction"])
            self.assertFalse(audit["reliable_for_spatial_analysis"])

        with tempfile.TemporaryDirectory(prefix="ads-nonfinite-derived-") as temporary:
            written = write_report(report, Path(temporary))
            document = json.loads(written.json_path.read_text(encoding="utf-8"))
        written_audit = document["series"][0]["points"][0]["temporal_pair"]["l2"]
        self.assertEqual(written_audit["temporal_indicator_status"], "nonfinite")
        self.assertIsNone(written_audit["temporal_error_fraction"])

    def test_equal_error_norms_cannot_hide_opposite_temporal_fields(self) -> None:
        results = h_results()
        coarse = next(
            result
            for result in results
            if result["configuration"]["mesh"]["elements"] == [2, 2, 2]
            and result["configuration"]["time"]["steps"] == 10
        )
        fine = next(
            result
            for result in results
            if result["configuration"]["mesh"]["elements"] == [2, 2, 2]
            and result["configuration"]["time"]["steps"] == 20
        )
        coarse["domain_result"]["l2_error"] = fine["domain_result"]["l2_error"]
        coarse["domain_result"]["linf_error"] = fine["domain_result"][
            "linf_error"
        ]
        _attach_field_samples(coarse, sign=-1.0)

        report = self.analyzer.analyze(results)
        self.assertEqual(report.status, "failed")
        point = report.series[0]["points"][0]
        for metric_name in ("l2", "linf"):
            audit = point["temporal_pair"][metric_name]
            self.assertGreater(audit["temporal_error_fraction"], 1.9)
            self.assertGreater(audit["temporal_field_difference"], 0.0)
            self.assertNotIn("temporal_error_change", audit)

    def test_field_sample_artifact_is_required_and_strictly_bound(self) -> None:
        missing = h_results()
        del missing[0]["analysis_artifacts"]
        with self.assertRaisesRegex(AnalysisError, "analysis_artifacts must be an object"):
            self.analyzer.analyze(missing)

        bad_order = h_results()
        artifact = bad_order[0]["analysis_artifacts"]["field_samples_csv"]
        lines = artifact.splitlines()
        lines[1], lines[2] = lines[2], lines[1]
        bad_order[0]["analysis_artifacts"]["field_samples_csv"] = (
            "\n".join(lines) + "\n"
        )
        with self.assertRaisesRegex(AnalysisError, "regular-grid order"):
            self.analyzer.analyze(bad_order)

        tampered = h_results()
        text = tampered[0]["analysis_artifacts"]["field_samples_csv"]
        rows = list(csv.reader(io.StringIO(text, newline="")))
        delta = 1.0e-3
        rows[1][3] = str(float(rows[1][3]) + delta)
        rows[1][5] = str(float(rows[1][5]) + delta)
        stream = io.StringIO(newline="")
        csv.writer(stream, lineterminator="\n").writerows(rows)
        tampered[0]["analysis_artifacts"]["field_samples_csv"] = stream.getvalue()
        with self.assertRaisesRegex(AnalysisError, "sample checksum differs"):
            self.analyzer.analyze(tampered)

    def test_h_plateau_is_detected_and_excluded_from_regression(self) -> None:
        meshes = (2, 4, 8, 16, 32)
        results = h_results(meshes, order=2.0)
        plateau_error = (1.0 / 8.0) ** 2
        for result in results:
            elements = result["configuration"]["mesh"]["elements"][0]
            if elements >= 16:
                is_coarse = result["configuration"]["time"]["steps"] == 10
                result["domain_result"]["l2_error"] = plateau_error * (
                    1.04 if is_coarse else 1.0
                )
                result["domain_result"]["linf_error"] = 2.5 * plateau_error * (
                    1.04 if is_coarse else 1.0
                )
                _attach_field_samples(result)
        report = self.analyzer.analyze(results)
        self.assertEqual(report.status, "passed")
        for metric in report.series[0]["metrics"].values():
            self.assertTrue(metric["plateau_detected"])
            self.assertEqual(metric["plateau_start_mesh"], [16, 16, 16])
            self.assertEqual(metric["regression"]["point_count"], 3)
            self.assertAlmostEqual(metric["regression"]["order"], 2.0, places=12)

    def test_nonmonotonic_h_error_outside_plateau_fails(self) -> None:
        results = h_results()
        for result in results:
            if (
                result["configuration"]["mesh"]["elements"] == [8, 8, 8]
                and result["configuration"]["time"]["steps"] == 20
            ):
                result["domain_result"]["l2_error"] = 1.0
                result["domain_result"]["linf_error"] = 2.5
                _attach_field_samples(result)
            if (
                result["configuration"]["mesh"]["elements"] == [8, 8, 8]
                and result["configuration"]["time"]["steps"] == 10
            ):
                result["domain_result"]["l2_error"] = 1.04
                result["domain_result"]["linf_error"] = 2.6
                _attach_field_samples(result)
        report = self.analyzer.analyze(results)
        self.assertEqual(report.status, "failed")
        self.assertTrue(
            report.series[0]["metrics"]["l2"]["fatal_nonmonotonic_anomalies"]
        )

    def test_temporal_pair_must_be_exact_dt_and_dt_over_two(self) -> None:
        results = h_results()
        item = results[1]
        item["configuration"]["time"] = {
            "final_time": "0.1",
            "time_step": "1/30",
            "steps": 3,
        }
        item["domain_result"]["time_step"] = 1.0 / 30.0
        item["domain_result"]["steps"] = 3
        with self.assertRaisesRegex(AnalysisError, "must use dt and dt/2"):
            self.analyzer.analyze(results)

    def test_missing_temporal_partner_is_rejected(self) -> None:
        with self.assertRaisesRegex(AnalysisError, "exactly two temporal runs"):
            self.analyzer.analyze(h_results()[:-1])

    def test_only_identical_controls_share_a_series(self) -> None:
        first = h_results(scheme="dg")
        second = h_results(scheme="pr")
        report = self.analyzer.analyze(first + second)
        self.assertEqual(report.summary["series_count"], 2)
        self.assertEqual({series["scheme"] for series in report.series}, {"dg", "pr"})

    def test_degree_change_is_not_hidden_inside_an_h_series(self) -> None:
        results = h_results()
        for result in results:
            if result["configuration"]["mesh"]["elements"] == [16, 16, 16]:
                result["configuration"]["spaces"] = {
                    "trial_degree": [2, 2, 2],
                    "test_degree": [3, 3, 3],
                }
        with self.assertRaisesRegex(AnalysisError, "h convergence requires at least"):
            self.analyzer.analyze(results)

    def test_non_spatial_manufactured_case_is_rejected(self) -> None:
        results = h_results()
        results[0]["configuration"]["exact_case"] = "temporal-polynomial"
        results[0]["domain_result"]["exact_case"] = "temporal-polynomial"
        with self.assertRaisesRegex(AnalysisError, "requires exact case"):
            self.analyzer.analyze(results)

    def test_nonfinite_error_and_duplicate_case_id_are_rejected(self) -> None:
        invalid = h_results()
        invalid[0]["domain_result"]["l2_error"] = math.inf
        with self.assertRaisesRegex(AnalysisError, "must be finite"):
            self.analyzer.analyze(invalid)

        duplicate = h_results()
        duplicate[-1]["case_id"] = duplicate[0]["case_id"]
        with self.assertRaisesRegex(AnalysisError, "duplicate case_id"):
            self.analyzer.analyze(duplicate)

    def test_configured_analyzer_rejects_results_outside_frozen_plan(self) -> None:
        results = h_results()
        plan = tuple(_planned_case(result) for result in results)
        configured = self.analyzer.configure(plan)
        configured.analyze(results)
        changed = deepcopy(results)
        changed[0]["configuration"]["build"]["profile"] = "release"
        with self.assertRaisesRegex(AnalysisError, "frozen spatial plan"):
            configured.analyze(changed)

    def test_h_json_and_csv_preserve_pair_audit_and_local_orders(self) -> None:
        report = self.analyzer.analyze(h_results())
        with tempfile.TemporaryDirectory(prefix="ads-h-analysis-") as temporary:
            written = write_report(report, Path(temporary))
            document = json.loads(written.json_path.read_text(encoding="utf-8"))
            point = document["series"][0]["points"][0]
            self.assertIn("coarse_case_id", point)
            self.assertIn("fine_case_id", point)
            self.assertIn("temporal_error_fraction", point["temporal_pair"]["l2"])
            self.assertIn(
                "temporal_field_difference", point["temporal_pair"]["l2"]
            )
            with written.csv_path.open(encoding="utf-8", newline="") as stream:
                rows = list(csv.DictReader(stream))
        self.assertEqual(
            sum(row["row_kind"] == "point" for row in rows), 8
        )
        self.assertEqual(
            sum(row["row_kind"] == "local-order" for row in rows), 6
        )
        self.assertEqual(
            sum(row["row_kind"] == "regression" for row in rows), 2
        )


class PConvergenceAnalysisTests(unittest.TestCase):
    def setUp(self) -> None:
        self.analyzer = PConvergenceAnalyzer()

    def test_p_reports_reduction_ratios_without_a_polynomial_order(self) -> None:
        report = self.analyzer.analyze(p_results())
        self.assertEqual(report.status, "passed")
        series = report.series[0]
        self.assertEqual(series["sequence_kind"], "isotropic-degree")
        self.assertEqual(
            [point["trial_degree"] for point in series["points"]],
            [[1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4]],
        )
        for metric in series["metrics"].values():
            self.assertIsNone(metric["regression"])
            self.assertEqual(len(metric["error_ratios"]), 3)
            for estimate in metric["error_ratios"]:
                self.assertAlmostEqual(estimate["error_ratio"], 5.0, places=12)
                self.assertNotIn("order", estimate)

    def test_p_time_dominated_endpoint_has_no_reduction_ratio(self) -> None:
        results = p_results()
        target = next(
            result
            for result in results
            if result["configuration"]["spaces"]["trial_degree"] == [2, 2, 2]
            and result["configuration"]["time"]["steps"] == 10
        )
        fine = next(
            result
            for result in results
            if result["configuration"]["spaces"]["trial_degree"] == [2, 2, 2]
            and result["configuration"]["time"]["steps"] == 20
        )
        target["domain_result"]["linf_error"] = 1.25 * fine["domain_result"][
            "linf_error"
        ]
        _attach_field_samples(target)
        report = self.analyzer.analyze(results)
        self.assertEqual(report.status, "failed")
        metric = report.series[0]["metrics"]["linf"]
        self.assertEqual(metric["status"], "unreliable")
        self.assertIsNone(metric["error_ratios"][0]["error_ratio"])
        self.assertIsNone(metric["error_ratios"][1]["error_ratio"])

    def test_p_plateau_is_detected_without_inventing_an_order(self) -> None:
        errors = (0.2, 0.04, 0.008, 0.008, 0.008)
        report = self.analyzer.analyze(p_results((1, 2, 3, 4, 5), errors=errors))
        self.assertEqual(report.status, "passed")
        for metric in report.series[0]["metrics"].values():
            self.assertTrue(metric["plateau_detected"])
            self.assertEqual(metric["plateau_start_trial_degree"], [4, 4, 4])
            self.assertEqual(metric["strict_pre_plateau_decrease_count"], 2)
            self.assertIsNone(metric["regression"])

    def test_flat_p_sequence_is_not_accepted_as_a_successful_plateau(self) -> None:
        report = self.analyzer.analyze(
            p_results((1, 2, 3, 4, 5), errors=(0.1, 0.1, 0.1, 0.1, 0.1))
        )
        self.assertEqual(report.status, "failed")
        for metric in report.series[0]["metrics"].values():
            self.assertTrue(metric["plateau_detected"])
            self.assertEqual(metric["strict_pre_plateau_decrease_count"], 0)
            self.assertEqual(metric["status"], "failed")

    def test_nonmonotonic_p_error_outside_plateau_fails(self) -> None:
        report = self.analyzer.analyze(
            p_results((1, 2, 3, 4), errors=(0.2, 0.04, 0.08, 0.01))
        )
        self.assertEqual(report.status, "failed")
        self.assertTrue(
            report.series[0]["metrics"]["l2"]["fatal_nonmonotonic_anomalies"]
        )

    def test_enrichment_branches_are_never_mixed(self) -> None:
        report = self.analyzer.analyze(
            p_results((1, 2, 3), enrichment=1)
            + p_results((1, 2, 3), enrichment=2)
        )
        self.assertEqual(report.summary["series_count"], 2)
        enrichments = {
            tuple(series["degree_cohort"]["enrichment"]) for series in report.series
        }
        self.assertEqual(enrichments, {(1, 1, 1), (2, 2, 2)})

    def test_mesh_change_is_not_hidden_inside_a_p_series(self) -> None:
        results = p_results()
        for result in results:
            if result["configuration"]["spaces"]["trial_degree"] == [4, 4, 4]:
                result["configuration"]["mesh"]["elements"] = [8, 8, 8]
        with self.assertRaisesRegex(AnalysisError, "p convergence requires at least"):
            self.analyzer.analyze(results)

    def test_anisotropic_rotations_preserve_vectors_and_are_not_ordered(self) -> None:
        rotations = ((3, 4, 5), (4, 5, 3), (5, 3, 4))
        results: list[dict[str, object]] = []
        for index, trial in enumerate(rotations):
            test = tuple(value + 1 for value in trial)
            error = 0.01 * (1.0 + 0.02 * index)
            results.extend(
                _paired_results(
                    family="p",
                    level_label=f"rotation-{index}",
                    mesh=(4, 4, 4),
                    trial=trial,
                    test=test,
                    fine_l2=error,
                    fine_linf=2.0 * error,
                )
            )
        report = self.analyzer.analyze(results)
        self.assertEqual(report.status, "passed")
        self.assertEqual(report.summary["sequence_kinds"], ["anisotropic-rotations"])
        series = report.series[0]
        self.assertEqual(series["acceptance_mode"], "audit-only")
        self.assertEqual(
            {tuple(point["trial_degree"]) for point in series["points"]},
            set(rotations),
        )
        for metric in series["metrics"].values():
            self.assertEqual(metric["acceptance_mode"], "audit-only")
            self.assertIn("no prescribed equality", metric["acceptance_reason"])
            self.assertEqual(metric["error_ratios"], [])
            self.assertIsNone(metric["regression"])
            self.assertGreater(metric["rotation_error_spread"], 1.0)

    def test_invalid_component_wise_enrichment_is_rejected(self) -> None:
        results = p_results((1, 2))
        for result in results[:2]:
            result["configuration"]["spaces"]["test_degree"] = [2, 3, 2]
        with self.assertRaisesRegex(AnalysisError, "p_test=p_trial"):
            self.analyzer.analyze(results)

    def test_p_csv_serializes_degree_vectors_without_flattening(self) -> None:
        report = self.analyzer.analyze(p_results((1, 2, 3)))
        with tempfile.TemporaryDirectory(prefix="ads-p-analysis-") as temporary:
            written = write_report(report, Path(temporary))
            with written.csv_path.open(encoding="utf-8", newline="") as stream:
                rows = list(csv.DictReader(stream))
        point_rows = [row for row in rows if row["row_kind"] == "point"]
        self.assertEqual({row["trial_degree"] for row in point_rows}, {
            "[1,1,1]",
            "[2,2,2]",
            "[3,3,3]",
        })


if __name__ == "__main__":
    unittest.main()
