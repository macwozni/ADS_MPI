from __future__ import annotations

import csv
from dataclasses import replace
from fractions import Fraction
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

from ads_benchmark.analysis import (
    AnalysisError,
    AnalysisPipeline,
    TemporalConvergenceAnalyzer,
    load_run_results,
    write_report,
)
from ads_benchmark.catalog import build_catalog
from ads_benchmark.components.manufactured import RESULT_PREFIX
from ads_benchmark.framework.config import load_profiles
from ads_benchmark.framework.executor import Executor
from ads_benchmark.framework.model import (
    ExecutionContext,
    PlannedCase,
    SamplingSpec,
)
from ads_benchmark.framework.planner import Planner
from ads_benchmark.framework.storage import ResultStore
from benchmark_paths import CONFIG_DIRECTORY


STEPS = (4, 8, 16, 32, 64, 128, 256, 512)


def synthetic_results(
    scheme: str,
    order: float,
    *,
    error_scale: float = 1.0,
) -> list[dict[str, object]]:
    results: list[dict[str, object]] = []
    for steps in STEPS:
        dt = 0.1 / steps
        configuration = {
            "family": "temporal",
            "problem": "synthetic_problem",
            "scheme": scheme,
            "exact_case": "temporal-polynomial",
            "time": {
                "final_time": "0.1",
                "time_step": format(dt, ".17g"),
                "steps": steps,
            },
            "mesh": {"elements": [3, 3, 3]},
            "spaces": {
                "test_degree": [4, 4, 4],
                "trial_degree": [3, 3, 3],
            },
            "mpi": {"ranks": 1, "process_grid": [1, 1, 1]},
            "openmp": {"threads": 1},
            "sampling": {"points_per_axis": 17, "write_samples": False},
            "measurement": {
                "warmups": 0,
                "samples": 1,
                "timeout_seconds": "60",
            },
            "build": {"profile": "debug"},
            "launcher": "direct",
        }
        # Use the exact decimal string from the configuration to keep the
        # synthetic T=N*dt contract exact, just as the planner does.
        configuration["time"]["time_step"] = str(0.1 / steps)  # type: ignore[index]
        l2_error = error_scale * dt**order
        results.append(
            {
                "schema_version": 1,
                "kind": "ads-benchmark-case-result",
                "case_id": f"synthetic-{scheme}-{steps}",
                "status": "passed",
                "configuration": configuration,
                "timing": {"wall_seconds": 0.01},
                "domain_result": {
                    "schema_version": 1,
                    "kind": "ads-manufactured-transient-result",
                    "problem": "synthetic_problem",
                    "scheme": scheme,
                    "exact_case": "temporal-polynomial",
                    "requested_final_time": 0.1,
                    "actual_final_time": 0.1,
                    "time_step": dt,
                    "steps": steps,
                    "l2_error": l2_error,
                    "linf_error": 2.5 * l2_error,
                    "solution_l2_norm": 0.2,
                    "solver_status": 0,
                },
            }
        )
    return results


class TemporalConvergenceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.analyzer = TemporalConvergenceAnalyzer()

    def test_second_order_has_seven_l2_and_linf_local_orders(self) -> None:
        report = self.analyzer.analyze(synthetic_results("dg", 2.0))
        self.assertEqual(report.status, "passed")
        self.assertEqual(report.summary["local_orders_per_metric"], 7)
        series = report.series[0]
        for metric_name in ("l2", "linf"):
            metric = series["metrics"][metric_name]
            self.assertEqual(len(metric["local_orders"]), 7)
            for estimate in metric["local_orders"]:
                self.assertAlmostEqual(estimate["order"], 2.0, places=12)
            self.assertAlmostEqual(metric["regression"]["order"], 2.0, places=12)

    def test_first_order_pr_and_be_pass_their_explicit_intervals(self) -> None:
        for scheme in ("pr", "be"):
            with self.subTest(scheme=scheme):
                report = self.analyzer.analyze(synthetic_results(scheme, 1.0))
                self.assertEqual(report.status, "passed")
                for metric in report.series[0]["metrics"].values():
                    self.assertEqual(metric["accepted_interval"], [0.8, 1.2])
                    self.assertAlmostEqual(
                        metric["regression"]["order"], 1.0, places=12
                    )

    def test_wrong_formal_order_produces_a_failed_report(self) -> None:
        report = self.analyzer.analyze(synthetic_results("dg", 1.0))
        self.assertEqual(report.status, "failed")
        self.assertEqual(report.summary["failed_series"], 1)

    def test_catastrophic_coarse_pr_point_cannot_hide_behind_clean_tail(self) -> None:
        results = synthetic_results("pr", 1.0)
        l2_errors = (
            2.501401e8,
            3.737280e-2,
            2.320375e-4,
            2.257949e-5,
            1.124023e-5,
            5.605335e-6,
            2.800000e-6,
            1.400000e-6,
        )
        linf_errors = (
            8.649282e9,
            2.030622e0,
            6.827116e-3,
            1.125050e-4,
            5.585438e-5,
            2.779671e-5,
            1.389000e-5,
            6.945000e-6,
        )
        for result, l2_error, linf_error in zip(
            results, l2_errors, linf_errors, strict=True
        ):
            result["domain_result"]["l2_error"] = l2_error
            result["domain_result"]["linf_error"] = linf_error

        report = self.analyzer.analyze(results)
        self.assertEqual(report.status, "failed")
        for metric in report.series[0]["metrics"].values():
            self.assertGreater(metric["regression"]["order"], 0.8)
            self.assertLess(metric["regression"]["order"], 1.2)
            self.assertEqual(metric["maximum_local_order"], 4.0)
            self.assertTrue(metric["implausible_local_orders"])

    def test_finest_plateau_is_reported_and_excluded_from_regression(self) -> None:
        results = synthetic_results("dg", 2.0)
        plateau_l2 = results[5]["domain_result"]["l2_error"]
        plateau_linf = results[5]["domain_result"]["linf_error"]
        for index in (6, 7):
            results[index]["domain_result"]["l2_error"] = plateau_l2
            results[index]["domain_result"]["linf_error"] = plateau_linf
        report = self.analyzer.analyze(results)
        self.assertEqual(report.status, "passed")
        for metric in report.series[0]["metrics"].values():
            self.assertTrue(metric["plateau_detected"])
            self.assertEqual(metric["plateau_start_steps"], 256)
            self.assertEqual(metric["fatal_nonmonotonic_anomalies"], [])
            self.assertEqual(len(metric["local_orders"]), 7)
            self.assertAlmostEqual(metric["regression"]["order"], 2.0, places=12)

    def test_missing_level_is_rejected(self) -> None:
        with self.assertRaisesRegex(AnalysisError, "incomplete temporal levels"):
            self.analyzer.analyze(synthetic_results("dg", 2.0)[:-1])

    def test_nan_and_infinity_are_rejected(self) -> None:
        for invalid in (float("nan"), float("inf"), float("-inf")):
            with self.subTest(invalid=invalid):
                results = synthetic_results("dg", 2.0)
                results[3]["domain_result"]["l2_error"] = invalid
                with self.assertRaisesRegex(AnalysisError, "must be finite"):
                    self.analyzer.analyze(results)

    def test_nonmonotonic_analytical_error_is_not_hidden_as_plateau(self) -> None:
        results = synthetic_results("be", 1.0)
        previous = results[-2]["domain_result"]["l2_error"]
        results[-1]["domain_result"]["l2_error"] = 2.0 * previous
        report = self.analyzer.analyze(results)
        self.assertEqual(report.status, "failed")
        l2 = report.series[0]["metrics"]["l2"]
        self.assertEqual(l2["status"], "failed")
        self.assertEqual(len(l2["local_orders"]), 7)
        self.assertEqual(
            l2["fatal_nonmonotonic_anomalies"],
            [{"coarse_steps": 256, "fine_steps": 512, "order": -1.0}],
        )
        self.assertTrue(l2["local_orders"][-1]["fatal_nonmonotonic"])
        self.assertTrue(
            any("fatal nonmonotonic" in message for message in report.diagnostics)
        )
        with tempfile.TemporaryDirectory(prefix="ads-nonmonotonic-report-") as temporary:
            written = write_report(report, Path(temporary))
            document = json.loads(written.json_path.read_text(encoding="utf-8"))
            self.assertEqual(document["status"], "failed")
            with written.csv_path.open(encoding="utf-8", newline="") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(
                sum(row["fatal_nonmonotonic"] == "True" for row in rows), 1
            )

    def test_wrong_final_time_is_rejected(self) -> None:
        results = synthetic_results("pr", 1.0)
        results[2]["domain_result"]["actual_final_time"] = 0.1001
        with self.assertRaisesRegex(AnalysisError, "does not match T"):
            self.analyzer.analyze(results)

    def test_nonzero_solver_status_is_rejected(self) -> None:
        results = synthetic_results("pr", 1.0)
        results[4]["domain_result"]["solver_status"] = 7
        with self.assertRaisesRegex(AnalysisError, "solver_status is nonzero"):
            self.analyzer.analyze(results)

    def test_duplicate_case_id_is_rejected(self) -> None:
        results = synthetic_results("be", 1.0)
        results[-1]["case_id"] = results[0]["case_id"]
        with self.assertRaisesRegex(AnalysisError, "duplicate case_id"):
            self.analyzer.analyze(results)

    def test_inconsistent_final_time_with_individually_valid_points_is_rejected(self) -> None:
        results = synthetic_results("be", 1.0)
        item = results[-1]
        item["configuration"]["time"] = {
            "final_time": "0.2",
            "time_step": str(0.2 / 512),
            "steps": 512,
        }
        item["domain_result"]["requested_final_time"] = 0.2
        item["domain_result"]["actual_final_time"] = 0.2
        item["domain_result"]["time_step"] = 0.2 / 512
        with self.assertRaisesRegex(AnalysisError, "inconsistent final time"):
            self.analyzer.analyze(results)


class AnalysisPipelineAndIoTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory(prefix="ads-temporal-analysis-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def test_analyzer_is_registered_without_family_branches(self) -> None:
        pipeline = AnalysisPipeline()
        pipeline.register(TemporalConvergenceAnalyzer())
        report = pipeline.analyze(
            "temporal-convergence",
            synthetic_results("dg", 2.0),
            source_run="synthetic-run",
        )
        self.assertEqual(report.source_run, "synthetic-run")
        self.assertEqual(report.analyzer, "temporal-convergence")

    def test_catalog_family_selects_registered_analyzer(self) -> None:
        catalog = build_catalog()
        family = catalog.families.get("temporal")
        self.assertEqual(family.analyzer, "temporal-convergence")
        self.assertIsInstance(
            catalog.analyzers.get(family.analyzer), TemporalConvergenceAnalyzer
        )

    def test_frozen_four_level_selection_configures_analyzer(self) -> None:
        catalog = build_catalog()
        plan = Planner(load_profiles(CONFIG_DIRECTORY), catalog).plan("temporal-full")
        exemplar = plan.cases[0].spec
        selected = tuple(
            case
            for case in plan.cases
            if case.spec.problem == exemplar.problem
            and case.spec.scheme == exemplar.scheme
            and case.spec.test_degree == exemplar.test_degree
            and case.spec.trial_degree == exemplar.trial_degree
            and case.spec.time.steps in STEPS[:4]
        )
        configured = catalog.analyzers.get("temporal-convergence").configure(selected)
        self.assertEqual(configured.expected_steps, STEPS[:4])
        self.assertEqual(len(configured.expected_groups), 1)
        with self.assertRaisesRegex(AnalysisError, "frozen planned groups"):
            AnalysisPipeline(catalog.analyzers).analyze(
                "temporal-convergence",
                synthetic_results(exemplar.scheme, 2.0)[:4],
                planned_cases=selected,
            )

    def test_frozen_three_level_selection_is_too_short(self) -> None:
        catalog = build_catalog()
        plan = Planner(load_profiles(CONFIG_DIRECTORY), catalog).plan("temporal-full")
        exemplar = plan.cases[0].spec
        selected = tuple(
            case
            for case in plan.cases
            if case.spec.problem == exemplar.problem
            and case.spec.scheme == exemplar.scheme
            and case.spec.test_degree == exemplar.test_degree
            and case.spec.trial_degree == exemplar.trial_degree
            and case.spec.time.steps in STEPS[:3]
        )
        with self.assertRaisesRegex(AnalysisError, "at least four planned levels"):
            catalog.analyzers.get("temporal-convergence").configure(selected)

    def test_json_and_csv_share_all_seven_local_estimates(self) -> None:
        report = TemporalConvergenceAnalyzer().analyze(
            synthetic_results("dg", 2.0), source_run="synthetic-run"
        )
        written = write_report(report, self.root / "analysis")
        document = json.loads(written.json_path.read_text(encoding="utf-8"))
        self.assertEqual(document["kind"], "ads-benchmark-analysis")
        self.assertEqual(document["status"], "passed")
        with written.csv_path.open(encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 14)
        self.assertEqual({row["metric"] for row in rows}, {"l2", "linf"})
        self.assertTrue(all(math.isfinite(float(row["local_order"])) for row in rows))

    def test_plot_is_the_only_artifact_omitted_without_matplotlib(self) -> None:
        report = TemporalConvergenceAnalyzer().analyze(
            synthetic_results("dg", 2.0)
        )
        with mock.patch.dict(sys.modules, {"matplotlib": None}):
            written = write_report(report, self.root / "analysis", plot=True)
        self.assertTrue(written.json_path.is_file())
        self.assertTrue(written.csv_path.is_file())
        self.assertIsNone(written.plot_path)
        self.assertIn("matplotlib is unavailable", written.plot_message)

    def _owned_case(self) -> tuple[ResultStore, PlannedCase, Path]:
        repository = self.root / "repository"
        repository.mkdir(exist_ok=True)
        catalog = build_catalog()
        case = Planner(load_profiles(CONFIG_DIRECTORY), catalog).plan("smoke").cases[0]
        store = ResultStore(repository)
        store.create_run(
            "analysis-run",
            {
                "schema_version": 1,
                "kind": "ads-benchmark-plan",
                "run_id": "analysis-run",
            },
        )
        case_directory = store.create_case_directory("analysis-run", case.case_id)
        return store, case, case_directory

    @staticmethod
    def _executor(store: ResultStore) -> Executor:
        return Executor(build_catalog(), store, store.repository_root)

    @staticmethod
    def _write_complete_case(
        store: ResultStore, case: PlannedCase, case_directory: Path
    ) -> None:
        catalog = build_catalog()
        context = ExecutionContext(
            repository_root=store.repository_root,
            case_directory=case_directory,
        )
        payload = catalog.adapters.get(case.spec.problem).build_payload_command(
            case.spec, context
        )
        command = catalog.launchers.get(case.spec.launcher).command(
            payload, case.spec
        )
        dt = float(Fraction(case.spec.time.time_step))
        domain_result = {
            "schema_version": 1,
            "kind": "ads-manufactured-transient-result",
            "exact_case": case.spec.exact_case,
            "problem": case.spec.problem,
            "scheme": case.spec.scheme,
            "requested_final_time": float(Fraction(case.spec.time.final_time)),
            "actual_final_time": float(Fraction(case.spec.time.final_time)),
            "time_step": dt,
            "steps": case.spec.time.steps,
            "initial_l2_error": 1.0e-15,
            "initial_linf_error": 2.0e-15,
            "l2_error": dt * dt,
            "linf_error": 2.0 * dt * dt,
            "solution_l2_norm": 0.2,
            "field_checksum": 1.0,
            "sample_points_per_axis": case.spec.sampling.points_per_axis,
            "field_samples_written": case.spec.sampling.write_samples,
            "physical_step_wall_seconds": 0.5,
            "solver_status": 0,
        }
        store.write_status(
            case_directory,
            {
                "schema_version": 1,
                "case_id": case.case_id,
                "state": "passed",
                "started_at": "2026-01-01T00:00:00+00:00",
                "finished_at": "2026-01-01T00:00:01+00:00",
                "duration_seconds": 1.0,
                "return_code": 0,
                "command": list(command),
            },
        )
        store.write_log(
            case_directory,
            "stdout.log",
            RESULT_PREFIX + json.dumps(domain_result, sort_keys=True) + "\n",
        )
        store.write_log(case_directory, "stderr.log", "")
        store.write_result(
            case_directory,
            {
                "schema_version": 1,
                "kind": "ads-benchmark-case-result",
                "case_id": case.case_id,
                "status": "passed",
                "configuration": case.spec.to_dict(),
                "timing": {"wall_seconds": 1.0},
                "domain_result": domain_result,
            },
        )

    def test_run_loader_requires_passed_status_both_logs_and_result(self) -> None:
        store, case, case_directory = self._owned_case()
        self._write_complete_case(store, case, case_directory)
        executor = self._executor(store)
        loaded = load_run_results(executor, "analysis-run", (case,))
        self.assertEqual(len(loaded), 1)
        self.assertNotIn("analysis_artifacts", loaded[0])

        status_path = case_directory / "status.json"
        status = json.loads(status_path.read_text(encoding="utf-8"))
        status["state"] = "running"
        status_path.write_text(json.dumps(status), encoding="utf-8")
        with self.assertRaisesRegex(AnalysisError, "fully verified completed"):
            load_run_results(executor, "analysis-run", (case,))

    def test_run_loader_attaches_requested_samples_lazily_only_in_memory(self) -> None:
        store, case, case_directory = self._owned_case()
        sampled_case = replace(
            case,
            spec=replace(
                case.spec,
                sampling=SamplingSpec(points_per_axis=3, write_samples=True),
            ),
        )
        self._write_complete_case(store, sampled_case, case_directory)
        sample_text = "x,y,z,numerical,exact,error\n0,0,0,1,1,0\n"
        (case_directory / "field_samples.csv").write_text(
            sample_text, encoding="utf-8"
        )
        result_path = case_directory / "result.json"
        persisted_result = result_path.read_bytes()

        loaded = load_run_results(
            self._executor(store), "analysis-run", (sampled_case,)
        )

        artifact = loaded[0]["analysis_artifacts"]["field_samples_csv"]
        self.assertNotIsInstance(artifact, str)
        self.assertEqual(artifact.read_text(), sample_text)
        self.assertEqual(result_path.read_bytes(), persisted_result)

    def test_run_loader_rejects_a_missing_log(self) -> None:
        store, case, case_directory = self._owned_case()
        self._write_complete_case(store, case, case_directory)
        (case_directory / "stderr.log").unlink()
        with self.assertRaisesRegex(AnalysisError, "fully verified completed"):
            load_run_results(self._executor(store), "analysis-run", (case,))

    def test_run_loader_uses_strict_json_reader(self) -> None:
        store, case, case_directory = self._owned_case()
        self._write_complete_case(store, case, case_directory)
        result_path = case_directory / "result.json"
        result_path.write_text('{"case_id":"one","case_id":"two"}', encoding="utf-8")
        with self.assertRaisesRegex(AnalysisError, "fully verified completed"):
            load_run_results(self._executor(store), "analysis-run", (case,))
        result_path.write_text('{"value":NaN}', encoding="utf-8")
        with self.assertRaisesRegex(AnalysisError, "fully verified completed"):
            load_run_results(self._executor(store), "analysis-run", (case,))

    def test_run_loader_rejects_result_forged_independently_of_stdout(self) -> None:
        store, case, case_directory = self._owned_case()
        self._write_complete_case(store, case, case_directory)
        result_path = case_directory / "result.json"
        result = json.loads(result_path.read_text(encoding="utf-8"))
        result["domain_result"]["l2_error"] *= 10.0
        result_path.write_text(json.dumps(result), encoding="utf-8")
        with self.assertRaisesRegex(AnalysisError, "fully verified completed"):
            load_run_results(self._executor(store), "analysis-run", (case,))


if __name__ == "__main__":
    unittest.main(verbosity=2)
