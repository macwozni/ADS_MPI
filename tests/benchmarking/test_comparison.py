from __future__ import annotations

import csv
from fractions import Fraction
import hashlib
import io
import math
import unittest

from ads_benchmark.comparison import (
    COMPARISON_KIND,
    COMPARISON_SCHEMA_VERSION,
    ComparisonError,
    ComparisonPolicy,
    RunSnapshot,
    compare_snapshots,
)
from ads_benchmark.framework.model import (
    CaseSpec,
    MeasurementSpec,
    MpiSpec,
    PlannedCase,
    RepositoryState,
    SamplingSpec,
    TimeSpec,
    WeakScalingSpec,
)
from ads_benchmark.framework.planner import FrozenPlan, canonical_json, case_identity


COMPATIBLE_PROVENANCE = {
    "schema_version": 1,
    "machine": {"cpu": "synthetic-cpu", "cores": 8},
    "toolchain": {"compiler": "fake-fortran 1", "flags": ["-O3"]},
    "runtime": {"mpi": "fake-mpi 1", "launcher": "mpiexec"},
}


def _case(
    *,
    scheme: str = "dg",
    samples: int = 5,
    family: str = "strong",
    weak_role: str | None = None,
) -> PlannedCase:
    weak_scaling = None
    if weak_role is not None:
        weak_scaling = WeakScalingSpec(
            local_elements=(2, 2, 2),
            workload_basis="per-rank",
            role=weak_role,
        )
    spec = CaseSpec(
        family=family,
        problem="igrm_l2",
        scheme=scheme,
        exact_case="spatial-cosine",
        time=TimeSpec(final_time="0.1", time_step="0.025", steps=4),
        mesh=(2, 2, 2),
        test_degree=(3, 3, 3),
        trial_degree=(2, 2, 2),
        mpi=MpiSpec(ranks=1, process_grid=(1, 1, 1)),
        openmp_threads=1,
        sampling=SamplingSpec(points_per_axis=3, write_samples=True),
        measurement=MeasurementSpec(
            warmups=2,
            samples=samples,
            timeout_seconds="60",
            minimum_sample_seconds="0.001",
        ),
        build_profile="release",
        launcher="mpi",
        openmp_dynamic=False,
        openmp_proc_bind="close",
        openmp_places="cores",
        weak_scaling=weak_scaling,
    )
    case_id, _ = case_identity(spec)
    return PlannedCase(case_id=case_id, spec=spec)


def _configuration_hash(cases: tuple[PlannedCase, ...]) -> str:
    encoded = canonical_json([case.spec.to_dict() for case in cases])
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _plan(
    run_id: str,
    cases: tuple[PlannedCase, ...],
    *,
    config_hash: str | None = None,
) -> FrozenPlan:
    return FrozenPlan(
        run_id=run_id,
        profile_name="synthetic-ab",
        profile_description="synthetic A/B comparison profile",
        filters={},
        config_hash=config_hash or _configuration_hash(cases),
        repository=RepositoryState(
            commit=("a" if run_id == "baseline" else "b") * 40,
            dirty=False,
            worktree_fingerprint=(
                "sha256:" + ("a" if run_id == "baseline" else "b") * 64
            ),
        ),
        cases=cases,
    )


def _checksum(values: list[float]) -> float:
    ordinary = math.fsum(values)
    weighted = math.fsum(
        (index + 1) * value for index, value in enumerate(values)
    )
    return ordinary + weighted / (len(values) + 1)


def _result(
    case: PlannedCase,
    samples: tuple[float, ...],
    *,
    perturbation: float = 0.0,
    reliable: bool = True,
) -> dict[str, object]:
    if len(samples) != case.spec.measurement.samples:
        raise AssertionError("test sample count differs from the synthetic plan")
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
    numerical[1] += perturbation
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

    return {
        "schema_version": 1,
        "kind": "ads-benchmark-case-result",
        "case_id": case.case_id,
        "status": "passed",
        "configuration": case.spec.to_dict(),
        "timing": {
            "wall_seconds": math.fsum(samples),
            "metric": "physical_step_wall_seconds",
            "warmup_samples": [samples[0], samples[0]],
            "measured_samples": list(samples),
            "warmup_process_wall_seconds": [samples[0], samples[0]],
            "measured_process_wall_seconds": list(samples),
            "minimum_reliable_seconds": 0.001,
            "reliable": reliable,
            "openmp_environment": {
                "OMP_NUM_THREADS": "1",
                "OMP_DYNAMIC": "FALSE",
                "OMP_PROC_BIND": "close",
                "OMP_PLACES": "cores",
            },
        },
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
        "analysis_artifacts": {"field_samples_csv": stream.getvalue()},
    }


def _snapshot(
    run_id: str,
    cases: tuple[PlannedCase, ...],
    results: tuple[dict[str, object], ...],
    *,
    provenance: dict[str, object] | None = COMPATIBLE_PROVENANCE,
    config_hash: str | None = None,
) -> RunSnapshot:
    return RunSnapshot.from_loaded(
        _plan(run_id, cases, config_hash=config_hash),
        results,
        execution_provenance=provenance,
    )


class ABComparisonTests(unittest.TestCase):
    def test_identical_runs_emit_versioned_report_and_robust_statistics(self) -> None:
        case = _case()
        samples = (0.9, 1.0, 1.1, 1.0, 1.0)
        baseline = _snapshot("baseline", (case,), (_result(case, samples),))
        candidate = _snapshot("candidate", (case,), (_result(case, samples),))

        report = compare_snapshots(baseline, candidate)
        document = report.to_dict()

        self.assertEqual(report.status, "no-regression-detected")
        self.assertEqual(document["schema_version"], COMPARISON_SCHEMA_VERSION)
        self.assertEqual(document["kind"], COMPARISON_KIND)
        self.assertEqual(document["summary"]["compared_timing_count"], 1)
        timing = document["cases"][0]["timing"]
        self.assertEqual(timing["status"], "no-regression-detected")
        self.assertEqual(timing["median_ratio_candidate_over_baseline"], 1.0)
        self.assertEqual(timing["baseline"]["count"], 5)
        self.assertEqual(timing["baseline"]["median_seconds"], 1.0)
        self.assertEqual(timing["baseline"]["median_absolute_deviation_seconds"], 0.0)
        self.assertAlmostEqual(timing["baseline"]["range_seconds"], 0.2)
        self.assertTrue(document["cases"][0]["numerical"]["passed"])

    def test_configurable_threshold_controls_candidate_over_baseline_ratio(self) -> None:
        case = _case()
        baseline_samples = (0.98, 0.99, 1.0, 1.01, 1.02)
        baseline = _snapshot(
            "baseline", (case,), (_result(case, baseline_samples),)
        )

        below = tuple(value * 1.04 for value in baseline_samples)
        below_report = compare_snapshots(
            baseline,
            _snapshot("candidate", (case,), (_result(case, below),)),
            ComparisonPolicy(regression_threshold=0.05, minimum_samples=5),
        )
        self.assertEqual(below_report.status, "no-regression-detected")

        above = tuple(value * 1.06 for value in baseline_samples)
        above_report = compare_snapshots(
            baseline,
            _snapshot("candidate", (case,), (_result(case, above),)),
            ComparisonPolicy(regression_threshold=0.05, minimum_samples=5),
        )
        self.assertEqual(above_report.status, "regression")
        self.assertEqual(above_report.summary["regression_count"], 1)
        self.assertAlmostEqual(
            above_report.cases[0]["timing"][
                "median_ratio_candidate_over_baseline"
            ],
            1.06,
        )

    def test_any_numerical_mismatch_blocks_all_timing_statistics(self) -> None:
        first = _case(scheme="dg")
        second = _case(scheme="pr")
        cases = (first, second)
        baseline_results = (
            _result(first, (1.0,) * 5),
            _result(second, (1.0,) * 5),
        )
        candidate_results = (
            _result(first, (2.0,) * 5),
            _result(second, (2.0,) * 5, perturbation=1.0e-4),
        )

        report = compare_snapshots(
            _snapshot("baseline", cases, baseline_results),
            _snapshot("candidate", cases, candidate_results),
        )

        self.assertEqual(report.status, "numerical-mismatch")
        self.assertEqual(report.summary["numerical_failure_count"], 1)
        self.assertEqual(report.summary["regression_count"], 0)
        self.assertEqual(
            {case["timing"]["status"] for case in report.cases},
            {"blocked-by-numerical-gate"},
        )
        self.assertTrue(
            all(
                "median_ratio_candidate_over_baseline" not in case["timing"]
                for case in report.cases
            )
        )

    def test_one_sample_is_inconclusive_and_never_a_regression(self) -> None:
        case = _case(samples=1)
        report = compare_snapshots(
            _snapshot("baseline", (case,), (_result(case, (1.0,)),)),
            _snapshot("candidate", (case,), (_result(case, (4.0,)),)),
            ComparisonPolicy(minimum_samples=3),
        )

        self.assertEqual(report.status, "inconclusive")
        self.assertEqual(report.summary["regression_count"], 0)
        timing = report.cases[0]["timing"]
        self.assertEqual(timing["status"], "inconclusive")
        self.assertEqual(timing["median_ratio_candidate_over_baseline"], 4.0)
        self.assertTrue(
            any("fewer than" in reason for reason in timing["reasons"])
        )

    def test_legacy_missing_execution_provenance_is_never_green_or_regression(self) -> None:
        case = _case()
        report = compare_snapshots(
            _snapshot(
                "baseline", (case,), (_result(case, (1.0,) * 5),), provenance=None
            ),
            _snapshot(
                "candidate", (case,), (_result(case, (2.0,) * 5),), provenance=None
            ),
        )

        self.assertEqual(report.status, "inconclusive")
        self.assertEqual(report.compatibility["execution_provenance"], "missing")
        self.assertEqual(report.summary["regression_count"], 0)
        self.assertEqual(report.cases[0]["timing"]["status"], "inconclusive")
        self.assertEqual(
            report.cases[0]["timing"]["median_ratio_candidate_over_baseline"],
            2.0,
        )

    def test_incompatible_hash_provenance_and_result_set_are_rejected(self) -> None:
        case = _case()
        result = _result(case, (1.0,) * 5)
        baseline = _snapshot("baseline", (case,), (result,))

        with self.subTest(reason="config hash"):
            report = compare_snapshots(
                baseline,
                _snapshot(
                    "candidate",
                    (case,),
                    (_result(case, (1.0,) * 5),),
                    config_hash="f" * 64,
                ),
            )
            self.assertEqual(report.status, "incompatible")
            self.assertIn("config_hash", " ".join(report.diagnostics))

        with self.subTest(reason="provenance"):
            incompatible = {
                **COMPATIBLE_PROVENANCE,
                "machine": {"cpu": "different-cpu", "cores": 8},
            }
            report = compare_snapshots(
                baseline,
                _snapshot(
                    "candidate",
                    (case,),
                    (_result(case, (1.0,) * 5),),
                    provenance=incompatible,
                ),
            )
            self.assertEqual(report.status, "incompatible")
            self.assertEqual(
                report.compatibility["execution_provenance"], "incompatible"
            )

        with self.subTest(reason="missing result"):
            report = compare_snapshots(
                baseline,
                _snapshot("candidate", (case,), ()),
            )
            self.assertEqual(report.status, "incompatible")
            self.assertIn("miss cases", " ".join(report.diagnostics))

    def test_unreliable_timing_is_descriptive_but_inconclusive(self) -> None:
        case = _case()
        report = compare_snapshots(
            _snapshot("baseline", (case,), (_result(case, (1.0,) * 5),)),
            _snapshot(
                "candidate",
                (case,),
                (_result(case, (2.0,) * 5, reliable=False),),
            ),
        )

        self.assertEqual(report.status, "inconclusive")
        self.assertEqual(report.summary["regression_count"], 0)
        self.assertIn(
            "candidate timing is marked unreliable",
            report.cases[0]["timing"]["reasons"],
        )

    def test_weak_field_reference_is_gated_but_excluded_from_timing(self) -> None:
        measurement = _case(
            scheme="dg", family="weak", weak_role="measurement"
        )
        helper = _case(
            scheme="pr", family="weak", weak_role="field-reference"
        )
        cases = (measurement, helper)
        baseline_results = tuple(
            _result(case, (1.0,) * 5) for case in cases
        )
        candidate_results = tuple(
            _result(case, (1.0,) * 5) for case in cases
        )

        report = compare_snapshots(
            _snapshot("baseline", cases, baseline_results),
            _snapshot("candidate", cases, candidate_results),
        )

        self.assertEqual(report.status, "no-regression-detected")
        self.assertEqual(report.summary["timing_case_count"], 1)
        self.assertEqual(report.summary["timing_not_applicable_count"], 1)
        by_id = {case["case_id"]: case for case in report.cases}
        self.assertTrue(by_id[helper.case_id]["numerical"]["passed"])
        self.assertEqual(
            by_id[helper.case_id]["timing"]["status"], "not-applicable"
        )

    def test_policy_and_snapshot_boundaries_reject_unsafe_values(self) -> None:
        for policy in (
            lambda: ComparisonPolicy(minimum_samples=2),
            lambda: ComparisonPolicy(regression_threshold=float("nan")),
            lambda: ComparisonPolicy(regression_threshold=-0.01),
            lambda: ComparisonPolicy(absolute_field_tolerance=-1.0),
        ):
            with self.subTest(policy=policy):
                with self.assertRaises(ComparisonError):
                    policy()

        case = _case()
        with self.assertRaisesRegex(ComparisonError, "iterable"):
            RunSnapshot.from_loaded(
                _plan("baseline", (case,)),
                {"case_id": case.case_id},
                execution_provenance=COMPATIBLE_PROVENANCE,
            )
        with self.assertRaisesRegex(ComparisonError, "strict JSON"):
            RunSnapshot.from_loaded(
                _plan("baseline", (case,)),
                (_result(case, (1.0,) * 5),),
                execution_provenance={"bad": float("nan")},
            )

    def test_configuration_mismatch_inside_a_loaded_result_is_incompatible(self) -> None:
        case = _case()
        baseline_result = _result(case, (1.0,) * 5)
        candidate_result = _result(case, (1.0,) * 5)
        candidate_result["configuration"] = {
            **case.spec.to_dict(),
            "scheme": "pr",
        }

        report = compare_snapshots(
            _snapshot("baseline", (case,), (baseline_result,)),
            _snapshot("candidate", (case,), (candidate_result,)),
        )

        self.assertEqual(report.status, "incompatible")
        self.assertIn("differs from its frozen plan", " ".join(report.diagnostics))


if __name__ == "__main__":
    unittest.main(verbosity=2)
