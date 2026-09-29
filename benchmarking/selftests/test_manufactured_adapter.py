from __future__ import annotations

from dataclasses import replace
import json
import math
from pathlib import Path
import unittest

from ads_benchmark.catalog import build_catalog
from ads_benchmark.components.manufactured import (
    ManufacturedTransientAdapter,
    RESULT_PREFIX,
)
from ads_benchmark.framework.config import load_profiles
from ads_benchmark.framework.errors import ExecutionError
from ads_benchmark.framework.model import ExecutionContext
from ads_benchmark.framework.planner import Planner


BENCHMARKING_ROOT = Path(__file__).resolve().parents[1]


def valid_result(problem: str = "igrm_l2", scheme: str = "dg") -> dict[str, object]:
    return {
        "schema_version": 1,
        "kind": "ads-manufactured-transient-result",
        "exact_case": "temporal-polynomial",
        "problem": problem,
        "scheme": scheme,
        "requested_final_time": 0.1,
        "actual_final_time": 0.1,
        "time_step": 0.025,
        "steps": 4,
        "initial_l2_error": 1.0e-15,
        "initial_linf_error": 2.0e-15,
        "l2_error": 1.0e-4,
        "linf_error": 2.0e-4,
        "solution_l2_norm": 0.204,
        "field_checksum": 12.5,
        "sample_points_per_axis": 17,
        "field_samples_written": False,
        "physical_step_wall_seconds": 0.5,
        "solver_status": 0,
    }


def tagged(document: object) -> str:
    return RESULT_PREFIX + json.dumps(document, separators=(",", ":")) + "\n"


class ManufacturedAdapterTests(unittest.TestCase):
    def setUp(self) -> None:
        self.catalog = build_catalog()
        profiles = load_profiles(BENCHMARKING_ROOT / "configs")
        self.plan = Planner(profiles, self.catalog).plan("smoke")
        self.adapter = ManufacturedTransientAdapter("igrm_l2")

    def test_all_registered_adapters_build_the_shared_harness_argv(self) -> None:
        repository = Path("/repository").resolve()
        context = ExecutionContext(repository, Path("/case").resolve())
        by_problem = {case.spec.problem: case.spec for case in self.plan.cases}

        for problem in ("igrm_l2", "igrm_heat", "pure_diffusion_igrm"):
            with self.subTest(problem=problem):
                adapter = self.catalog.adapters.get(problem)
                case = replace(by_problem[problem], scheme="dg")
                command = tuple(adapter.build_payload_command(case, context))
                self.assertEqual(
                    command,
                    (
                        str(
                            repository
                            / "benchmarking"
                            / "build"
                            / "debug"
                            / "EXEC"
                            / f"{problem}_manufactured"
                        ),
                        "dg",
                        "0.1",
                        "4",
                        "3",
                        "3",
                        "3",
                        "4",
                        "4",
                        "4",
                        "3",
                        "3",
                        "3",
                        "1",
                        "1",
                        "1",
                        "17",
                        "0",
                    ),
                )

    def test_fractional_final_time_is_converted_for_fortran(self) -> None:
        case = replace(self.plan.cases[0].spec, problem="igrm_l2")
        case = replace(
            case,
            time=replace(case.time, final_time="1/3", time_step="1/12"),
        )
        command = self.adapter.build_payload_command(
            case, ExecutionContext(Path("/repo"), Path("/case"))
        )
        self.assertNotIn("/", command[2])
        self.assertAlmostEqual(float(command[2]), 1.0 / 3.0)

    def test_payload_rejects_non_cubic_trial_space(self) -> None:
        case = next(
            item.spec for item in self.plan.cases if item.spec.problem == "igrm_l2"
        )
        case = replace(
            case,
            test_degree=(3, 3, 3),
            trial_degree=(2, 2, 2),
        )
        with self.assertRaisesRegex(ValueError, "trial degree at least 3"):
            self.adapter.build_payload_command(
                case, ExecutionContext(Path("/repo"), Path("/case"))
            )

    def test_spatial_case_accepts_linear_anisotropic_trial_space(self) -> None:
        case = next(
            item.spec for item in self.plan.cases if item.spec.problem == "igrm_l2"
        )
        case = replace(
            case,
            family="p",
            exact_case="spatial-cosine",
            test_degree=(2, 3, 4),
            trial_degree=(1, 2, 3),
        )
        command = self.adapter.build_payload_command(
            case, ExecutionContext(Path("/repo"), Path("/case"))
        )
        self.assertEqual(command[-1], "spatial-cosine")
        self.assertEqual(len(command), 19)
        self.assertEqual(command[7:13], ("2", "3", "4", "1", "2", "3"))

    def test_experiment_family_requires_its_registered_exact_case(self) -> None:
        case = next(
            item.spec for item in self.plan.cases if item.spec.problem == "igrm_l2"
        )
        with self.assertRaisesRegex(ValueError, "temporal family requires"):
            self.adapter.validate_case(replace(case, exact_case="spatial-cosine"))
        with self.assertRaisesRegex(ValueError, "h family requires"):
            self.adapter.validate_case(replace(case, family="h"))

    def test_valid_result_is_normalized_and_json_serializable(self) -> None:
        parsed = self.adapter.parse_result(
            "solver prelude\n" + tagged(valid_result()), "diagnostic\n"
        )
        self.assertEqual(parsed["problem"], "igrm_l2")
        self.assertEqual(parsed["steps"], 4)
        self.assertIsInstance(parsed["l2_error"], float)
        json.dumps(parsed, allow_nan=False)

        spatial = valid_result()
        spatial["exact_case"] = "spatial-cosine"
        self.assertEqual(
            self.adapter.parse_result(tagged(spatial), "")["exact_case"],
            "spatial-cosine",
        )

    def test_result_validation_binds_every_case_identity_field(self) -> None:
        case = next(
            item.spec
            for item in self.plan.cases
            if item.spec.problem == "igrm_l2" and item.spec.scheme == "dg"
        )
        parsed = dict(self.adapter.parse_result(tagged(valid_result()), ""))
        self.adapter.validate_result(case, parsed)

        within_round_trip_tolerance = dict(parsed)
        within_round_trip_tolerance["requested_final_time"] = math.nextafter(
            parsed["requested_final_time"], math.inf
        )
        self.adapter.validate_result(case, within_round_trip_tolerance)

        mutations = (
            ("problem", "igrm_heat"),
            ("scheme", "pr"),
            ("exact_case", "another-case"),
            ("requested_final_time", 0.100000000001),
            ("actual_final_time", 0.100000000001),
            ("time_step", 0.025000000001),
            ("steps", 8),
            ("sample_points_per_axis", 9),
            ("field_samples_written", True),
        )
        for field, value in mutations:
            with self.subTest(field=field):
                candidate = dict(parsed)
                candidate[field] = value
                with self.assertRaisesRegex(
                    ExecutionError, f"result {field} does not match"
                ):
                    self.adapter.validate_result(case, candidate)

    def test_result_requires_one_stdout_tag_and_one_strict_object(self) -> None:
        cases = (
            ("", "", "exactly one"),
            (tagged(valid_result()) * 2, "", "exactly one"),
            ("", tagged(valid_result()), "must be written to stdout"),
            (RESULT_PREFIX + "[]\n", "", "JSON object"),
            (RESULT_PREFIX + "{not-json}\n", "", "invalid tagged"),
            (
                RESULT_PREFIX + '{"schema_version":1,"schema_version":1}\n',
                "",
                "duplicate result key",
            ),
        )
        for stdout, stderr, message in cases:
            with self.subTest(message=message):
                with self.assertRaisesRegex(ExecutionError, message):
                    self.adapter.parse_result(stdout, stderr)

    def test_result_rejects_key_type_range_and_identity_errors(self) -> None:
        mutations: list[tuple[dict[str, object], str]] = []
        missing = valid_result()
        del missing["l2_error"]
        mutations.append((missing, "missing l2_error"))
        extra = valid_result()
        extra["surprise"] = 1
        mutations.append((extra, "unknown surprise"))
        wrong_problem = valid_result(problem="igrm_heat")
        mutations.append((wrong_problem, "result problem"))
        wrong_scheme = valid_result(scheme="fe")
        mutations.append((wrong_scheme, "result scheme"))
        bad_status = valid_result()
        bad_status["solver_status"] = 1
        mutations.append((bad_status, "must be zero"))
        bool_steps = valid_result()
        bool_steps["steps"] = True
        mutations.append((bool_steps, "positive integer"))
        negative_error = valid_result()
        negative_error["l2_error"] = -1.0
        mutations.append((negative_error, "must be nonnegative"))
        bad_samples = valid_result()
        bad_samples["sample_points_per_axis"] = 1
        mutations.append((bad_samples, "between 2 and 257"))
        bad_write_flag = valid_result()
        bad_write_flag["field_samples_written"] = 0
        mutations.append((bad_write_flag, "must be a boolean"))
        inconsistent_step = valid_result()
        inconsistent_step["time_step"] = 0.02
        mutations.append((inconsistent_step, "time_step\\*steps"))
        inconsistent_final = valid_result()
        inconsistent_final["actual_final_time"] = 0.2
        mutations.append((inconsistent_final, "actual_final_time"))

        for document, message in mutations:
            with self.subTest(message=message):
                with self.assertRaisesRegex(ExecutionError, message):
                    self.adapter.parse_result(tagged(document), "")

    def test_result_rejects_every_nonfinite_number_spelling(self) -> None:
        for spelling in ("NaN", "Infinity", "-Infinity", "1e10000"):
            with self.subTest(spelling=spelling):
                document = json.dumps(valid_result(), separators=(",", ":"))
                document = document.replace(
                    '"l2_error":0.0001', f'"l2_error":{spelling}'
                )
                with self.assertRaisesRegex(ExecutionError, "non-finite|finite number"):
                    self.adapter.parse_result(RESULT_PREFIX + document + "\n", "")


if __name__ == "__main__":
    unittest.main(verbosity=2)
