from __future__ import annotations

import copy
from dataclasses import replace
from decimal import Decimal
import json
from pathlib import Path
import tempfile
import unittest

from ads_benchmark.catalog import build_catalog
from ads_benchmark.framework.config import (
    canonical_decimal,
    load_json,
    load_profiles,
    parse_profile,
)
from ads_benchmark.framework.errors import ConfigurationError, ValidationError
from ads_benchmark.framework.model import MeasurementSpec, MpiSpec, TimeSpec
from ads_benchmark.framework.planner import Planner
from ads_benchmark.framework.validation import validate_case


BENCHMARKING_ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIRECTORY = BENCHMARKING_ROOT / "configs"


class ValidationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.catalog = build_catalog()
        profiles = load_profiles(CONFIG_DIRECTORY)
        self.case = Planner(profiles, self.catalog).plan("smoke").cases[0].spec
        self.document = json.loads(
            (CONFIG_DIRECTORY / "smoke.json").read_text(encoding="utf-8")
        )

    def assert_invalid_case(self, candidate, message: str) -> None:
        with self.assertRaisesRegex(ValidationError, message):
            validate_case(candidate, self.catalog)

    def test_unknown_registered_axes_are_rejected(self) -> None:
        fields = {
            "family": "missing-family",
            "problem": "missing-problem",
            "scheme": "missing-scheme",
            "exact_case": "missing-case",
            "build_profile": "missing-build",
            "launcher": "missing-launcher",
        }
        for field, value in fields.items():
            with self.subTest(field=field):
                self.assert_invalid_case(
                    replace(self.case, **{field: value}), "unknown"
                )

    def test_time_validation(self) -> None:
        candidates = (
            (TimeSpec("0", "0.025", 4), "final_time must be positive"),
            (TimeSpec("0.1", "0", 4), "time_step must be positive"),
            (TimeSpec("0.1", "0.025", 0), "steps must be a positive integer"),
            (TimeSpec("0.1", "0.03", 3), "integral step count"),
            (TimeSpec("NaN", "0.025", 4), "not an exact number"),
        )
        for time_spec, message in candidates:
            with self.subTest(time=time_spec):
                self.assert_invalid_case(replace(self.case, time=time_spec), message)

    def test_degree_mesh_and_parallel_validation(self) -> None:
        candidates = (
            (replace(self.case, test_degree=(3, 4, 4), trial_degree=(3, 3, 3)), "must exceed"),
            (replace(self.case, test_degree=(10, 4, 4)), "maximum 9"),
            (replace(self.case, trial_degree=(0, 3, 3)), "three positive"),
            (replace(self.case, mesh=(0, 4, 4)), "three positive"),
            (replace(self.case, mpi=MpiSpec(2, (1, 1, 1))), "must equal"),
            (
                replace(self.case, mpi=MpiSpec(8, (8, 1, 1))),
                "trial-space DOFs",
            ),
            (replace(self.case, openmp_threads=0), "OpenMP threads"),
            (
                replace(
                    self.case,
                    measurement=replace(self.case.measurement, warmups=-1),
                ),
                "warmups",
            ),
            (
                replace(
                    self.case,
                    measurement=replace(self.case.measurement, samples=0),
                ),
                "samples",
            ),
            (
                replace(
                    self.case,
                    measurement=replace(
                        self.case.measurement, timeout_seconds="0"
                    ),
                ),
                "timeout_seconds",
            ),
            (
                replace(
                    self.case,
                    measurement=replace(
                        self.case.measurement, timeout_seconds="1e400"
                    ),
                ),
                "must not exceed",
            ),
            (replace(self.case, openmp_threads=True), "OpenMP threads"),
            (replace(self.case, mpi=MpiSpec(True, (1, 1, 1))), "MPI ranks"),
        )
        for candidate, message in candidates:
            with self.subTest(message=message):
                self.assert_invalid_case(candidate, message)

        # The library partitions basis-function DOFs, not geometric elements.
        # With four elements and cubic trial splines, five x-ranks are valid.
        validate_case(replace(self.case, mpi=MpiSpec(5, (5, 1, 1))), self.catalog)
        with self.assertRaisesRegex(ValidationError, "direct launcher supports"):
            validate_case(
                replace(
                    self.case,
                    launcher="direct",
                    mpi=MpiSpec(2, (2, 1, 1)),
                ),
                self.catalog,
            )

    def test_time_normalization_is_exact_for_long_and_recurring_values(self) -> None:
        first = Decimal("0.12345678901234567890123456789")
        second = Decimal("0.12345678901234567890123456788")
        self.assertEqual(canonical_decimal(first), str(first))
        self.assertEqual(canonical_decimal(second), str(second))
        self.assertNotEqual(canonical_decimal(first), canonical_decimal(second))

        recurring = copy.deepcopy(self.document)
        recurring["time_discretizations"] = [{"final_time": "0.1", "steps": 3}]
        with tempfile.TemporaryDirectory(prefix="ads-recurring-time-") as temporary:
            path = Path(temporary) / "profile.json"
            path.write_text(json.dumps(recurring), encoding="utf-8")
            plan = Planner(load_profiles(path.parent), self.catalog).plan("smoke")
        self.assertEqual({case.spec.time.time_step for case in plan.cases}, {"1/30"})

        inconsistent = copy.deepcopy(self.document)
        inconsistent["time_discretizations"] = [
            {
                "final_time": "1",
                "steps": 3,
                "time_step": "0.3333333333333333333333333334",
            }
        ]
        with self.assertRaisesRegex(ConfigurationError, "inconsistent"):
            self._parse_document(inconsistent)

    def _parse_document(self, document: dict[str, object]) -> None:
        with tempfile.TemporaryDirectory(prefix="ads-invalid-profile-") as temporary:
            path = Path(temporary) / "profile.json"
            path.write_text(json.dumps(document), encoding="utf-8")
            parse_profile(load_json(path), path)

    def test_profile_rejects_unknown_missing_and_bad_shapes(self) -> None:
        mutations = []
        extra = copy.deepcopy(self.document)
        extra["surprise"] = True
        mutations.append((extra, "unknown surprise"))
        missing = copy.deepcopy(self.document)
        del missing["launcher"]
        mutations.append((missing, "missing launcher"))
        short_vector = copy.deepcopy(self.document)
        short_vector["meshes"] = [[4, 4]]
        mutations.append((short_vector, "exactly three"))
        bool_threads = copy.deepcopy(self.document)
        bool_threads["thread_counts"] = [True]
        mutations.append((bool_threads, "positive integer"))
        bool_schema = copy.deepcopy(self.document)
        bool_schema["schema_version"] = True
        mutations.append((bool_schema, "schema_version"))
        bad_np = copy.deepcopy(self.document)
        bad_np["process_layouts"] = [{"ranks": 2, "grid": [1, 1, 1]}]
        mutations.append((bad_np, "must equal"))
        nonintegral = copy.deepcopy(self.document)
        nonintegral["time_discretizations"] = [
            {"final_time": "0.1", "time_step": "0.03"}
        ]
        mutations.append((nonintegral, "positive integer"))

        for document, message in mutations:
            with self.subTest(message=message):
                if message == "must equal":
                    with tempfile.TemporaryDirectory(
                        prefix="ads-invalid-expanded-"
                    ) as temporary:
                        path = Path(temporary) / "profile.json"
                        path.write_text(json.dumps(document), encoding="utf-8")
                        planner = Planner(load_profiles(path.parent), self.catalog)
                        with self.assertRaisesRegex(ValidationError, message):
                            planner.plan("smoke")
                else:
                    with self.assertRaisesRegex(ConfigurationError, message):
                        self._parse_document(document)

    def test_json_rejects_duplicate_keys_and_nonfinite_numbers(self) -> None:
        invalid_documents = (
            ('{"schema_version": 1, "schema_version": 1}', "duplicate JSON key"),
            ('{"schema_version": NaN}', "non-finite JSON number"),
            ('{"schema_version": Infinity}', "non-finite JSON number"),
        )
        for content, message in invalid_documents:
            with self.subTest(content=content):
                with tempfile.TemporaryDirectory(prefix="ads-invalid-json-") as temporary:
                    path = Path(temporary) / "profile.json"
                    path.write_text(content, encoding="utf-8")
                    with self.assertRaisesRegex(ConfigurationError, message):
                        load_json(path)


if __name__ == "__main__":
    unittest.main(verbosity=2)
