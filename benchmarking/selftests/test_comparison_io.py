from __future__ import annotations

import csv
import io
import json
from pathlib import Path
import tempfile
import unittest

from ads_benchmark.comparison import ComparisonReport
from ads_benchmark.comparison_io import (
    comparison_csv_text,
    comparison_json_text,
    write_comparison_file,
)
from ads_benchmark.framework.errors import StorageError


class ComparisonIoTests(unittest.TestCase):
    def _report(self) -> ComparisonReport:
        return ComparisonReport(
            status="no-regression-detected",
            baseline={"run_id": "a"},
            candidate={"run_id": "b"},
            compatibility={"passed": True, "reasons": []},
            policy={"regression_threshold": 0.05},
            summary={"case_count": 1},
            cases=(
                {
                    "case_id": "case-one",
                    "numerical": {"passed": True},
                    "timing": {
                        "status": "no-regression-detected",
                        "baseline": {
                            "count": 5,
                            "median_seconds": 1.0,
                            "median_absolute_deviation_seconds": 0.1,
                            "range_seconds": 0.4,
                        },
                        "candidate": {
                            "count": 5,
                            "median_seconds": 1.02,
                            "median_absolute_deviation_seconds": 0.11,
                            "range_seconds": 0.5,
                        },
                        "median_ratio_candidate_over_baseline": 1.02,
                    },
                },
            ),
            diagnostics=(),
        )

    def test_json_and_csv_are_deterministic_and_complete(self) -> None:
        report = self._report()
        document = json.loads(comparison_json_text(report))
        self.assertEqual(document["schema_version"], 1)
        self.assertEqual(document["status"], "no-regression-detected")
        rows = list(csv.DictReader(io.StringIO(comparison_csv_text(report))))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["case_id"], "case-one")
        self.assertEqual(rows[0]["baseline_median_seconds"], "1.0")
        self.assertEqual(rows[0]["candidate_over_baseline_median_ratio"], "1.02")

    def test_explicit_output_is_atomic_and_rejects_symlink(self) -> None:
        with tempfile.TemporaryDirectory(prefix="ads-comparison-io-") as temporary:
            root = Path(temporary)
            output = root / "comparison.json"
            write_comparison_file(output, "first\n")
            write_comparison_file(output, "second\n")
            self.assertEqual(output.read_text(encoding="utf-8"), "second\n")
            outside = root / "outside"
            outside.write_text("keep\n", encoding="utf-8")
            output.unlink()
            output.symlink_to(outside)
            with self.assertRaisesRegex(StorageError, "not a regular"):
                write_comparison_file(output, "bad\n")
            self.assertEqual(outside.read_text(encoding="utf-8"), "keep\n")


if __name__ == "__main__":
    unittest.main(verbosity=2)
