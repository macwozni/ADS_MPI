#!/usr/bin/env python3
"""Unit tests for the LCOV denominator and threshold validator."""

from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
import io
import json
from pathlib import Path
import stat
import tempfile
import unittest

import check_coverage


class CoverageCheckerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="ads-coverage-check-")
        self.root = Path(self.temporary.name)
        self.source_root = self.root / "src"
        self.source_root.mkdir()
        self.sources = [self.source_root / "one.F90", self.source_root / "two.F90"]
        for source in self.sources:
            source.write_text(f"module {source.stem}\nend module\n", encoding="utf-8")
        self.manifest = self.source_root / "sources.mk"
        self.manifest.write_text(
            "SRC_FILES := \\\n\tone.F90 \\\n\ttwo.F90\n", encoding="utf-8"
        )
        self.tracefile = self.root / "coverage.info"
        self.summary = self.root / "summary.json"
        self.object_root = self.root / "objects"
        self.object_root.mkdir()
        self.fake_gcov = self.root / "fake-gcov"
        self.fake_gcov.write_text(
            "#!/bin/sh\nprintf '%s\\n' 'No executable lines'\n", encoding="utf-8"
        )
        self.fake_gcov.chmod(
            self.fake_gcov.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def write_tracefile(self, include_second: bool = True, include_extra: bool = False) -> None:
        records = [
            f"SF:{self.sources[0]}\n"
            "FN:10,11,one_fn\nFNDA:2,one_fn\n"
            "DA:10,2\nDA:11,0\n"
            "BRDA:10,0,0,2\nBRDA:10,0,1,0\nend_of_record\n"
        ]
        if include_second:
            records.append(
                f"SF:{self.sources[1]}\n"
                "FN:20,two_fn\nFNDA:0,two_fn\n"
                "DA:20,0\nDA:21,0\n"
                "BRDA:20,0,0,-\nBRDA:20,0,1,0\nend_of_record\n"
            )
        if include_extra:
            extra = self.source_root / "extra.F90"
            extra.write_text("module extra\nend module\n", encoding="utf-8")
            records.append(f"SF:{extra}\nFN:1,extra\nFNDA:1,extra\nDA:1,1\nBRDA:1,0,0,1\nend_of_record\n")
        self.tracefile.write_text("".join(records), encoding="utf-8")

    def arguments(self, **thresholds: float) -> list[str]:
        values = {"lines": 20.0, "functions": 40.0, "branches": 20.0}
        values.update(thresholds)
        return [
            "--tracefile", str(self.tracefile),
            "--source-manifest", str(self.manifest),
            "--source-root", str(self.source_root),
            "--repository-root", str(self.root),
            "--object-root", str(self.object_root),
            "--gcov", str(self.fake_gcov),
            "--summary-json", str(self.summary),
            "--min-lines", str(values["lines"]),
            "--min-functions", str(values["functions"]),
            "--min-branches", str(values["branches"]),
        ]

    def invoke(self, arguments: list[str]) -> tuple[int, str, str]:
        stdout = io.StringIO()
        stderr = io.StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            status = check_coverage.main(arguments)
        return status, stdout.getvalue(), stderr.getvalue()

    def test_complete_tracefile_passes_and_writes_machine_readable_summary(self) -> None:
        self.write_tracefile()
        status, stdout, stderr = self.invoke(self.arguments())
        self.assertEqual(status, 0, stderr)
        self.assertIn("2 active src files", stdout)
        payload = json.loads(self.summary.read_text(encoding="utf-8"))
        self.assertTrue(payload["passed"])
        self.assertEqual(payload["totals"]["lines"]["hit"], 1)
        self.assertEqual(payload["totals"]["lines"]["found"], 4)
        self.assertEqual(payload["totals"]["functions"]["percent"], 50.0)
        self.assertEqual(payload["totals"]["branches"]["percent"], 25.0)
        self.assertEqual(payload["instrumentable_source_files"], 2)
        self.assertEqual(payload["noninstrumentable_sources"], [])
        self.assertEqual(sorted(payload["files"]), ["src/one.F90", "src/two.F90"])

    def test_each_threshold_is_enforced_independently(self) -> None:
        self.write_tracefile()
        for name, thresholds in (
            ("lines", {"lines": 25.01}),
            ("functions", {"functions": 50.01}),
            ("branches", {"branches": 25.01}),
        ):
            with self.subTest(metric=name):
                status, stdout, stderr = self.invoke(self.arguments(**thresholds))
                self.assertEqual(status, 1)
                self.assertIn("[FAIL]", stdout)
                self.assertIn("thresholds were not met", stderr)

    def test_lcov_2_function_ids_are_supported(self) -> None:
        self.tracefile.write_text(
            f"SF:{self.sources[0]}\n"
            "FNL:0,10,11\nFNA:0,2,one_fn\nFNA:0,1,one_generated_fn\n"
            "DA:10,2\nBRDA:10,0,0,2\nend_of_record\n"
            f"SF:{self.sources[1]}\n"
            "FNL:0,20,21\nFNA:0,0,two_fn\n"
            "DA:20,0\nBRDA:20,0,0,0\nend_of_record\n",
            encoding="utf-8",
        )
        status, _, stderr = self.invoke(self.arguments())
        self.assertEqual(status, 0, stderr)
        payload = json.loads(self.summary.read_text(encoding="utf-8"))
        self.assertEqual(payload["totals"]["functions"]["hit"], 2)
        self.assertEqual(payload["totals"]["functions"]["found"], 3)
        self.assertAlmostEqual(payload["totals"]["functions"]["percent"], 66.666667)

    def test_missing_active_source_is_rejected_before_thresholds(self) -> None:
        self.write_tracefile(include_second=False)
        status, _, stderr = self.invoke(self.arguments(lines=0, functions=0, branches=0))
        self.assertEqual(status, 2)
        self.assertIn("active sources missing", stderr)
        self.assertIn("src/two.F90", stderr)

    def test_absent_source_is_accepted_only_with_gcov_no_code_proof(self) -> None:
        self.write_tracefile(include_second=False)
        (self.object_root / "two.gcno").write_text("metadata\n", encoding="utf-8")
        status, stdout, stderr = self.invoke(
            self.arguments(lines=0, functions=0, branches=0)
        )
        self.assertEqual(status, 0, stderr)
        self.assertIn("2 active src files (1 instrumentable)", stdout)
        payload = json.loads(self.summary.read_text(encoding="utf-8"))
        self.assertEqual(payload["noninstrumentable_sources"], ["src/two.F90"])
        self.assertFalse(payload["files"]["src/two.F90"]["instrumentable"])

    def test_undeclared_source_is_rejected(self) -> None:
        self.write_tracefile(include_extra=True)
        status, _, stderr = self.invoke(self.arguments(lines=0, functions=0, branches=0))
        self.assertEqual(status, 2)
        self.assertIn("undeclared sources present", stderr)
        self.assertIn("src/extra.F90", stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
