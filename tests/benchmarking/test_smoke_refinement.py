from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from ads_benchmark.smoke_refinement import (
    PROBLEMS,
    SCHEMES,
    SmokeRefinementError,
    verify_refinement,
)


class SmokeRefinementTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory(prefix="ads-smoke-refinement-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def _write_run(
        self,
        name: str,
        *,
        steps: int,
        l2_error: float,
        initial_l2_error: float = 1.0e-14,
        initial_linf_error: float = 2.0e-14,
    ) -> Path:
        run = self.root / name
        for index, (problem, scheme) in enumerate(
            (problem, scheme) for problem in PROBLEMS for scheme in SCHEMES
        ):
            case = run / "cases" / f"case-{index:02d}"
            case.mkdir(parents=True)
            result = {
                "status": "passed",
                "domain_result": {
                    "problem": problem,
                    "scheme": scheme,
                    "steps": steps,
                    "initial_l2_error": initial_l2_error,
                    "initial_linf_error": initial_linf_error,
                    "l2_error": l2_error + index * 1.0e-6,
                },
            }
            (case / "result.json").write_text(
                json.dumps(result), encoding="utf-8"
            )
        return run

    def test_accepts_complete_matrix_with_strict_l2_improvement(self) -> None:
        coarse = self._write_run("coarse", steps=4, l2_error=0.02)
        refined = self._write_run("refined", steps=8, l2_error=0.01)
        rows = verify_refinement(coarse, refined)
        self.assertEqual(len(rows), 9)
        self.assertEqual((rows[0][0], rows[0][1]), ("igrm_l2", "dg"))

    def test_rejects_equal_or_larger_refined_error(self) -> None:
        coarse = self._write_run("coarse", steps=4, l2_error=0.02)
        refined = self._write_run("refined", steps=8, l2_error=0.01)
        target = refined / "cases" / "case-00" / "result.json"
        document = json.loads(target.read_text(encoding="utf-8"))
        document["domain_result"]["l2_error"] = 0.02
        target.write_text(json.dumps(document), encoding="utf-8")
        with self.assertRaisesRegex(SmokeRefinementError, "did not improve"):
            verify_refinement(coarse, refined)

    def test_rejects_initial_projection_above_limit(self) -> None:
        coarse = self._write_run(
            "coarse", steps=4, l2_error=0.02, initial_linf_error=1.1e-10
        )
        refined = self._write_run("refined", steps=8, l2_error=0.01)
        with self.assertRaisesRegex(SmokeRefinementError, "initial projection"):
            verify_refinement(coarse, refined)

    def test_rejects_missing_matrix_member(self) -> None:
        coarse = self._write_run("coarse", steps=4, l2_error=0.02)
        refined = self._write_run("refined", steps=8, l2_error=0.01)
        (refined / "cases" / "case-08" / "result.json").unlink()
        with self.assertRaisesRegex(SmokeRefinementError, "missing smoke results"):
            verify_refinement(coarse, refined)


if __name__ == "__main__":
    unittest.main(verbosity=2)
