from __future__ import annotations

import unittest

from benchmark_paths import BENCHMARKING_ROOT


class BenchmarkRepositoryLayoutTests(unittest.TestCase):
    def test_benchmarking_tree_contains_no_automated_tests(self) -> None:
        misplaced: list[str] = []
        for path in BENCHMARKING_ROOT.rglob("*"):
            relative = path.relative_to(BENCHMARKING_ROOT)
            if relative.parts and relative.parts[0] == "build":
                continue
            if path.is_dir() and path.name in {"tests", "selftests"}:
                misplaced.append(relative.as_posix())
            if path.is_file() and path.suffix == ".py" and (
                path.name.startswith("test_") or path.name.endswith("_test.py")
            ):
                misplaced.append(relative.as_posix())
        self.assertEqual(
            misplaced,
            [],
            "automated tests belong under tests/benchmarking, not benchmarking/",
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
