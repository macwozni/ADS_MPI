from __future__ import annotations

from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from benchmark_paths import BENCHMARKING_ROOT, REPOSITORY_ROOT


SOURCE_MAKEFILE = BENCHMARKING_ROOT / "GNUmakefile"


class BenchmarkMakeCleanupTests(unittest.TestCase):
    def test_build_dry_run_refreshes_core_once_before_adapter_compilation(self) -> None:
        completed = subprocess.run(
            [
                "make",
                "--no-print-directory",
                "-j1",
                "-B",
                "-n",
                "-C",
                str(BENCHMARKING_ROOT),
                "build",
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        output = completed.stdout + completed.stderr
        self.assertEqual(completed.returncode, 0, output)
        core_invocation = f"-C '{REPOSITORY_ROOT / 'src'}'"
        self.assertEqual(output.count(core_invocation), 1, output)
        core_position = output.index(core_invocation)
        for adapter in ("igrm_l2", "igrm_heat", "pure_diffusion_igrm"):
            with self.subTest(adapter=adapter):
                object_fragment = f"benchmark_OBJ/{adapter}/benchmark_contract.o"
                self.assertIn(object_fragment, output)
                self.assertLess(core_position, output.index(object_fragment))

    def test_unowned_build_is_unchanged_and_owned_cleanup_preserves_results(self) -> None:
        with tempfile.TemporaryDirectory(prefix="ads-benchmark-make-") as temporary:
            repository = Path(temporary) / "repository"
            benchmarking = repository / "benchmarking"
            benchmarking.mkdir(parents=True)
            shutil.copy2(SOURCE_MAKEFILE, benchmarking / "GNUmakefile")

            build_cache = benchmarking / "build" / "__pycache__" / "foreign.pyc"
            build_cache.parent.mkdir(parents=True)
            build_cache.write_bytes(b"foreign-build-data")
            framework_cache = benchmarking / "package" / "__pycache__" / "cache.pyc"
            framework_cache.parent.mkdir(parents=True)
            framework_cache.write_bytes(b"generated-cache")
            result = repository / "benchmarks" / "legacy" / "result.csv"
            result.parent.mkdir(parents=True)
            result.write_text("user-result\n", encoding="utf-8")

            command = [
                "make",
                "--no-print-directory",
                "-j1",
                "-C",
                str(benchmarking),
                "clean-build",
            ]
            refused = subprocess.run(
                command, check=False, capture_output=True, text=True
            )
            self.assertNotEqual(refused.returncode, 0)
            self.assertIn("Refusing unowned", refused.stdout + refused.stderr)
            self.assertEqual(build_cache.read_bytes(), b"foreign-build-data")
            self.assertEqual(framework_cache.read_bytes(), b"generated-cache")
            self.assertEqual(result.read_text(encoding="utf-8"), "user-result\n")

            marker = benchmarking / "build" / ".ads-benchmark-build"
            marker.write_text(
                f"ADS_MPI_BENCHMARK_BUILD={marker.parent.resolve()}\n",
                encoding="utf-8",
            )
            subprocess.run(command, check=True, capture_output=True, text=True)
            self.assertFalse((benchmarking / "build").exists())
            self.assertFalse(framework_cache.exists())
            self.assertEqual(result.read_text(encoding="utf-8"), "user-result\n")


if __name__ == "__main__":
    unittest.main(verbosity=2)
