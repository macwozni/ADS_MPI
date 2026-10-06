"""Small real-MPI integration gate for the manufactured benchmark adapters.

This is a correctness test, not a performance benchmark: it uses tiny fixed
cases, makes no timing assertion, and writes all runtime artifacts below a
temporary test directory.
"""

from __future__ import annotations

from collections.abc import Mapping
import math
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import tempfile
import unittest

from ads_benchmark.analysis.validation import (
    PARALLEL_ABSOLUTE_TOLERANCE,
    PARALLEL_RELATIVE_TOLERANCE,
)
from ads_benchmark.components.manufactured import ManufacturedTransientAdapter
from ads_benchmark.framework.model import (
    CaseSpec,
    ExecutionContext,
    MeasurementSpec,
    MpiSpec,
    SamplingSpec,
    TimeSpec,
)
from ads_benchmark.validation.fields import (
    RegularGridField,
    compare_fields,
    parse_regular_grid_csv,
    scalar_component,
)

from benchmark_paths import BENCHMARKING_ROOT, REPOSITORY_ROOT


PROBLEMS = ("igrm_l2", "igrm_heat", "pure_diffusion_igrm")
SAMPLE_POINTS = 5


class RealAdapterIntegrationTests(unittest.TestCase):
    """Compile, launch, parse, and validate a bounded seven-case matrix."""

    @classmethod
    def setUpClass(cls) -> None:
        if os.environ.get("SKIP_MPI_CASES", "0") == "1":
            raise unittest.SkipTest("SKIP_MPI_CASES=1")

        cls._temporary = tempfile.TemporaryDirectory(
            prefix="ads-benchmark-adapter-integration."
        )
        cls.addClassCleanup(cls._temporary.cleanup)
        cls._work = Path(cls._temporary.name)
        cls._runs: dict[
            tuple[str, int, int, tuple[int, int, int], bool],
            tuple[Mapping[str, object], RegularGridField | None],
        ] = {}

        config = Path(
            os.environ.get(
                "ADS_BENCHMARK_TEST_CONFIG", str(REPOSITORY_ROOT / "m_options")
            )
        ).resolve()
        if not config.is_file():
            raise RuntimeError(f"benchmark integration config is missing: {config}")

        build = subprocess.run(
            (
                "make",
                "--no-print-directory",
                "-j1",
                "-C",
                str(BENCHMARKING_ROOT),
                f"CONFIG={config}",
                "BUILD=release",
                "build-adapters",
            ),
            cwd=REPOSITORY_ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=300,
            check=False,
        )
        if build.returncode != 0:
            raise RuntimeError(
                "real adapter integration build failed with status "
                f"{build.returncode}:\n{build.stdout[-12000:]}"
            )

        cls._mpi_command = tuple(
            shlex.split(
                os.environ.get(
                    "MPIEXEC", "/opt/lib/mpich-5.0.0/bin/mpiexec"
                )
            )
        )
        cls._mpi_flags = tuple(shlex.split(os.environ.get("MPIEXEC_FLAGS", "")))
        cls._mpi_rank_flag = os.environ.get("MPI_NP_FLAG", "-n")
        if (
            not cls._mpi_command
            or not cls._mpi_rank_flag
            or "\0" in cls._mpi_rank_flag
        ):
            raise RuntimeError("MPI launcher command and rank flag must be nonempty")
        executable = cls._mpi_command[0]
        if os.sep in executable:
            launcher_available = Path(executable).is_file() and os.access(
                executable, os.X_OK
            )
        else:
            launcher_available = shutil.which(executable) is not None
        if not launcher_available:
            raise RuntimeError(f"MPI launcher is unavailable: {executable}")

        try:
            cls._timeout = float(
                os.environ.get("BENCHMARK_INTEGRATION_TIMEOUT", "60")
            )
        except ValueError as error:
            raise RuntimeError(
                "BENCHMARK_INTEGRATION_TIMEOUT must be numeric"
            ) from error
        if not math.isfinite(cls._timeout) or cls._timeout <= 0.0:
            raise RuntimeError(
                "BENCHMARK_INTEGRATION_TIMEOUT must be finite and positive"
            )

    @staticmethod
    def _case_spec(
        problem: str,
        *,
        steps: int,
        ranks: int,
        process_grid: tuple[int, int, int],
        write_samples: bool,
    ) -> CaseSpec:
        return CaseSpec(
            family="temporal",
            problem=problem,
            scheme="dg",
            exact_case="temporal-polynomial",
            time=TimeSpec(
                final_time="0.1",
                time_step={4: "0.025", 8: "0.0125"}[steps],
                steps=steps,
            ),
            mesh=(4, 4, 4),
            test_degree=(4, 4, 4),
            trial_degree=(3, 3, 3),
            mpi=MpiSpec(ranks=ranks, process_grid=process_grid),
            openmp_threads=1,
            sampling=SamplingSpec(
                points_per_axis=SAMPLE_POINTS,
                write_samples=write_samples,
            ),
            measurement=MeasurementSpec(
                warmups=0,
                samples=1,
                timeout_seconds="60",
            ),
            build_profile="release",
            launcher="mpi",
            openmp_dynamic=False,
            openmp_proc_bind="close",
            openmp_places="cores",
        )

    @staticmethod
    def _stop_process_group(
        process: subprocess.Popen[str],
    ) -> tuple[str, str]:
        if process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        try:
            return process.communicate(timeout=5)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            return process.communicate()

    def _run_case(
        self,
        problem: str,
        *,
        steps: int,
        ranks: int,
        process_grid: tuple[int, int, int],
        write_samples: bool,
    ) -> tuple[Mapping[str, object], RegularGridField | None]:
        key = (problem, steps, ranks, process_grid, write_samples)
        cached = self._runs.get(key)
        if cached is not None:
            return cached

        case_directory = Path(
            tempfile.mkdtemp(
                prefix=f"{problem}-n{steps}-r{ranks}.", dir=self._work
            )
        )
        spec = self._case_spec(
            problem,
            steps=steps,
            ranks=ranks,
            process_grid=process_grid,
            write_samples=write_samples,
        )
        adapter = ManufacturedTransientAdapter(name=problem)
        payload = tuple(
            adapter.build_payload_command(
                spec,
                ExecutionContext(
                    repository_root=REPOSITORY_ROOT,
                    case_directory=case_directory,
                ),
            )
        )
        command = (
            *self._mpi_command,
            *self._mpi_flags,
            self._mpi_rank_flag,
            str(ranks),
            *payload,
        )
        environment = os.environ.copy()
        environment.update(
            {
                "OMP_NUM_THREADS": "1",
                "OMP_DYNAMIC": "FALSE",
                "OMP_PROC_BIND": "close",
                "OMP_PLACES": "cores",
            }
        )
        process = subprocess.Popen(
            command,
            cwd=case_directory,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
            start_new_session=True,
        )
        try:
            stdout, stderr = process.communicate(timeout=self._timeout)
        except subprocess.TimeoutExpired:
            stdout, stderr = self._stop_process_group(process)
            self.fail(
                f"{problem}/N={steps}/MPI={ranks} timed out after "
                f"{self._timeout:g}s\nstdout:\n{stdout}\nstderr:\n{stderr}"
            )
        except BaseException:
            self._stop_process_group(process)
            raise

        self.assertEqual(
            process.returncode,
            0,
            f"{problem}/N={steps}/MPI={ranks} failed\n"
            f"stdout:\n{stdout}\nstderr:\n{stderr}",
        )
        result = adapter.parse_result(stdout, stderr)
        adapter.validate_result(spec, result)

        artifact = case_directory / "field_samples.csv"
        field: RegularGridField | None = None
        if write_samples:
            self.assertTrue(artifact.is_file(), f"missing field artifact: {artifact}")
            field = parse_regular_grid_csv(
                artifact.read_text(encoding="utf-8"),
                shape=SAMPLE_POINTS,
                components=("numerical", "exact", "error"),
                label=f"{problem}/N={steps}/MPI={ranks}",
            )
            numerical = field.component_values("numerical")
            exact = field.component_values("exact")
            error = field.component_values("error")
            for actual, expected, reported_error in zip(
                numerical, exact, error, strict=True
            ):
                self.assertAlmostEqual(
                    actual - expected,
                    reported_error,
                    delta=2.0e-15,
                )
        else:
            self.assertFalse(
                artifact.exists(), f"unexpected field artifact: {artifact}"
            )

        completed = (result, field)
        self._runs[key] = completed
        return completed

    def test_all_adapters_execute_and_improve_from_n4_to_n8(self) -> None:
        for problem in PROBLEMS:
            with self.subTest(problem=problem):
                coarse, field = self._run_case(
                    problem,
                    steps=4,
                    ranks=1,
                    process_grid=(1, 1, 1),
                    write_samples=True,
                )
                refined, _ = self._run_case(
                    problem,
                    steps=8,
                    ranks=1,
                    process_grid=(1, 1, 1),
                    write_samples=False,
                )
                self.assertIsNotNone(field)
                self.assertLess(
                    float(refined["l2_error"]), float(coarse["l2_error"])
                )
                self.assertLessEqual(float(coarse["initial_l2_error"]), 1.0e-10)
                self.assertLessEqual(float(refined["initial_l2_error"]), 1.0e-10)

    def test_two_rank_field_matches_serial_field(self) -> None:
        serial_result, serial_field = self._run_case(
            "igrm_l2",
            steps=4,
            ranks=1,
            process_grid=(1, 1, 1),
            write_samples=True,
        )
        parallel_result, parallel_field = self._run_case(
            "igrm_l2",
            steps=4,
            ranks=2,
            process_grid=(2, 1, 1),
            write_samples=True,
        )
        self.assertIsNotNone(serial_field)
        self.assertIsNotNone(parallel_field)
        assert serial_field is not None
        assert parallel_field is not None
        comparison = compare_fields(
            scalar_component(serial_field, "numerical", label="MPI=1"),
            scalar_component(parallel_field, "numerical", label="MPI=2"),
            absolute_tolerance=PARALLEL_ABSOLUTE_TOLERANCE,
            relative_tolerance=PARALLEL_RELATIVE_TOLERANCE,
        )
        self.assertTrue(comparison.passed, comparison.diagnostic())
        self.assertTrue(
            math.isclose(
                float(serial_result["l2_error"]),
                float(parallel_result["l2_error"]),
                rel_tol=PARALLEL_RELATIVE_TOLERANCE,
                abs_tol=PARALLEL_ABSOLUTE_TOLERANCE,
            )
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
