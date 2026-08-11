#!/usr/bin/env python3
"""End-to-end numerical checks for the real MPI/OpenMP problem drivers."""

from __future__ import annotations

from dataclasses import dataclass
import math
import os
from pathlib import Path
import re
import shlex
import signal
import subprocess
import tempfile
import xml.etree.ElementTree as ET


EXPECTED_EXTENT = (0, 30, 0, 30, 0, 30)
EXPECTED_POINT_COUNT = 31**3
VTK_ABS_TOL = 1.0e-9
VTK_REL_TOL = 1.0e-9
SCALAR_ABS_TOL = 1.0e-12
SCALAR_REL_TOL = 1.0e-9
ZERO_SOLUTION_TOL = 1.0e-10


@dataclass(frozen=True)
class DriverCase:
    label: str
    directory: Path
    output: str
    returncode: int

    @property
    def ok(self) -> bool:
        return self.returncode == 0


@dataclass(frozen=True)
class VTIData:
    extent: tuple[int, ...]
    origin: tuple[float, ...]
    spacing: tuple[float, ...]
    values: tuple[float, ...]


class IntegrationSuite:
    def __init__(self) -> None:
        mpi_command = shlex.split(
            os.environ.get("MPIEXEC", "/opt/lib/mpich-5.0.0/bin/mpiexec")
        )
        if not mpi_command:
            raise ValueError("MPIEXEC must name an MPI launcher")
        self.mpi_command = mpi_command
        self.exec_dir = Path(os.environ.get("EXEC_DIR", "../../mymake/EXEC")).resolve()
        self.timeout = float(os.environ.get("DRIVER_INTEGRATION_TIMEOUT", "90"))
        self.work = tempfile.TemporaryDirectory(prefix="mpi-ads-driver-integration.")
        self.work_dir = Path(self.work.name)
        self.checks = 0
        self.failures = 0
        self._vti_cache: dict[Path, VTIData] = {}

    def close(self) -> None:
        self.work.cleanup()

    def report(self, passed: bool, label: str, details: str = "") -> bool:
        self.checks += 1
        if passed:
            print(f"PASS {label}")
            return True
        self.failures += 1
        print(f"FAIL {label}{': ' + details if details else ''}")
        return False

    @staticmethod
    def _safe_case_name(label: str) -> str:
        return re.sub(r"[^A-Za-z0-9_.-]+", "-", label).strip("-")

    @staticmethod
    def _asan_options(environment: dict[str, str]) -> str:
        current = environment.get("ASAN_OPTIONS", "")
        return f"{current}:detect_leaks=0" if current else "detect_leaks=0"

    @staticmethod
    def _stop_process_group(process: subprocess.Popen[str]) -> str:
        """Stop mpiexec and every rank, then collect the remaining output."""
        if process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        try:
            output, _ = process.communicate(timeout=5)
            return output
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            output, _ = process.communicate()
            return output

    def run_driver(
        self,
        label: str,
        ranks: int,
        threads: int,
        driver: str,
        arguments: list[object],
        extra_environment: dict[str, str] | None = None,
    ) -> DriverCase:
        case_dir = self.work_dir / self._safe_case_name(label)
        case_dir.mkdir()
        executable = self.exec_dir / driver
        command = [
            *self.mpi_command,
            "-n",
            str(ranks),
            str(executable),
            *(str(argument) for argument in arguments),
        ]
        environment = os.environ.copy()
        environment.update(
            {
                "OMP_DYNAMIC": "FALSE",
                "OMP_NUM_THREADS": str(threads),
                "OMP_PROC_BIND": "close",
            }
        )
        if extra_environment:
            environment.update(extra_environment)
        environment["ASAN_OPTIONS"] = self._asan_options(environment)

        process = subprocess.Popen(
            command,
            cwd=case_dir,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            start_new_session=True,
        )
        timed_out = False
        try:
            output, _ = process.communicate(timeout=self.timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            output = self._stop_process_group(process)
        except BaseException:
            self._stop_process_group(process)
            raise

        returncode = 124 if timed_out else process.returncode
        (case_dir / "driver.log").write_text(output, encoding="utf-8")
        case = DriverCase(label, case_dir, output, returncode)
        if timed_out:
            details = f"timeout after {self.timeout:g}s"
        elif returncode != 0:
            details = f"status {returncode}"
        elif "ERROR: AddressSanitizer" in output or "Fortran runtime error" in output:
            details = "runtime diagnostic in output"
        else:
            details = ""
        passed = returncode == 0 and not details
        if not self.report(passed, f"{label} run", details) and output:
            print("  " + output.rstrip().replace("\n", "\n  "))
        return case

    def expect_no_l2_failure(self, case: DriverCase) -> None:
        self.report(
            case.ok and "not OK" not in case.output,
            f"{case.label} constant-projection oracle",
            "driver printed 'not OK' despite returning success"
            if "not OK" in case.output
            else "driver did not complete",
        )

    def expect_iterations(self, case: DriverCase, expected: list[int]) -> None:
        actual = [
            int(line)
            for line in case.output.splitlines()
            if re.fullmatch(r"\s*[+-]?\d+\s*", line)
        ]
        self.report(
            case.ok and actual == expected,
            f"{case.label} physical iterations",
            f"expected {expected}, got {actual}",
        )

    def expect_zero_solution(self, case: DriverCase, expected: list[int]) -> None:
        measurements: list[tuple[int, float]] = []
        pattern = re.compile(r"^\s*solution max abs iter\s+(\d+):\s+(\S+)\s*$")
        for line in case.output.splitlines():
            match = pattern.match(line)
            if match:
                value = match.group(2).replace("D", "E").replace("d", "e")
                try:
                    measurements.append((int(match.group(1)), float(value)))
                except ValueError:
                    measurements.append((int(match.group(1)), math.nan))
        iterations = [iteration for iteration, _ in measurements]
        values = [value for _, value in measurements]
        passed = (
            case.ok
            and iterations == expected
            and all(
                math.isfinite(value) and abs(value) <= ZERO_SOLUTION_TOL
                for value in values
            )
        )
        self.report(
            passed,
            f"{case.label} zero-solution oracle",
            f"expected finite global maxima <= {ZERO_SOLUTION_TOL:g} at "
            f"{expected}, got {measurements}",
        )

    @staticmethod
    def _tag(element: ET.Element) -> str:
        return element.tag.rsplit("}", 1)[-1]

    @classmethod
    def _first_tag(cls, root: ET.Element, tag: str) -> ET.Element:
        for element in root.iter():
            if cls._tag(element) == tag:
                return element
        raise ValueError(f"missing {tag} element")

    def read_vti(self, path: Path) -> VTIData:
        path = path.resolve()
        if path in self._vti_cache:
            return self._vti_cache[path]
        if not path.is_file():
            raise ValueError("file does not exist")

        try:
            root = ET.parse(path).getroot()
        except ET.ParseError as error:
            raise ValueError(f"invalid XML: {error}") from error
        if self._tag(root) != "VTKFile" or root.get("type") != "ImageData":
            raise ValueError("root is not a VTK ImageData file")

        image = self._first_tag(root, "ImageData")
        try:
            extent = tuple(int(item) for item in image.attrib["WholeExtent"].split())
            origin = tuple(float(item) for item in image.attrib["Origin"].split())
            spacing = tuple(float(item) for item in image.attrib["Spacing"].split())
        except (KeyError, ValueError) as error:
            raise ValueError(f"invalid ImageData metadata: {error}") from error
        if extent != EXPECTED_EXTENT:
            raise ValueError(f"extent is {extent}, expected {EXPECTED_EXTENT}")
        if len(origin) != 3 or len(spacing) != 3:
            raise ValueError("origin and spacing must have three components")
        if not all(math.isfinite(item) for item in (*origin, *spacing)):
            raise ValueError("origin or spacing contains NaN/Inf")

        result = None
        for element in root.iter():
            if self._tag(element) == "DataArray" and element.get("Name") == "Result":
                result = element
                break
        if result is None:
            raise ValueError("missing Result DataArray")
        if result.get("type") != "Float64" or result.get("format") != "ascii":
            raise ValueError("Result must be an ASCII Float64 DataArray")
        try:
            values = tuple(float(item) for item in (result.text or "").split())
        except ValueError as error:
            raise ValueError(f"invalid Result value: {error}") from error
        if len(values) != EXPECTED_POINT_COUNT:
            raise ValueError(
                f"Result contains {len(values)} values, expected {EXPECTED_POINT_COUNT}"
            )
        if not all(math.isfinite(item) for item in values):
            raise ValueError("Result contains NaN/Inf")
        if not any(item != 0.0 for item in values):
            raise ValueError("Result is identically zero")

        data = VTIData(extent, origin, spacing, values)
        self._vti_cache[path] = data
        return data

    def expect_valid_vti(self, case: DriverCase, filename: str) -> VTIData | None:
        path = case.directory / filename
        try:
            data = self.read_vti(path)
        except ValueError as error:
            self.report(False, f"{case.label} {filename} VTI", str(error))
            return None
        self.report(True, f"{case.label} {filename} VTI")
        return data

    def expect_matching_vti(
        self, reference: DriverCase, candidate: DriverCase, filename: str
    ) -> None:
        try:
            left = self.read_vti(reference.directory / filename)
            right = self.read_vti(candidate.directory / filename)
        except ValueError as error:
            self.report(False, f"{candidate.label} {filename} matches serial", str(error))
            return
        if left.extent != right.extent or left.origin != right.origin or left.spacing != right.spacing:
            self.report(
                False,
                f"{candidate.label} {filename} matches serial",
                "VTK grid metadata differs",
            )
            return
        for index, (expected, actual) in enumerate(zip(left.values, right.values)):
            if not math.isclose(
                actual, expected, rel_tol=VTK_REL_TOL, abs_tol=VTK_ABS_TOL
            ):
                self.report(
                    False,
                    f"{candidate.label} {filename} matches serial",
                    f"value {index}: expected {expected:.12g}, got {actual:.12g}",
                )
                return
        self.report(True, f"{candidate.label} {filename} matches serial")

    def expect_solution_changed(self, case: DriverCase) -> None:
        try:
            initial = self.read_vti(case.directory / "step0.vti")
            final = self.read_vti(case.directory / "step1.vti")
        except ValueError as error:
            self.report(False, f"{case.label} physical solution changed", str(error))
            return
        changed = any(
            not math.isclose(before, after, rel_tol=1.0e-12, abs_tol=1.0e-12)
            for before, after in zip(initial.values, final.values)
        )
        self.report(
            changed,
            f"{case.label} physical solution changed",
            "step1 is identical to the initialization at step0",
        )

    def expect_sample_norm_decreased(self, case: DriverCase) -> None:
        try:
            initial = self.read_vti(case.directory / "step0.vti")
            final = self.read_vti(case.directory / "step1.vti")
        except ValueError as error:
            self.report(False, f"{case.label} dissipative heat oracle", str(error))
            return
        initial_norm = math.sqrt(sum(value * value for value in initial.values))
        final_norm = math.sqrt(sum(value * value for value in final.values))
        self.report(
            final_norm < initial_norm,
            f"{case.label} dissipative heat oracle",
            f"sample norm did not decrease ({initial_norm:.8g} -> {final_norm:.8g})",
        )

    @staticmethod
    def final_scalar(case: DriverCase) -> float:
        for line in reversed(case.output.splitlines()):
            text = line.strip().replace("D", "E").replace("d", "e")
            if not text:
                continue
            try:
                return float(text)
            except ValueError:
                continue
        raise ValueError("no numeric result in driver output")

    def expect_oil_result(
        self, reference: DriverCase, candidate: DriverCase | None = None
    ) -> None:
        label = reference.label if candidate is None else candidate.label
        try:
            expected = self.final_scalar(reference)
            actual = expected if candidate is None else self.final_scalar(candidate)
        except ValueError as error:
            self.report(False, f"{label} drained result", str(error))
            return
        if candidate is None:
            passed = math.isfinite(actual) and actual > 0.0
            details = f"expected a finite positive result, got {actual!r}"
        else:
            passed = math.isfinite(actual) and math.isclose(
                actual,
                expected,
                rel_tol=SCALAR_REL_TOL,
                abs_tol=SCALAR_ABS_TOL,
            )
            details = f"serial {expected:.16g}, candidate {actual:.16g}"
        self.report(passed, f"{label} drained result", details)


def run_suite(suite: IntegrationSuite) -> None:
    # Constant projection: real MUMPS, a 2x2x2 MPI grid, uneven DOF ownership,
    # and degree-four support crossing rank boundaries in every direction.
    l2_serial = suite.run_driver(
        "L2 serial p4", 1, 1, "l2", [3, 3, 3, 4, 1, 1, 1]
    )
    suite.expect_no_l2_failure(l2_serial)
    l2_hybrid = suite.run_driver(
        "L2 hybrid p4 2x2x2", 8, 4, "l2", [3, 3, 3, 4, 2, 2, 2]
    )
    suite.expect_no_l2_failure(l2_hybrid)

    # A real transient plus distributed gather/broadcast and VTK output.
    heat_serial = suite.run_driver(
        "heat serial", 1, 1, "heat", [3, 1, 1, 0.01, 1, 1, 1]
    )
    suite.expect_iterations(heat_serial, [0, 1])
    suite.expect_valid_vti(heat_serial, "step0.vti")
    suite.expect_valid_vti(heat_serial, "step1.vti")
    suite.expect_solution_changed(heat_serial)
    heat_hybrid = suite.run_driver(
        "heat hybrid 2x2x2", 8, 4, "heat", [3, 1, 1, 0.01, 2, 2, 2]
    )
    suite.expect_iterations(heat_hybrid, [0, 1])
    suite.expect_matching_vti(heat_serial, heat_hybrid, "step0.vti")
    suite.expect_matching_vti(heat_serial, heat_hybrid, "step1.vti")

    eriksson_serial = suite.run_driver(
        "Eriksson serial", 1, 1, "eriksson", [2, 1, 1, 0.01, 1, 1, 1]
    )
    suite.expect_iterations(eriksson_serial, [0, 1])
    suite.expect_valid_vti(eriksson_serial, "step0.vti")
    suite.expect_valid_vti(eriksson_serial, "step1.vti")
    suite.expect_solution_changed(eriksson_serial)
    eriksson_hybrid = suite.run_driver(
        "Eriksson hybrid", 2, 4, "eriksson", [2, 1, 1, 0.01, 2, 1, 1]
    )
    suite.expect_iterations(eriksson_hybrid, [0, 1])
    suite.expect_matching_vti(eriksson_serial, eriksson_hybrid, "step0.vti")
    suite.expect_matching_vti(eriksson_serial, eriksson_hybrid, "step1.vti")

    # The pure-diffusion driver has no result file, so require every real
    # scheme to complete both initialization and the t>0 step on MPI ranks.
    for scheme in ("dg", "pr", "be"):
        pure = suite.run_driver(
            f"pure diffusion {scheme.upper()}",
            2,
            4,
            "pure_diffusion_igrm",
            [3, 1, 2, 1, 1, 1, 0.1, scheme],
        )
        suite.expect_iterations(pure, [0, 1])
        suite.expect_zero_solution(pure, [0, 1])

    # All public iGRM scheme dispatches with a serial/distributed result oracle.
    for scheme in ("dg", "pr", "be"):
        serial = suite.run_driver(
            f"iGRM L2 {scheme.upper()} serial",
            1,
            1,
            "igrm_l2",
            [2, 2, 2, 3, 3, 3, 1, 1, 1, 1, 1, 1, scheme],
        )
        suite.expect_valid_vti(serial, "step0.vti")
        hybrid = suite.run_driver(
            f"iGRM L2 {scheme.upper()} hybrid",
            2,
            4,
            "igrm_l2",
            [2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 1, 1, scheme],
        )
        suite.expect_matching_vti(serial, hybrid, "step0.vti")

    for scheme in ("dg", "pr", "be"):
        serial = suite.run_driver(
            f"iGRM heat {scheme.upper()} serial",
            1,
            1,
            "igrm_heat",
            [2, 2, 2, 3, 3, 3, 2, 2, 2, 1, 1, 1, 1, 0.001, scheme],
        )
        suite.expect_iterations(serial, [0, 1])
        suite.expect_valid_vti(serial, "step0.vti")
        suite.expect_valid_vti(serial, "step1.vti")
        suite.expect_solution_changed(serial)
        if scheme == "dg":
            suite.expect_sample_norm_decreased(serial)
        hybrid = suite.run_driver(
            f"iGRM heat {scheme.upper()} hybrid",
            2,
            4,
            "igrm_heat",
            [2, 2, 2, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 0.001, scheme],
        )
        suite.expect_iterations(hybrid, [0, 1])
        suite.expect_matching_vti(serial, hybrid, "step0.vti")
        suite.expect_matching_vti(serial, hybrid, "step1.vti")

    # The opt-in seed makes the random oil geometry identical without changing
    # normal runs. This turns the former status-only check into a race/MPI oracle.
    oil_environment = {"ADS_OIL_RANDOM_SEED": "20260811"}
    oil_tail = [1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5]
    oil_serial = suite.run_driver(
        "oil serial OMP1",
        1,
        1,
        "oil",
        [3, 2, 1, 1, 1, 1, 0.05, *oil_tail],
        oil_environment,
    )
    suite.expect_iterations(oil_serial, [0, 1])
    suite.expect_oil_result(oil_serial)
    oil_openmp = suite.run_driver(
        "oil serial OMP4",
        1,
        4,
        "oil",
        [3, 2, 1, 1, 1, 1, 0.05, *oil_tail],
        oil_environment,
    )
    suite.expect_iterations(oil_openmp, [0, 1])
    suite.expect_oil_result(oil_serial, oil_openmp)
    oil_hybrid = suite.run_driver(
        "oil hybrid",
        2,
        2,
        "oil",
        [3, 2, 2, 1, 1, 1, 0.05, *oil_tail],
        oil_environment,
    )
    suite.expect_iterations(oil_hybrid, [0, 1])
    suite.expect_oil_result(oil_serial, oil_hybrid)


def main() -> int:
    def exit_on_termination(signum: int, _frame: object) -> None:
        raise SystemExit(128 + signum)

    signal.signal(signal.SIGTERM, exit_on_termination)
    suite = IntegrationSuite()
    try:
        run_suite(suite)
        if suite.failures:
            print(f"FAILED ({suite.failures}/{suite.checks} integration checks)")
            return 1
        print(f"OK ({suite.checks} integration checks)")
        return 0
    finally:
        suite.close()


if __name__ == "__main__":
    raise SystemExit(main())
