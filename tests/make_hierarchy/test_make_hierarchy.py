#!/usr/bin/env python3
"""Hermetic regression tests for the hierarchical GNU Make interface.

The production tree is copied to a temporary directory and built with this
file acting as a deterministic compiler, archiver, and MPI launcher.  No real
Fortran compiler, MPI installation, external library, or repository artifact
directory is used by this suite.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import time
import unittest


SCRIPT_PATH = Path(__file__).resolve()
REPOSITORY_ROOT = SCRIPT_PATH.parents[2]

CORE_SOURCES = [
    "Setup.F90",
    "Interfaces.F90",
    "argument_parser.F90",
    "knot_vector.F90",
    "parallelism.F90",
    "communicators.F90",
    "gauss.F90",
    "basis.F90",
    "math.F90",
    "reorderRHS.F90",
    "my_mpi.F90",
    "sparse.F90",
    "RHS_eq.F90",
    "igrm_space.F90",
    "operator_assembly.F90",
    "rhs_assembly.F90",
    "solution_reconstruction.F90",
    "projection_engine.F90",
    "ads_lifecycle.F90",
    "mumps_solver.F90",
    "ads_directional_solve.F90",
    "utils.F90",
    "plot.F90",
    "gnuplot.F90",
    "vtk.F90",
    "solution_output.F90",
    "ADS.F90",
    "time_scheme.F90",
]

CORE_MODULES = [
    "setup",
    "interfaces",
    "argument_parser",
    "knot_vector",
    "parallelism",
    "communicators",
    "gauss",
    "basis",
    "math",
    "reorderrhs",
    "my_mpi",
    "sparse",
    "rhs_eq",
    "igrm_space",
    "operator_assembly",
    "rhs_assembly",
    "solution_reconstruction",
    "projection_engine",
    "ads_lifecycle",
    "mumps_solver",
    "ads_directional_solve",
    "utils",
    "plot",
    "gnuplot",
    "vtk",
    "solution_output",
    "adss",
    "time_scheme",
]

# These files remain in src as historical references and are deliberately not
# part of the active library.  Keeping the exception list explicit makes a new
# unregistered production source fail this regression test.
INACTIVE_CORE_SOURCES = ["_sparse.F90", "debug.F90"]

PROBLEMS = {
    "l2": ("L2Projection", ["input_data.F90", "RHS_fun.F90", "main.F90"]),
    "heat": ("heat", ["input_data.F90", "RHS_fun.F90", "main.F90"]),
    "eriksson": ("eriksson", ["input_data.F90", "RHS_fun.F90", "main.F90"]),
    "pure_diffusion_igrm": (
        "pure_diffusion_igrm",
        ["input_data.F90", "RHS_fun.F90", "main.F90"],
    ),
    "oil": ("oil", ["input_data.F90", "RHS_fun.F90", "main.F90"]),
    "igrm_l2": (
        "igrm_L2Projection",
        ["input_data.F90", "RHS_fun.F90", "main.F90"],
    ),
    "igrm_heat": (
        "igrm_heat",
        ["input_data.F90", "RHS_fun.F90", "main.F90"],
    ),
    "igrm_eirksson": (
        "igrm_eirksson",
        ["input_data.F90", "RHS_fun.F90", "igrm_mumps_solver.F90", "main.F90"],
    ),
    "igrm_stokes": (
        "igrm_stokes",
        ["input_data.F90", "RHS_fun.F90", "igrm_stokes_solver.F90", "main.F90"],
    ),
    "igrm_pollution": (
        "igrm_pollution",
        ["input_data.F90", "RHS_fun.F90", "pollution_dpg_solver.F90", "main.F90"],
    ),
}


def _append_json(path_variable: str, payload: dict[str, object]) -> None:
    path_text = os.environ.get(path_variable)
    if not path_text:
        return
    with Path(path_text).open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(payload, sort_keys=True) + "\n")


def _module_output_directory(arguments: list[str], output: Path) -> Path:
    for index, argument in enumerate(arguments):
        if argument.startswith("-J") and len(argument) > 2:
            return Path(argument[2:])
        if argument == "-J" and index + 1 < len(arguments):
            return Path(arguments[index + 1])
        if argument == "-module" and index + 1 < len(arguments):
            return Path(arguments[index + 1])
    return output.parent


def _fake_compiler(arguments: list[str]) -> int:
    try:
        output = Path(arguments[arguments.index("-o") + 1])
    except (ValueError, IndexError):
        print("fake compiler: missing -o OUTPUT", file=sys.stderr)
        return 2

    output.parent.mkdir(parents=True, exist_ok=True)
    if "-c" in arguments:
        try:
            source = Path(arguments[arguments.index("-c") + 1])
        except IndexError:
            print("fake compiler: missing source after -c", file=sys.stderr)
            return 2
        output.write_text(f"object={source.resolve()}\n", encoding="utf-8")
        module_directory = _module_output_directory(arguments, output)
        module_directory.mkdir(parents=True, exist_ok=True)
        module_pattern = re.compile(
            r"^\s*module\s+(?!procedure\b|function\b|subroutine\b)([A-Za-z]\w*)",
            re.IGNORECASE | re.MULTILINE,
        )
        for module_name in module_pattern.findall(source.read_text(encoding="utf-8")):
            (module_directory / f"{module_name.lower()}.mod").write_text(
                f"module={module_name.lower()}\nsource={source.resolve()}\n",
                encoding="utf-8",
            )
        kind = "compile"
        source_text: str | None = str(source.resolve())
    else:
        program = """#!/usr/bin/env python3
import json
import os
from pathlib import Path
import sys

payload = {
    "argv": sys.argv[1:],
    "cwd": os.getcwd(),
    "environment": {
        name: os.environ.get(name)
        for name in (
            "ADS_OIL_RANDOM_SEED",
            "FORWARDED",
            "OMP_DYNAMIC",
            "OMP_NUM_THREADS",
            "OMP_PROC_BIND",
        )
    },
}
with Path(os.environ["FAKE_RUN_LOG"]).open("a", encoding="utf-8") as stream:
    stream.write(json.dumps(payload, sort_keys=True) + "\\n")
"""
        output.write_text(program, encoding="utf-8")
        output.chmod(output.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        kind = "link"
        source_text = None

    _append_json(
        "FAKE_TOOL_LOG",
        {
            "argv": arguments,
            "kind": kind,
            "output": str(output.resolve()),
            "source": source_text,
            "tool": "compiler",
        },
    )
    return 0


def _fake_archiver(arguments: list[str]) -> int:
    if len(arguments) < 2:
        print("fake archiver: expected FLAGS ARCHIVE [OBJECT ...]", file=sys.stderr)
        return 2
    archive = Path(arguments[1])
    archive.parent.mkdir(parents=True, exist_ok=True)
    archive.write_text(
        "\n".join(str(Path(item).resolve()) for item in arguments[2:]) + "\n",
        encoding="utf-8",
    )
    _append_json(
        "FAKE_TOOL_LOG",
        {
            "argv": arguments,
            "kind": "archive",
            "output": str(archive.resolve()),
            "tool": "archiver",
        },
    )
    return 0


def _fake_mpiexec(arguments: list[str]) -> int:
    _append_json(
        "FAKE_TOOL_LOG",
        {"argv": arguments, "kind": "launch", "tool": "mpiexec"},
    )
    for index, argument in enumerate(arguments):
        candidate = Path(argument)
        if candidate.is_file() and os.access(candidate, os.X_OK):
            completed = subprocess.run(
                [str(candidate), *arguments[index + 1 :]],
                check=False,
                env=os.environ.copy(),
            )
            return completed.returncode
    print("fake mpiexec: executable argument not found", file=sys.stderr)
    return 2


def _dispatch_fake_tool() -> int | None:
    tool_name = Path(sys.argv[0]).name
    if tool_name == "fake-ff":
        return _fake_compiler(sys.argv[1:])
    if tool_name == "fake-ar":
        return _fake_archiver(sys.argv[1:])
    if tool_name == "fake-mpiexec":
        return _fake_mpiexec(sys.argv[1:])
    return None


class MakeFixture:
    """A private, disposable repository with deterministic build tools."""

    def __init__(self) -> None:
        self._temporary_directory = tempfile.TemporaryDirectory(
            prefix="ads-make-hierarchy-"
        )
        self.root = Path(self._temporary_directory.name) / "repository"
        self.root.mkdir()
        self.build_root = self.root / "artifacts"
        self.config = self.root / "fake-config.mk"
        self.tool_log = self.root / "fake-tools.jsonl"
        self.run_log = self.root / "fake-runs.jsonl"
        self._copy_repository_inputs()
        self._install_fake_tools()
        self._write_configuration()
        self.environment = os.environ.copy()
        # A suite launched by recursive make receives every command-line
        # override as both a make flag and an exported environment variable.
        # Remove those inputs so the disposable fixture is governed solely by
        # its fake configuration and by overrides made explicitly in a test.
        for variable in (
            "ARGS",
            "BUILD",
            "BUILD_ROOT",
            "COMPILER",
            "CONFIG",
            "CORE_OBJ_DIR",
            "EXEC",
            "EXEC_DIR",
            "FC",
            "FF",
            "FFLAGS",
            "GNUMAKEFLAGS",
            "MAKEFLAGS",
            "MAKELEVEL",
            "MAKEOVERRIDES",
            "MFLAGS",
            "MPIEXEC",
            "MPIEXEC_FLAGS",
            "MPIFC",
            "OIL_SEED",
            "OMP_DYNAMIC",
            "OMP_NUM_THREADS",
            "OMP_PROC_BIND",
            "PROBLEM",
            "PROBLEM_OBJ_DIR",
            "RUN_DIR",
            "RUN_ENV",
            "SOURCE_ALL",
            "USER_LIB",
        ):
            self.environment.pop(variable, None)
        self.environment.update(
            {
                "FAKE_RUN_LOG": str(self.run_log),
                "FAKE_TOOL_LOG": str(self.tool_log),
                "LC_ALL": "C",
            }
        )

    def close(self) -> None:
        self._temporary_directory.cleanup()

    def _copy_repository_inputs(self) -> None:
        shutil.copy2(REPOSITORY_ROOT / "Makefile", self.root / "Makefile")

        source_destination = self.root / "src"
        source_destination.mkdir()
        for name in ("GNUmakefile", "sources.mk"):
            shutil.copy2(REPOSITORY_ROOT / "src" / name, source_destination / name)
        for source_path in sorted((REPOSITORY_ROOT / "src").glob("*.F90")):
            shutil.copy2(source_path, source_destination / source_path.name)

        problems_destination = self.root / "problems"
        problems_destination.mkdir()
        for name in ("GNUmakefile", "problems.mk", "problem.mk"):
            shutil.copy2(
                REPOSITORY_ROOT / "problems" / name,
                problems_destination / name,
            )
        problem_directories = sorted(
            path
            for path in (REPOSITORY_ROOT / "problems").iterdir()
            if path.is_dir() and (path / "GNUmakefile").is_file()
        )
        for source_directory in problem_directories:
            directory_name = source_directory.name
            destination = problems_destination / directory_name
            destination.mkdir()
            shutil.copy2(
                source_directory / "GNUmakefile",
                destination / "GNUmakefile",
            )
            for source_path in sorted(source_directory.glob("*.F90")):
                shutil.copy2(
                    source_path,
                    destination / source_path.name,
                )

        mymake_destination = self.root / "mymake"
        mymake_destination.mkdir()
        for name in ("makefile", "m_files", "legacy-source-build.mk"):
            shutil.copy2(REPOSITORY_ROOT / "mymake" / name, mymake_destination / name)

        # Root `clean` delegates test cleanup.  A tiny owned fixture keeps that
        # orchestration hermetic; the real tests hierarchy is validated by the
        # outer repository before this suite is launched.
        tests_destination = self.root / "tests"
        tests_destination.mkdir()
        (tests_destination / "GNUmakefile").write_text(
            ".PHONY: clean\nclean:\n\t$(RM) -- generated-test-artifact\n",
            encoding="utf-8",
        )

    def _install_fake_tools(self) -> None:
        tools_directory = self.root / "fake-tools"
        tools_directory.mkdir()
        for tool_name in ("fake-ff", "fake-ar", "fake-mpiexec"):
            destination = tools_directory / tool_name
            shutil.copy2(SCRIPT_PATH, destination)
            destination.chmod(
                destination.stat().st_mode
                | stat.S_IXUSR
                | stat.S_IXGRP
                | stat.S_IXOTH
            )
        self.fake_compiler = tools_directory / "fake-ff"
        self.fake_archiver = tools_directory / "fake-ar"
        self.fake_mpiexec = tools_directory / "fake-mpiexec"

    def _write_configuration(self) -> None:
        self.config.write_text(
            f"""COMPILER = {self.fake_compiler}
MODULE_OUTPUT = -J.
BUILD ?= debug
ifeq ($(BUILD),debug)
BUILD_OPTS = --profile=debug
else ifeq ($(BUILD),release)
BUILD_OPTS = --profile=release
else
$(error unsupported test BUILD='$(BUILD)')
endif
FF = $(COMPILER) $(BUILD_OPTS)
FFLAGS =
AR = {self.fake_archiver}
USER_LIB =
MPIEXEC = {self.fake_mpiexec}
MPIFC = $(COMPILER)
PYTHON = {sys.executable}
PFUNIT_ROOT = unused
MUMPS_DIR = unused
DOXYGEN = false
SUITE_TIMEOUT = 30s
DRIVER_CLI_TIMEOUT = 1s
DRIVER_SMOKE_TIMEOUT = 1s
DRIVER_INTEGRATION_TIMEOUT = 1
SKIP_MPI_CASES = 1
""",
            encoding="utf-8",
        )

    def make(
        self,
        *targets: str,
        directory: Path | None = None,
        variables: dict[str, str | Path | int] | None = None,
        expect_success: bool = True,
    ) -> subprocess.CompletedProcess[str]:
        effective_variables: dict[str, str | Path | int] = {
            "CONFIG": self.config,
            "BUILD_ROOT": self.build_root,
        }
        if variables:
            effective_variables.update(variables)
        command = [
            "make",
            "--no-print-directory",
            "-j1",
            "-C",
            str(directory or self.root),
            *(f"{name}={value}" for name, value in effective_variables.items()),
            *targets,
        ]
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            encoding="utf-8",
            env=self.environment,
            timeout=30,
        )
        if expect_success and completed.returncode != 0:
            raise AssertionError(
                f"command failed ({completed.returncode}): {' '.join(command)}\n"
                f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
            )
        if not expect_success and completed.returncode == 0:
            raise AssertionError(
                f"command unexpectedly succeeded: {' '.join(command)}\n"
                f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
            )
        return completed

    def clear_logs(self) -> None:
        self.tool_log.unlink(missing_ok=True)
        self.run_log.unlink(missing_ok=True)

    def tool_records(self) -> list[dict[str, object]]:
        return self._json_lines(self.tool_log)

    def run_records(self) -> list[dict[str, object]]:
        return self._json_lines(self.run_log)

    @staticmethod
    def _json_lines(path: Path) -> list[dict[str, object]]:
        if not path.exists():
            return []
        return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]


class HierarchicalMakeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.fixture = MakeFixture()

    def tearDown(self) -> None:
        self.fixture.close()

    def assert_failed_with(
        self, result: subprocess.CompletedProcess[str], message: str
    ) -> None:
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(message, result.stdout + result.stderr)

    @staticmethod
    def _records_of_kind(
        records: list[dict[str, object]], kind: str
    ) -> list[dict[str, object]]:
        return [record for record in records if record.get("kind") == kind]

    def test_manifests_incremental_recovery_full_build_run_and_clean(self) -> None:
        copied_core_sources = sorted(
            path.name for path in (self.fixture.root / "src").glob("*.F90")
        )
        self.assertEqual(
            copied_core_sources,
            sorted([*CORE_SOURCES, *INACTIVE_CORE_SOURCES]),
        )
        listed_sources = self.fixture.make(
            "-s", "list-sources", directory=self.fixture.root / "src"
        ).stdout.splitlines()
        self.assertEqual(listed_sources, CORE_SOURCES)

        listed_problems = self.fixture.make("-s", "list-problems").stdout.splitlines()
        self.assertEqual(listed_problems, list(PROBLEMS))
        copied_problem_directories = sorted(
            path.name
            for path in (self.fixture.root / "problems").iterdir()
            if path.is_dir() and (path / "GNUmakefile").is_file()
        )
        self.assertEqual(
            copied_problem_directories,
            sorted(directory_name for directory_name, _ in PROBLEMS.values()),
        )
        directory_lines = self.fixture.make(
            "-s", "list-directories", directory=self.fixture.root / "problems"
        ).stdout.splitlines()
        self.assertEqual(
            dict(line.split() for line in directory_lines),
            {name: description[0] for name, description in PROBLEMS.items()},
        )
        for problem_name, (directory_name, expected_sources) in PROBLEMS.items():
            self.assertTrue(
                (self.fixture.root / "problems" / directory_name / "GNUmakefile").is_file()
            )
            copied_problem_sources = sorted(
                path.name
                for path in (self.fixture.root / "problems" / directory_name).glob(
                    "*.F90"
                )
            )
            self.assertEqual(
                copied_problem_sources,
                sorted(expected_sources),
                problem_name,
            )
            actual_sources = self.fixture.make(
                "-s",
                "list-sources",
                directory=self.fixture.root / "problems",
                variables={"PROBLEM": problem_name},
            ).stdout.splitlines()
            self.assertEqual(actual_sources, expected_sources, problem_name)

        self.fixture.clear_logs()
        self.fixture.make("build", variables={"PROBLEM": "heat"})
        records = self.fixture.tool_records()
        compiles = self._records_of_kind(records, "compile")
        self.assertEqual(len(compiles), len(CORE_SOURCES) + 3)
        self.assertEqual(len(self._records_of_kind(records, "archive")), 1)
        self.assertEqual(len(self._records_of_kind(records, "link")), 1)
        self.assertEqual(
            [Path(str(record["source"])).name for record in compiles[: len(CORE_SOURCES)]],
            CORE_SOURCES,
        )

        self.fixture.clear_logs()
        self.fixture.make("build", variables={"PROBLEM": "heat"})
        self.assertEqual(self.fixture.tool_records(), [], "no-op build ran a tool")

        self.fixture.clear_logs()
        self.fixture.make(
            "build", variables={"PROBLEM": "heat", "BUILD": "release"}
        )
        records = self.fixture.tool_records()
        self.assertEqual(len(self._records_of_kind(records, "compile")), 31)
        self.assertEqual(len(self._records_of_kind(records, "archive")), 1)
        self.assertEqual(len(self._records_of_kind(records, "link")), 1)
        for record in self._records_of_kind(records, "compile"):
            self.assertIn("--profile=release", record["argv"])

        self.fixture.clear_logs()
        self.fixture.make(
            "build", variables={"PROBLEM": "heat", "BUILD": "release"}
        )
        self.assertEqual(self.fixture.tool_records(), [])

        time.sleep(0.02)
        (self.fixture.build_root / "setup.mod").unlink()
        self.fixture.clear_logs()
        self.fixture.make(
            "build", variables={"PROBLEM": "heat", "BUILD": "release"}
        )
        records = self.fixture.tool_records()
        self.assertEqual(len(self._records_of_kind(records, "compile")), 31)
        self.assertEqual(len(self._records_of_kind(records, "archive")), 1)
        self.assertEqual(len(self._records_of_kind(records, "link")), 1)
        self.assertTrue((self.fixture.build_root / "setup.mod").is_file())

        heat_objects = self.fixture.root / "problems" / "heat" / "_OBJ"
        time.sleep(0.02)
        (heat_objects / "input_data.mod").unlink()
        self.fixture.clear_logs()
        self.fixture.make(
            "build", variables={"PROBLEM": "heat", "BUILD": "release"}
        )
        records = self.fixture.tool_records()
        self.assertEqual(len(self._records_of_kind(records, "compile")), 3)
        self.assertEqual(len(self._records_of_kind(records, "archive")), 0)
        self.assertEqual(len(self._records_of_kind(records, "link")), 1)
        self.assertTrue((heat_objects / "input_data.mod").is_file())

        self.fixture.clear_logs()
        self.fixture.make("all", variables={"BUILD": "release"})
        records = self.fixture.tool_records()
        expected_new_problem_sources = sum(
            len(sources) for _, sources in PROBLEMS.values()
        ) - len(PROBLEMS["heat"][1])
        self.assertEqual(
            len(self._records_of_kind(records, "compile")),
            expected_new_problem_sources,
        )
        self.assertEqual(len(self._records_of_kind(records, "archive")), 0)
        self.assertEqual(len(self._records_of_kind(records, "link")), 9)
        for problem_name in PROBLEMS:
            self.assertTrue(
                (self.fixture.build_root / "EXEC" / problem_name).is_file(),
                problem_name,
            )

        self.fixture.clear_logs()
        self.fixture.make("all", variables={"BUILD": "release"})
        self.assertEqual(self.fixture.tool_records(), [], "no-op all ran a tool")

        heat_executable = self.fixture.build_root / "EXEC" / "heat"
        heat_executable.unlink()
        self.fixture.clear_logs()
        self.fixture.make(
            "build", variables={"PROBLEM": "heat", "BUILD": "release"}
        )
        records = self.fixture.tool_records()
        self.assertEqual(len(self._records_of_kind(records, "compile")), 0)
        self.assertEqual(len(self._records_of_kind(records, "archive")), 0)
        self.assertEqual(len(self._records_of_kind(records, "link")), 1)

        time.sleep(0.02)
        (self.fixture.build_root / "LIB" / "libads.a").touch()
        self.fixture.clear_logs()
        self.fixture.make(
            "build", variables={"PROBLEM": "heat", "BUILD": "release"}
        )
        records = self.fixture.tool_records()
        self.assertEqual(len(self._records_of_kind(records, "compile")), 3)
        self.assertEqual(len(self._records_of_kind(records, "archive")), 0)
        self.assertEqual(len(self._records_of_kind(records, "link")), 1)

        self.fixture.clear_logs()
        self.fixture.make("all", variables={"BUILD": "release"})
        records = self.fixture.tool_records()
        self.assertEqual(
            len(self._records_of_kind(records, "compile")),
            expected_new_problem_sources,
        )
        self.assertEqual(len(self._records_of_kind(records, "archive")), 0)
        self.assertEqual(len(self._records_of_kind(records, "link")), 9)

        runtime_variables = {
            "ARGS": "arg-one arg-two",
            "BUILD": "release",
            "MPIEXEC_FLAGS": "--fake-option",
            "MPI_NP_FLAG": "--ranks",
            "NP": 3,
            "OIL_SEED": 4242,
            "OMP_DYNAMIC": "TRUE",
            "OMP_NUM_THREADS": 6,
            "OMP_PROC_BIND": "spread",
            "PROBLEM": "oil",
            "RUN_DIR": "custom-run",
            "RUN_ENV": "FORWARDED=visible",
        }
        shown = self.fixture.make("show-run", variables=runtime_variables).stdout
        for expected_line in (
            r"(?m)^PROBLEM\s+oil$",
            r"(?m)^ARGS\s+arg-one arg-two$",
            r"(?m)^NP\s+3$",
            r"(?m)^MPI_NP_FLAG\s+--ranks$",
            r"(?m)^OMP_NUM_THREADS\s+6$",
            r"(?m)^OMP_DYNAMIC\s+TRUE$",
            r"(?m)^OMP_PROC_BIND\s+spread$",
            r"(?m)^OIL_SEED\s+4242$",
            r"(?m)^PROBLEM_RUN_ENV\s+ADS_OIL_RANDOM_SEED=4242$",
        ):
            self.assertRegex(shown, expected_line)
        self.assertRegex(
            shown,
            rf"(?m)^RUN_DIR\s+{re.escape(str(self.fixture.root / 'custom-run'))}$",
        )

        self.fixture.clear_logs()
        self.fixture.make("run", variables=runtime_variables)
        tool_records = self.fixture.tool_records()
        self.assertEqual(self._records_of_kind(tool_records, "compile"), [])
        launches = self._records_of_kind(tool_records, "launch")
        self.assertEqual(len(launches), 1)
        launch_arguments = launches[0]["argv"]
        self.assertEqual(launch_arguments[0:3], ["--fake-option", "--ranks", "3"])
        self.assertEqual(Path(launch_arguments[3]).name, "oil")
        self.assertEqual(launch_arguments[4:], ["arg-one", "arg-two"])
        run_records = self.fixture.run_records()
        self.assertEqual(len(run_records), 1)
        self.assertEqual(run_records[0]["argv"], ["arg-one", "arg-two"])
        self.assertEqual(run_records[0]["cwd"], str(self.fixture.root / "custom-run"))
        self.assertEqual(
            run_records[0]["environment"],
            {
                "ADS_OIL_RANDOM_SEED": "4242",
                "FORWARDED": "visible",
                "OMP_DYNAMIC": "TRUE",
                "OMP_NUM_THREADS": "6",
                "OMP_PROC_BIND": "spread",
            },
        )

        core_sentinel = self.fixture.build_root / "_OBJ" / "foreign.keep"
        executable_sentinel = self.fixture.build_root / "EXEC" / "foreign.keep"
        module_sentinel = self.fixture.build_root / "foreign.mod"
        problem_sentinel = heat_objects / "foreign.keep"
        for sentinel in (
            core_sentinel,
            executable_sentinel,
            module_sentinel,
            problem_sentinel,
        ):
            sentinel.parent.mkdir(parents=True, exist_ok=True)
            sentinel.write_text("foreign\n", encoding="utf-8")

        documentation_artifact = self.fixture.root / "doxygen" / "owned"
        documentation_artifact.parent.mkdir()
        documentation_artifact.write_text("generated\n", encoding="utf-8")
        test_artifact = self.fixture.root / "tests" / "generated-test-artifact"
        test_artifact.write_text("generated\n", encoding="utf-8")

        self.fixture.make("clean", variables={"BUILD": "release"})
        self.fixture.make("clean", variables={"BUILD": "release"})
        for sentinel in (
            core_sentinel,
            executable_sentinel,
            module_sentinel,
            problem_sentinel,
        ):
            self.assertEqual(sentinel.read_text(encoding="utf-8"), "foreign\n")
        self.assertTrue(self.fixture.config.is_file())
        self.assertFalse(documentation_artifact.exists())
        self.assertFalse(test_artifact.exists())
        self.assertFalse((self.fixture.build_root / "LIB" / "libads.a").exists())
        for source_name in CORE_SOURCES:
            self.assertFalse(
                (self.fixture.build_root / "_OBJ" / Path(source_name).with_suffix(".o")).exists()
            )
        for module_name in CORE_MODULES:
            self.assertFalse((self.fixture.build_root / f"{module_name}.mod").exists())
        for problem_name, (directory_name, sources) in PROBLEMS.items():
            self.assertFalse((self.fixture.build_root / "EXEC" / problem_name).exists())
            object_directory = self.fixture.root / "problems" / directory_name / "_OBJ"
            for source_name in sources:
                self.assertFalse((object_directory / Path(source_name).with_suffix(".o")).exists())

    def test_unsafe_path_guards_fail_before_invoking_tools(self) -> None:
        checks = [
            (
                self.fixture.root / "src",
                {"BUILD_ROOT": self.fixture.root},
                "Unsafe BUILD_ROOT",
            ),
            (
                self.fixture.root / "src",
                {"CORE_OBJ_DIR": self.fixture.root / "src"},
                "Unsafe CORE_OBJ_DIR",
            ),
            (
                self.fixture.root / "problems" / "heat",
                {"BUILD_ROOT": self.fixture.root},
                "Unsafe BUILD_ROOT",
            ),
            (
                self.fixture.root / "problems" / "heat",
                {"PROBLEM_OBJ_DIR": self.fixture.root / "problems" / "heat"},
                "Unsafe PROBLEM_OBJ_DIR",
            ),
            (
                self.fixture.root / "problems" / "heat",
                {"EXEC_DIR": self.fixture.root / "problems"},
                "Unsafe EXEC_DIR",
            ),
        ]
        self.fixture.clear_logs()
        for directory, variables, expected_message in checks:
            with self.subTest(message=expected_message, variables=variables):
                result = self.fixture.make(
                    "build",
                    directory=directory,
                    variables=variables,
                    expect_success=False,
                )
                self.assert_failed_with(result, expected_message)

        legacy_source = self.fixture.root / "src" / "Setup.F90"
        result = self.fixture.make(
            "build",
            directory=self.fixture.root / "mymake",
            variables={
                "BUILD_ROOT": self.fixture.root,
                "EXEC": "legacy_safe",
                "SOURCE_ALL": legacy_source,
            },
            expect_success=False,
        )
        self.assert_failed_with(result, "Unsafe legacy BUILD_ROOT")
        result = self.fixture.make(
            "build",
            directory=self.fixture.root / "mymake",
            variables={"EXEC": "../escape", "SOURCE_ALL": legacy_source},
            expect_success=False,
        )
        self.assert_failed_with(result, "Unsafe legacy EXEC")
        self.assertEqual(self.fixture.tool_records(), [])

    def test_legacy_build_validation_and_atomic_marker_cleanup(self) -> None:
        source_all = " ".join(
            str(self.fixture.root / "src" / source_name)
            for source_name in ("Setup.F90", "Interfaces.F90")
        )
        legacy_root = self.fixture.build_root / "_LEGACY_OBJ"
        legacy_root.mkdir(parents=True)
        unowned = legacy_root / "foreign.keep"
        unowned.write_text("foreign\n", encoding="utf-8")
        result = self.fixture.make(
            "build",
            directory=self.fixture.root / "mymake",
            variables={"EXEC": "legacy_a", "SOURCE_ALL": source_all},
            expect_success=False,
        )
        self.assert_failed_with(result, "Refusing to adopt non-empty unmarked directory")
        self.assertEqual(unowned.read_text(encoding="utf-8"), "foreign\n")
        shutil.rmtree(legacy_root)

        self.fixture.clear_logs()
        self.fixture.make(
            "build",
            directory=self.fixture.root / "mymake",
            variables={"EXEC": "legacy_a", "SOURCE_ALL": source_all},
        )
        records = self.fixture.tool_records()
        self.assertEqual(len(self._records_of_kind(records, "compile")), 2)
        self.assertEqual(len(self._records_of_kind(records, "link")), 1)
        self.assertEqual(len(self._records_of_kind(records, "archive")), 0)
        legacy_a_directory = legacy_root / "legacy_a"
        legacy_a_executable = self.fixture.build_root / "EXEC" / "legacy_a"
        self.assertTrue(legacy_a_directory.is_dir())
        self.assertTrue(legacy_a_executable.is_file())

        self.fixture.clear_logs()
        self.fixture.make(
            "build",
            directory=self.fixture.root / "mymake",
            variables={"EXEC": "legacy_a", "SOURCE_ALL": source_all},
        )
        self.assertEqual(self.fixture.tool_records(), [])

        inherited_manifest = self.fixture.make(
            "-s",
            "list-problems",
            directory=self.fixture.root / "mymake",
            variables={"EXEC": "legacy_a", "SOURCE_ALL": source_all},
        ).stdout.splitlines()
        self.assertEqual(inherited_manifest, list(PROBLEMS))

        invalid_cases: list[tuple[dict[str, str | Path], str]] = [
            ({"EXEC": "empty", "SOURCE_ALL": ""}, "SOURCE_ALL was supplied but is empty"),
            (
                {
                    "EXEC": "lowercase",
                    "SOURCE_ALL": self.fixture.root / "unsupported.f90",
                },
                "ordered .F90 files only",
            ),
        ]
        (self.fixture.root / "unsupported.f90").write_text(
            "program unsupported\nend program unsupported\n", encoding="utf-8"
        )
        duplicate_one = self.fixture.root / "duplicate-one" / "same.F90"
        duplicate_two = self.fixture.root / "duplicate-two" / "same.F90"
        duplicate_one.parent.mkdir()
        duplicate_two.parent.mkdir()
        duplicate_one.write_text("module same_one\nend module\n", encoding="utf-8")
        duplicate_two.write_text("module same_two\nend module\n", encoding="utf-8")
        invalid_cases.append(
            (
                {
                    "EXEC": "duplicates",
                    "SOURCE_ALL": f"{duplicate_one} {duplicate_two}",
                },
                "duplicate basenames",
            )
        )
        for variables, expected_message in invalid_cases:
            with self.subTest(legacy_validation=expected_message):
                result = self.fixture.make(
                    "build",
                    directory=self.fixture.root / "mymake",
                    variables=variables,
                    expect_success=False,
                )
                self.assert_failed_with(result, expected_message)

        self.fixture.make(
            "build",
            directory=self.fixture.root / "mymake",
            variables={"EXEC": "legacy_z", "SOURCE_ALL": source_all},
        )
        legacy_z_directory = legacy_root / "legacy_z"
        legacy_z_executable = self.fixture.build_root / "EXEC" / "legacy_z"
        corrupt_marker = legacy_z_directory / ".ads-legacy-object-dir"
        corrupt_marker.write_text("corrupt\n", encoding="utf-8")

        def snapshot(directory: Path) -> dict[str, bytes]:
            return {
                str(path.relative_to(directory)): path.read_bytes()
                for path in directory.rglob("*")
                if path.is_file()
            }

        before = snapshot(legacy_root)
        result = self.fixture.make(
            "legacy-clean-all",
            directory=self.fixture.root / "mymake",
            expect_success=False,
        )
        self.assert_failed_with(result, "Refusing unmarked legacy child")
        self.assertEqual(snapshot(legacy_root), before)
        self.assertTrue(legacy_a_executable.is_file())
        self.assertTrue(legacy_z_executable.is_file())

        corrupt_marker.write_text(
            "ADS_MPI_LEGACY_EXEC=legacy_z\n"
            f"ADS_MPI_LEGACY_DIR={legacy_z_directory}\n",
            encoding="utf-8",
        )
        executable_sentinel = self.fixture.build_root / "EXEC" / "foreign.keep"
        executable_sentinel.write_text("foreign\n", encoding="utf-8")
        self.fixture.make(
            "legacy-clean-all", directory=self.fixture.root / "mymake"
        )
        self.assertFalse(legacy_root.exists())
        self.assertFalse(legacy_a_executable.exists())
        self.assertFalse(legacy_z_executable.exists())
        self.assertEqual(executable_sentinel.read_text(encoding="utf-8"), "foreign\n")
        self.assertTrue(
            (self.fixture.build_root / "EXEC" / ".ads-legacy-exec-dir").is_file()
        )
        self.fixture.make(
            "legacy-clean-all", directory=self.fixture.root / "mymake"
        )
        self.assertEqual(executable_sentinel.read_text(encoding="utf-8"), "foreign\n")


def main() -> int:
    fake_tool_result = _dispatch_fake_tool()
    if fake_tool_result is not None:
        return fake_tool_result
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(HierarchicalMakeTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
