#!/usr/bin/env python3
"""Verify one-to-one primary tests for library and problem source files."""

from pathlib import Path
import os
import subprocess
import sys


REPO = Path(__file__).resolve().parent.parent
RUNNER = Path(__file__).resolve().parent / "GNUmakefile"
SOURCE_RUNNER = REPO / "src"
PROBLEM_RUNNER = REPO / "problems"
SOURCE_TEST_MAP = Path(__file__).resolve().parent / "test-map.tsv"
PROBLEM_TEST_MAP = Path(__file__).resolve().parent / "problem-test-map.tsv"
GROUP_RUNNERS = {
    "src": RUNNER.parent / "src" / "GNUmakefile",
    "problems": RUNNER.parent / "problems" / "GNUmakefile",
    "driver": RUNNER.parent / "driver" / "GNUmakefile",
    "build": RUNNER.parent / "build" / "GNUmakefile",
}


def active_sources():
    """Ask the src build owner for its public list of production sources."""
    result = subprocess.run(
        [
            "make",
            "--no-print-directory",
            "-s",
            "-C",
            str(SOURCE_RUNNER),
            "list-sources",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    sources = []
    for raw_line in result.stdout.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        path = Path(line)
        if path.is_absolute():
            try:
                line = path.relative_to(REPO).as_posix()
            except ValueError as error:
                raise ValueError(
                    f"src list-sources returned a path outside the repository: {line}"
                ) from error
        elif line.startswith("../src/"):
            line = line.removeprefix("../")
        elif not line.startswith("src/"):
            line = f"src/{line}"
        sources.append(line)
    return sources


def active_problem_sources():
    """Return every manifested problem source that is not a program driver."""
    directory_result = subprocess.run(
        [
            "make",
            "--no-print-directory",
            "-s",
            "-C",
            str(PROBLEM_RUNNER),
            "list-directories",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    directories = {}
    for raw_line in directory_result.stdout.splitlines():
        fields = raw_line.split()
        if len(fields) != 2:
            raise ValueError(
                "problems list-directories returned an invalid line: " + raw_line
            )
        problem, directory = fields
        directories[problem] = directory

    problem_result = subprocess.run(
        [
            "make",
            "--no-print-directory",
            "-s",
            "-C",
            str(PROBLEM_RUNNER),
            "list-problems",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    problems = [line.strip() for line in problem_result.stdout.splitlines() if line.strip()]
    if problems != list(directories):
        raise ValueError("problem names and directory mappings are inconsistent")

    sources = []
    for problem in problems:
        source_result = subprocess.run(
            [
                "make",
                "--no-print-directory",
                "-s",
                "-C",
                str(PROBLEM_RUNNER),
                f"PROBLEM={problem}",
                "list-sources",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        source_names = [
            line.strip() for line in source_result.stdout.splitlines() if line.strip()
        ]
        if sum(Path(name).name == "main.F90" for name in source_names) != 1:
            raise ValueError(f"problem {problem} must manifest exactly one main.F90")
        problem_directory = PROBLEM_RUNNER / directories[problem]
        for source_name in source_names:
            source_path = Path(source_name)
            candidate = (
                source_path if source_path.is_absolute() else problem_directory / source_path
            )
            try:
                relative = candidate.resolve().relative_to(REPO.resolve()).as_posix()
            except ValueError as error:
                raise ValueError(
                    f"problem {problem} returned a source outside the repository: "
                    f"{source_name}"
                ) from error
            if Path(relative).name != "main.F90":
                sources.append(relative)
    return sources


def declared_pairs(test_map):
    pairs = []
    for line_number, raw_line in enumerate(
        test_map.read_text(encoding="utf-8").splitlines(), start=1
    ):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 2:
            raise ValueError(f"{test_map}:{line_number}: expected source and test path")
        pairs.append(tuple(fields))
    return pairs


def duplicates(values):
    seen = set()
    repeated = set()
    for value in values:
        if value in seen:
            repeated.add(value)
        seen.add(value)
    return sorted(repeated)


def make_list(path, variable):
    """Read a simple, backslash-continued make variable assignment."""
    entries = []
    collecting = False
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not collecting:
            prefix = f"{variable} :="
            if not line.startswith(prefix):
                continue
            line = line[len(prefix) :].strip()
            collecting = True

        continued = line.endswith("\\")
        if continued:
            line = line[:-1].strip()
        entries.extend(line.split())
        if not continued:
            break
    return entries


def primary_tests(suite_dir):
    """Return authored primary-test candidates, ignoring generated pFUnit F90."""
    pf_tests = sorted(suite_dir.glob("test_*.pf"))
    generated_stems = {path.stem for path in pf_tests}
    fortran_tests = sorted(
        path
        for path in suite_dir.glob("test_*.F90")
        if path.stem not in generated_stems
    )
    return pf_tests + fortran_tests


def expanded_make_database(suite_dir):
    """Return GNU Make's expanded database without building the suite."""
    result = subprocess.run(
        [
            "make",
            "--no-print-directory",
            "-pn",
            "-C",
            str(suite_dir),
            "--eval=.PHONY: __layout_database__",
            "--eval=__layout_database__:",
            "__layout_database__",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout


def validate_primary_mapping(active, pairs, map_name, errors):
    """Validate one source/test map and return its ordered suite names."""
    mapped_sources = [source for source, _ in pairs]
    mapped_tests = [test for _, test in pairs]
    mapped_suites = []

    for source in sorted(set(active) - set(mapped_sources)):
        errors.append(f"active {map_name} source has no primary test: {source}")
    for source in sorted(set(mapped_sources) - set(active)):
        errors.append(f"{map_name} test map contains inactive source: {source}")
    for source in duplicates(mapped_sources):
        errors.append(f"{map_name} source is mapped more than once: {source}")
    for test in duplicates(mapped_tests):
        errors.append(f"{map_name} test file is mapped more than once: {test}")

    for source, test in pairs:
        source_path = REPO / source
        test_path = REPO / test
        test_parts = Path(test).parts

        if not source_path.is_file():
            errors.append(f"mapped source does not exist: {source}")
        if not test_path.is_file():
            errors.append(f"mapped test does not exist: {test}")
        if not test.startswith("tests/"):
            errors.append(f"primary test is outside tests/: {test}")
            continue
        if len(test_parts) != 3:
            errors.append(f"primary test must be directly inside its suite: {test}")
            continue

        suite = test_parts[1]
        mapped_suites.append(suite)
        suite_dir = REPO / "tests" / suite
        makefile = suite_dir / "GNUmakefile"
        if not makefile.is_file():
            errors.append(f"suite has no GNUmakefile: {suite}")
            continue

        make_text = makefile.read_text(encoding="utf-8")
        try:
            make_database = expanded_make_database(suite_dir)
        except subprocess.CalledProcessError as error:
            errors.append(f"cannot expand suite GNUmakefile for {suite}: {error}")
            continue
        source_reference = Path(
            os.path.relpath(source_path.resolve(), suite_dir.resolve())
        ).as_posix()
        if source_reference not in make_database:
            errors.append(
                f"suite GNUmakefile does not reference the mapped source path "
                f"{source}: {suite}"
            )
        if "run:" not in make_text:
            errors.append(f"suite has no run target: {suite}")
        if test_path.stem not in make_text:
            errors.append(f"suite GNUmakefile does not reference primary test: {test}")

        candidates = [path.relative_to(REPO).as_posix() for path in primary_tests(suite_dir)]
        if candidates != [test]:
            errors.append(
                f"suite {suite} must have exactly the mapped primary test; found: "
                + (", ".join(candidates) if candidates else "none")
            )

    return mapped_suites, mapped_tests


def main():
    errors = []
    try:
        sources = active_sources()
        problem_sources = active_problem_sources()
        source_pairs = declared_pairs(SOURCE_TEST_MAP)
        problem_pairs = declared_pairs(PROBLEM_TEST_MAP)
    except (OSError, subprocess.CalledProcessError, ValueError) as error:
        print(
            f"ERROR: cannot obtain or parse the active test manifests: {error}",
            file=sys.stderr,
        )
        return 1
    mapped_suites, mapped_tests = validate_primary_mapping(
        sources, source_pairs, "src", errors
    )
    problem_mapped_suites, problem_mapped_tests = validate_primary_mapping(
        problem_sources, problem_pairs, "problem", errors
    )
    for test in duplicates(mapped_tests + problem_mapped_tests):
        errors.append(f"primary test is mapped across multiple source groups: {test}")

    runner_groups = make_list(RUNNER, "GROUPS")
    expected_groups = list(GROUP_RUNNERS)
    if runner_groups != expected_groups:
        errors.append(
            "tests/GNUmakefile GROUPS must be exactly: " + " ".join(expected_groups)
        )

    missing_group_runners = [
        group for group, path in GROUP_RUNNERS.items() if not path.is_file()
    ]
    for group in missing_group_runners:
        errors.append(f"test group has no GNUmakefile: {group}")

    runner_suites = make_list(GROUP_RUNNERS["src"], "SUITES")
    problem_suites = make_list(GROUP_RUNNERS["problems"], "SUITES")
    driver_suites = make_list(GROUP_RUNNERS["driver"], "SUITES")
    build_suites = make_list(GROUP_RUNNERS["build"], "SUITES")
    all_runner_suites = (
        runner_suites + problem_suites + driver_suites + build_suites
    )

    for suite in sorted(set(mapped_suites) - set(runner_suites)):
        errors.append(f"mapped suite is missing from tests/src SUITES: {suite}")
    for suite in sorted(set(runner_suites) - set(mapped_suites)):
        errors.append(f"tests/src SUITES contains an unmapped suite: {suite}")
    for suite in duplicates(mapped_suites):
        errors.append(f"suite contains more than one mapped source: {suite}")
    for suite in duplicates(runner_suites):
        errors.append(f"suite is listed more than once in tests/src SUITES: {suite}")

    for suite in sorted(set(problem_mapped_suites) - set(problem_suites)):
        errors.append(f"mapped suite is missing from tests/problems SUITES: {suite}")
    for suite in sorted(set(problem_suites) - set(problem_mapped_suites)):
        errors.append(f"tests/problems SUITES contains an unmapped suite: {suite}")
    for suite in duplicates(problem_mapped_suites):
        errors.append(f"problem suite contains more than one mapped source: {suite}")
    for suite in duplicates(problem_suites):
        errors.append(f"suite is listed more than once in tests/problems SUITES: {suite}")

    for suite in duplicates(all_runner_suites):
        errors.append(f"suite is registered in more than one runner group: {suite}")

    registered = set(all_runner_suites)
    suite_directories = {
        path.name
        for path in RUNNER.parent.iterdir()
        if path.is_dir()
        and path.name not in GROUP_RUNNERS
        and (path / "GNUmakefile").is_file()
    }
    for suite in sorted(registered - suite_directories):
        errors.append(f"registered suite directory does not exist: {suite}")
    for suite in sorted(suite_directories - registered):
        errors.append(f"suite directory is not registered in the runner: {suite}")

    for suite in sorted(registered & suite_directories):
        make_text = (RUNNER.parent / suite / "GNUmakefile").read_text(encoding="utf-8")
        for target in ("all", "run", "clean"):
            if f"{target}:" not in make_text:
                errors.append(f"suite has no {target} target: {suite}")

    group_implementation = RUNNER.parent / "make" / "group.mk"
    if not group_implementation.is_file():
        errors.append("shared test-group implementation is missing: tests/make/group.mk")
    else:
        implementation_text = group_implementation.read_text(encoding="utf-8")
        for target in ("all", "run", "run-suite", "list", "clean"):
            if f"{target}:" not in implementation_text:
                errors.append(f"test groups have no {target} target")

    for group, group_runner in GROUP_RUNNERS.items():
        if not group_runner.is_file():
            continue
        group_text = group_runner.read_text(encoding="utf-8")
        if "include ../make/group.mk" not in group_text:
            errors.append(f"test group does not use the shared runner: {group}")

    runner_text = RUNNER.read_text(encoding="utf-8")
    for target in (
        "all",
        "run",
        "run-src",
        "run-problems",
        "run-driver",
        "run-build-system",
        "run-integration",
        "run-suite",
        "list",
        "clean",
    ):
        if f"{target}:" not in runner_text:
            errors.append(f"tests/GNUmakefile has no {target} target")

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1

    print(
        f"OK ({len(sources)} active src files and {len(problem_sources)} "
        "non-driver problem files mapped one-to-one to primary tests; "
        f"{len(all_runner_suites)} suites in {len(GROUP_RUNNERS)} groups)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
