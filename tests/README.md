# Test suite

[Repository README](../README.md)

All commands in this guide are run from the repository root.

The root Makefile delegates to `tests/GNUmakefile`, which delegates to the
`src`, `problems`, `driver`, and `build` group GNUmakefiles; those in turn
delegate to the individual suites. The active configuration is forwarded
through every level. These effective defaults come from root `m_options` and
the hierarchical test makefiles and may be overridden on the command line:

```text
PFUNIT_ROOT=/opt/lib/pfunit/PFUNIT-4.16
MPIEXEC=/opt/lib/mpich-5.0.0/bin/mpiexec
MPIEXEC_FLAGS=
MPI_NP_FLAG=-n
MPIFC=/opt/lib/mpich-5.0.0/bin/mpif90
MUMPS_DIR=/opt/lib/MUMPS_5.8.2
SUITE_TIMEOUT=600s
DRIVER_CLI_TIMEOUT=20s
DRIVER_SMOKE_TIMEOUT=60s
DRIVER_INTEGRATION_TIMEOUT=90
SKIP_MPI_CASES=0
PERFORMANCE_BUILD_ROOT=build/openmp-performance
PERFORMANCE_TIMEOUT=300
PERFORMANCE_SUITE_TIMEOUT=3600s
PERFORMANCE_WARMUPS=1
PERFORMANCE_SAMPLES=3
PERFORMANCE_MIN_SPEEDUP=1.10
PERFORMANCE_MAX_REGRESSION=1.15
PERFORMANCE_BASELINE=
COVERAGE_ROOT=build/coverage
COVERAGE_BUILD_ROOT=build/coverage/build
COVERAGE_FLAGS=-O0 -g --coverage -fprofile-abs-path
COVERAGE_MIN_LINES=90.0
COVERAGE_MIN_FUNCTIONS=90.0
COVERAGE_MIN_BRANCHES=50.0
GCOV=gcov
LCOV=lcov
GENHTML=genhtml
OMP_PROC_BIND=close
OMP_PLACES=cores
```

## One source file, one primary test file

Every active library source in `src/sources.mk` and every non-driver problem
source in the per-problem `SOURCES` manifests has exactly one primary,
authored test file. The tab-separated mappings are stored in
`tests/test-map.tsv` and `tests/problem-test-map.tsv`. Each problem's
`main.F90` is exercised by the driver CLI, smoke, and integration layers. A
unit-test suite may still use fixtures, probes, stubs, generated pFUnit
sources, or link other production modules; those support files are not
additional primary tests for the mapped source.

Validate this invariant before changing or running the suites:

```bash
make test-layout
```

`check-layout` asks the `src` and `problems` owners for their active source
lists and rejects missing mappings, inactive sources, duplicate sources or
test files, and paths that do not exist. It keeps both the `tests/src` and
`tests/problems` suite manifests synchronized with their maps, validates the
four group runners, and verifies that every unit suite references its
production source and primary test. All registered library, problem, driver,
and framework/build-system suites must have `all`, `run`, and `clean` targets;
unregistered suite directories are rejected. The four group manifests
currently register 55 suites: 28 library, 23 problem, one driver, and three
framework/build-system suites. Problem modules are kept in separate suites because
several drivers deliberately use the same Fortran module names (`input_data`
and `RHS_fun`).

## Test targets

The runner separates library tests, problem-specific callback tests,
full-driver tests, and build-system tests:

```bash
# Check only the one-to-one layout.
make test-layout

# Run the 28 suites mapped to src/*.F90.
make test-src

# Run all 23 problem-specific input, RHS, and solver suites.
make test-problems

# Build all ten problem executables and run CLI, smoke, numerical integration,
# and the isolated release performance gate.
make test-driver

# Exercise the hierarchical Make interface in an isolated build tree.
make test-build-system

# Individual driver layers.
make test-cli
make test-smoke
make test-integration

# Build an isolated release executable and gate OMP1/OMP4 scaling.
make test-performance

# Exercise only the performance-gate logic, without compiling or launching MPI.
make test-performance-self-test

# Instrument and run the functional suite, build HTML/JSON/LCOV reports, and
# enforce independent source-line, function, and branch thresholds.
make test-coverage

# Run the complete regression above in one command.
make test
make check

# Clean-build every test without executing it, list suites, or clean-run one.
make test-build
make test-list
make test-suite TEST_SUITE=rhs_assembly
```

`test-coverage` performs a clean, isolated GNU build below `build/coverage`,
runs the build-system, library, problem, CLI, smoke, and numerical integration
tests, and measures only the active core sources declared by
`src/sources.mk`. It excludes the timing-based performance suite because gcov
instrumentation changes runtime. Every manifest source must be represented in
the LCOV denominator; the two declaration-only modules currently omitted by
gcov (`Interfaces.F90` and `projection_engine.F90`) are accepted only after a
fresh `.gcno` file and `gcov` itself confirm that they contain no executable
lines.

The default gates are 90% executable source lines, 90% functions, and 50%
branches. Override `COVERAGE_MIN_LINES`, `COVERAGE_MIN_FUNCTIONS`, or
`COVERAGE_MIN_BRANCHES` on the command line when deliberately changing the
policy. The run writes:

```text
build/coverage/coverage.info
build/coverage/coverage-summary.json
build/coverage/html/index.html
```

Use `make clean-coverage` to remove only owned coverage builds and reports.
The target refuses unsafe, overlapping, symlinked, nonempty unowned, or
incorrectly marked roots; unrelated files later placed in an owned coverage
root are retained. Per-suite `.gcda` and `.gcno` files are removed after every
attempted run, including a test or threshold failure; metadata inside the
isolated owned build tree stays there until `clean-coverage`.

Every test level can be invoked directly. Driver executables built in
`mymake/EXEC` are deliberately retained by `clean-tests`; use `clean-build` or
`clean` to remove them.

```bash
make -j1 -C tests run-src
make -j1 -C tests/problems run-suite TEST_SUITE=heat_rhs_fun
make -j1 -C tests/build run-suite TEST_SUITE=make_hierarchy
make -j1 -C tests/rhs_assembly run
```

The aggregate `run-src`, `run-problems`, `run-driver`, and `run-build-system`
targets clean each suite before running it, so compiler or flag changes cannot
silently reuse a stale test executable. Driver executables are rebuilt
unconditionally.

Pass non-default tool locations once at the root; they are forwarded to every
suite:

```bash
make test \
  PFUNIT_ROOT=/path/to/pfunit \
  MPIEXEC=/path/to/mpiexec \
  MPIFC=/path/to/mpif90 \
  MUMPS_DIR=/path/to/mumps
```

The runner is deliberately serialized to keep diagnostics deterministic and
avoid oversubscribing MPI/OpenMP test jobs. Each problem build nevertheless
has a private `_OBJ` directory, so identically named problem modules cannot be
reused accidentally. Each suite is protected by `SUITE_TIMEOUT`; the driver
suite that includes the timed gate uses `PERFORMANCE_SUITE_TIMEOUT`. MPI suites
exercise up to eight ranks. The end-to-end problem matrix holds the MPI
topology fixed while comparing one and four OpenMP threads. Lower-level
parallel tests additionally exercise thread counts 2 and 8. The top-level
runner requires a POSIX environment with Bash and the coreutils `timeout`
command; selected error-path probes additionally use POSIX process primitives.

## Detailed numerical gates

- [Positive smoke and numerical integration matrix](integration.md)
- [OpenMP performance regression](performance.md)
- [Benchmark framework tests](benchmarking/README.md)
