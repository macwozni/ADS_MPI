# ADS MPI IGA

ADS MPI IGA is a Fortran/MPI implementation of the alternating direction
solver (ADS) for isogeometric analysis, including the current iGRM workflow,
ten problem drivers, an automated test hierarchy, and reproducible benchmark
workflows.

The repository-root `Makefile` is the public entry point. Source ownership
remains in `src/`, problem ownership in `problems/`, tests in `tests/`,
and benchmark implementation in `benchmarking/`.

> Current development scope is GCC/GFortran, MPI, and MUMPS. Intel compiler
> configurations and ParMETIS are outside the supported scope for ongoing
> work; legacy configuration entries must not be treated as qualified paths.

## Quick start

Run these commands from the repository root:

```bash
make show-config
make                         # core library and all problem drivers
make run-help                # arguments for every problem
make test-layout             # fast test/source ownership audit
make benchmark-self-test     # benchmark contracts + short real-MPI integration
make docs-check
```

Build profiles are selected through the active `m_options` file:

```bash
make BUILD=debug
make BUILD=release
make CONFIG=makeconfig/gnu-debug.mk all
make CONFIG=makeconfig/gnu-release.mk all
```

Generated executables are written below `mymake/EXEC/`. A small problem can
be launched through the root interface, for example:

```bash
make run-heat
make run-heat ARGS='4 2 3 0.01 2 1 1' NP=2 OMP_NUM_THREADS=4
```

## Documentation

The landing pages stay intentionally short. Detailed material is organized by
topic:

| Topic | Document |
|---|---|
| Documentation index | [`docs/README.md`](docs/README.md) |
| Build, dependencies, cleanup, and Doxygen | [`docs/building.md`](docs/building.md) |
| Core API, time schemes, and iGRM mesh contract | [`docs/architecture.md`](docs/architecture.md) |
| Running and runtime controls | [`docs/running.md`](docs/running.md) |
| Classic ADS problem drivers | [`docs/problems-classic.md`](docs/problems-classic.md) |
| iGRM and DPG problem drivers | [`docs/problems-igrm.md`](docs/problems-igrm.md) |
| Test hierarchy and commands | [`tests/README.md`](tests/README.md) |
| Benchmark framework | [`benchmarking/README.md`](benchmarking/README.md) |
| Benchmark framework tests | [`tests/benchmarking/README.md`](tests/benchmarking/README.md) |

## Repository layout

```text
Makefile              Public build, run, test, benchmark, and docs interface
m_options             Active local GCC/GFortran/MPI/MUMPS configuration
makeconfig/           GNU debug and release configuration examples
src/                  Core ADS library and ordered source manifest
problems/             Ten problem drivers with local build/run ownership
tests/                Unit, integration, performance, and framework tests
benchmarking/         Benchmark implementation, profiles, adapters, and docs
benchmarks/           Ignored generated benchmark runs and preserved user data
docs/                 Tracked topic-oriented project documentation
mymake/               Compatibility entry point and generated build artifacts
doxygen/              Ignored generated HTML/PDF documentation
```

The CMake files are retained for historical compatibility, but the maintained
workflow documented here is Make-based.

## Build and run

The normal build is serialized where necessary because several problem-local
Fortran modules intentionally reuse names such as `input_data`, `RHS_fun`,
and `main`.

```bash
make
make library
make problems
make build PROBLEM=heat
make build-igrm_l2
make list-problems
make show-run PROBLEM=heat
```

MPI rank count must match the process grid encoded in a problem's arguments:

```text
NP = procx * procy * procz
```

See the [build guide](docs/building.md) for configuration and cleanup, and the
[running guide](docs/running.md) for shared MPI/OpenMP controls.

## Testing

The root test runner delegates through four groups: core sources, problem
sources, drivers, and build/framework tests.

```bash
make test-layout
make test-src
make test-problems
make test-driver
make test-build-system
make test
```

The complete `make test` workflow is substantial and includes real MPI,
numerical integration, and a release OpenMP performance gate. The benchmark
test suite adds only a bounded seven-case correctness integration; it does not
execute any full benchmark profile. See [`tests/README.md`](tests/README.md)
for exact scope and smaller targets.

## Benchmarking

The benchmark framework covers manufactured temporal, h/p convergence,
MPI/OpenMP field validation, strong and weak scaling, sharding, provenance,
resume, and numerics-first A/B comparison.

```bash
make benchmark-plan
make benchmark-build BUILD=release
make benchmark-smoke BENCHMARK_RUN_ID=smoke-local
make benchmark-self-test
```

Full profiles remain explicit workflows. Merely planning a profile or passing
the small integration test is not a scientific qualification. The latest
temporal status is documented in the
[convergence guide](benchmarking/docs/convergence.md) and its
[reproducers](benchmarking/reproducers/).

Every generated run belongs below ignored `benchmarks/<run-id>/`. Ordinary
`make clean` and `make clean-benchmark-build` never remove run data.

## Documentation build

```bash
make docs-check
make docs-html
make docs-pdf
make clean-docs
```

Doxygen inputs are explicitly limited to tracked source and documentation
areas; local prompts, benchmark results, and unrelated scratch files are not
documentation inputs.
