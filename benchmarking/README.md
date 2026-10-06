# ADS MPI benchmarking framework

This directory owns the reusable benchmark implementation: strict profiles,
the neutral planner and executor, three manufactured-problem adapters,
field/numerical validation, scaling analysis, provenance, resume, sharding,
and A/B comparison.

Generated manifests, logs, fields, and reports belong only below the ignored
`benchmarks/<run-id>/` tree. The framework never adopts or deletes unrelated
user data.

> Current development scope is GCC/GFortran, MPI, and MUMPS. Intel compiler
> configurations and ParMETIS are outside the supported benchmark scope.

## Quick start

Run public commands from the repository root:

```bash
make benchmark-plan
make benchmark-build BUILD=release
make benchmark-smoke BENCHMARK_RUN_ID=smoke-local
make benchmark-self-test
```

Planning performs no execution and does not create a run unless a run workflow
is explicitly selected. `benchmark-smoke` builds and runs the bounded
three-problem × three-scheme matrix at `N=4` and `N=8`, then requires every
pair to improve.

Typical explicit workflows are:

```bash
make benchmark-convergence RUN_ID=temporal-local \
  BENCHMARK_CONVERGENCE_PROFILE=temporal-validation

make benchmark-validate RUN_ID=field-local \
  BENCHMARK_PLAN_ARGS='--available-mpi-slots 2 --available-cpu-slots 8'

make benchmark-strong RUN_ID=strong-local \
  BENCHMARK_STRONG_PROFILE=strong-scaling-smoke \
  BENCHMARK_PLAN_ARGS='--available-mpi-slots 2 --available-cpu-slots 2'

make benchmark-weak RUN_ID=weak-local \
  BENCHMARK_WEAK_PROFILE=weak-scaling-smoke \
  BENCHMARK_PLAN_ARGS='--available-mpi-slots 2 --available-cpu-slots 2'
```

A run ID is immutable. Continue an interrupted compatible run with
`benchmark-resume`; do not start a new run with an existing identifier.

## Design guarantees

- strict JSON profiles expand deterministically into immutable cases;
- process execution, timeout, logging, retry, and status transitions have one
  owner in `Executor`;
- atomic, path-contained persistence has one owner in `ResultStore`;
- adapters own problem-specific command construction, parsing, and validation;
- complete sampled fields gate timing eligibility;
- configuration, build, launcher, machine, MPI, and OpenMP provenance is
  recorded explicitly;
- numerical agreement is checked before A/B timing statistics;
- ordinary nonzero payload exits are `numerical`, including under an MPI
  launcher, and are never retried;
- results are never versioned and ordinary cleanup never removes them.

Full convergence and scaling matrices are deliberately separate from
`make test`. Passing the framework tests or a smoke run does not imply that
an unexecuted full profile passed.

## Documentation

| Topic | Document |
|---|---|
| Architecture and extension API | [`docs/architecture.md`](docs/architecture.md) |
| Manufactured cases and adapters | [`docs/cases-and-adapters.md`](docs/cases-and-adapters.md) |
| Profiles, matrix sizes, and case IDs | [`docs/profiles.md`](docs/profiles.md) |
| Planning, execution, provenance, and resume | [`docs/running-and-resume.md`](docs/running-and-resume.md) |
| Temporal, h, and p convergence | [`docs/convergence.md`](docs/convergence.md) |
| Full-field MPI/OpenMP validation | [`docs/validation.md`](docs/validation.md) |
| Strong scaling | [`docs/strong-scaling.md`](docs/strong-scaling.md) |
| Weak and hybrid scaling | [`docs/weak-scaling.md`](docs/weak-scaling.md) |
| Sharded runs and merge | [`docs/sharding.md`](docs/sharding.md) |
| A/B numerical and performance comparison | [`docs/comparison.md`](docs/comparison.md) |
| Measurements, schemas, gates, and safety | [`docs/results-and-safety.md`](docs/results-and-safety.md) |
| Framework and real-MPI integration tests | [`../tests/benchmarking/README.md`](../tests/benchmarking/README.md) |

Scientific failure records and focused reproduction commands live in
[`reproducers/`](reproducers/).

## Current qualification status

`temporal-full` expands to 792 cases. Core fix `5e1161b` removed several
previously recorded transient failures, but post-fix spot checks still contain
a deterministic degree-dependent non-monotone sequence. A fresh complete run
of both temporal profiles and resolution of remaining outliers are required
before full scientific qualification. The oracle has not been relaxed.

The complete strong and weak cluster matrices are configuration artifacts
until they are run on a declared allocation and pass both field and statistical
gates. No documentation page should describe an unexecuted matrix as passed.

## Result location and cleanup

```text
benchmarks/<run-id>/
  manifest.json
  execution.json
  cases/<case-id>/
    status.json
    result.json
    stdout.log
    stderr.log
    field_samples.csv    # only when requested
  analysis/
```

`make clean-benchmark-build` removes only marker-owned build/cache data.
Neither it nor `make clean` removes `benchmarks/<run-id>/`. Inspect and
remove one exact run manually only when its data is no longer needed.
