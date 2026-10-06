# Architecture and extension contract

[Benchmarking index](../README.md) · [Repository README](../../README.md)

```text
configs/*.json
       |
       v
catalog.py --> neutral planner/executor/storage --> benchmarks/<run-id>/
       |                         |
       |                         +--> MPI launcher
       |                         +--> registered analysis pipeline
       v
ManufacturedTransientAdapter (Python command/parser/validation)
       |
       v
three thin Fortran entry points and adapters
       |
       v
shared contract + manufactured case + lifecycle + measurement
       |
       v
public ADS library and problem-specific production RHS path
```

The modules in `ads_benchmark/framework/` never import concrete problems.
They receive extension points through `Catalog` registries.
`ads_benchmark/catalog.py` is the composition root. Adding a component means
implementing its protocol and registering it there; the generic planner and
executor do not branch on problem names.

The registered extension axes are:

- problem adapter;
- experiment family (`temporal`, `h`, `p`, `validation`, `strong`,
  or `weak`);
- exact/manufactured case;
- build profile;
- launcher;
- analyzer and the experiment-family-to-analyzer binding.

The generic executor owns the case working directory, environment, OpenMP
settings, process group, timeout, logs, status transitions, and common result
envelope. A Python problem adapter validates one planned case, constructs the
payload argv, parses the tagged solver record, and binds it back to the exact
planned configuration.

Analysis is another registry-backed extension point. The generic pipeline
loads results against their frozen manifest and dispatches them to a family
analyzer. The temporal and spatial analyzers own only their scientific
grouping, qualification, plateau handling, and reports; they do not reimplement
process execution, resume, or result storage.

There is one implementation owner for each cross-cutting responsibility:

| Responsibility | Single implementation owner |
| --- | --- |
| Process lifecycle, process groups, timeout, retry, and failure classification | `framework/executor.py` (`Executor`) |
| Capturing and publishing per-case stdout/stderr | `framework/executor.py`, through the store API |
| Resume orchestration and exact identity checks | `cli.py`, delegating completed-case revalidation to `Executor` |
| Contained atomic run/case/analysis storage | `framework/storage.py` (`ResultStore`) |
| MPI/direct/scheduler argv and resource validation | registered launchers in `components/planning.py` |
| Complete-field comparison | `validation/fields.py` |
| Median, MAD, range, speedup, and efficiency primitives | `analysis/statistics.py` and `analysis/scaling.py` |

Problem and experiment code depends on these owners; the owners do not import
concrete problems. In particular, no analyzer starts processes, no adapter
writes logs, and no comparison command implements a second result loader.

The Fortran side follows the same separation. A registered
`BenchmarkAdapter` supplies a manufactured-case descriptor and procedures
for initialize, initial projection, one physical step, measurement, and
cleanup. `benchmark_harness.F90` owns their ordering and the `1..N` time
loop. The three main programs only register one adapter and invoke that shared
harness.

## How to add a benchmark

A new benchmark over an already supported problem normally needs data and a
small analysis extension, not a new runner:

1. Add a strict profile in `configs/<profile>.json`, using the existing case
   schema for family, exact case, time package, mesh, degrees, MPI/OpenMP,
   sampling, measurement, build profile, and launcher.
2. If it is a new experiment family, add one `FamilyDefinition` and its
   analyzer to `catalog.py`. Reuse `CaseSpec`, `Planner`, `Executor`,
   `ResultStore`, field validation, and statistics; do not add a parallel
   process/log/resume path.
3. Add focused tests under `tests/benchmarking/` for exact expansion,
   validation, stable case IDs, synthetic analysis, and a no-write `make benchmark-plan
   BENCHMARK_PROFILE=<profile>` check.
4. Expose a root Make target only when the workflow needs more than the generic
   plan, run/resume, analyze, or compare operations. Generated data must remain
   under `benchmarks/<run-id>/`.

## How to add a problem adapter

1. Implement the structural `ProblemAdapter` protocol from
   `framework/protocols.py`: `name`, `execution_ready`, `validate_case`,
   `build_payload_command`, `parse_result`, and `validate_result`. Return an
   argv vector, never shell syntax; leave process management, MPI wrapping,
   OpenMP environment, timeout, logs, retry, and storage to `Executor`.
2. Register the adapter once in `catalog.py`. Add an exact-case definition or
   analyzer registration there only if the new problem actually needs one.
3. Add the benchmark-local Fortran entry point/lifecycle adapter and its
   isolated build rule in `benchmarking/GNUmakefile`; do not change a public
   problem CLI merely to serve the benchmark.
4. Test command construction, strict tagged-result parsing, planned/result
   binding, domain validation, and registration with a fake or minimal
   adapter. The generic planner and executor must remain unchanged.
