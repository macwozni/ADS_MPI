# Planning, execution, provenance, and resume

[Benchmarking index](../README.md) · [Repository README](../../README.md)

## Plan, build, and run

From the repository root:

```bash
# No-write planning.
make benchmark-plan
make benchmark-plan BENCHMARK_PROFILE=temporal-full

# Isolated ADS library plus all three adapters.
make benchmark-build BUILD=debug

# Build and run N=4 and N=8, then require improvement for all nine pairs.
make benchmark-smoke BENCHMARK_RUN_ID=smoke-local

# Plan the complete and diagnostic temporal matrices without running them.
make benchmark-plan BENCHMARK_PROFILE=temporal-full
make benchmark-plan BENCHMARK_PROFILE=temporal-validation

# Plan the spatial matrices without running them.
make benchmark-plan BENCHMARK_PROFILE=h-convergence-full
make benchmark-plan BENCHMARK_PROFILE=p-convergence-full
make benchmark-plan BENCHMARK_PROFILE=p-anisotropic-full

# Execute frozen full workflows; select the matching *-smoke profile for smoke.
make benchmark-h-convergence RUN_ID=h-local
make benchmark-p-convergence RUN_ID=p-local
make benchmark-p-convergence RUN_ID=p-anisotropic-local \
  BENCHMARK_P_CONVERGENCE_PROFILE=p-anisotropic-full

make benchmark-h-convergence RUN_ID=h-smoke-local \
  BENCHMARK_H_CONVERGENCE_PROFILE=h-convergence-smoke
make benchmark-p-convergence RUN_ID=p-smoke-local \
  BENCHMARK_P_CONVERGENCE_PROFILE=p-convergence-smoke

# Weak-scaling matrices: planning writes nothing and needs no allocation.
make benchmark-plan BENCHMARK_PROFILE=local-weak-scaling
make benchmark-plan BENCHMARK_PROFILE=cluster-weak-scaling

# Real two-level weak-scaling check on a declared local allocation.
make benchmark-weak RUN_ID=weak-smoke-local \
  BENCHMARK_WEAK_PROFILE=weak-scaling-smoke \
  BENCHMARK_PLAN_ARGS='--available-mpi-slots 2 --available-cpu-slots 2'
```

`benchmark-smoke` requires a new base ID. It creates
`benchmarks/smoke-local-n4/` and `benchmarks/smoke-local-n8/`; existing
run directories are never overwritten. The lower-level equivalents include:

```bash
make -C benchmarking build BUILD=release
make -C benchmarking build-igrm_heat BUILD=debug
make -C benchmarking smoke BENCHMARK_RUN_ID=smoke-local
```

To run a profile or filtered subset explicitly:

```bash
cd benchmarking
MPIEXEC=/path/to/mpiexec python3 -m ads_benchmark run \
  --profile smoke --run-id smoke-one \
  --problem igrm_heat --scheme dg
```

The standalone Fortran form, useful for diagnosis, is:

```bash
OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close \
  mpiexec -n 1 \
  ./build/debug/EXEC/igrm_l2_manufactured \
  dg 0.1 4 3 3 3 4 4 4 3 3 3 1 1 1 17 0 temporal-polynomial
```

Repeatable planner/runner filters include `--problem`, `--scheme`,
`--degree-pair`, `--mesh`, `--mpi-grid`, `--mpi-ranks`, `--omp`, and
`--steps`. Optional `--available-mpi-slots` and `--available-cpu-slots`
declare an actual allocation. The planner then rejects a rank count above the
first limit or `MPI ranks * OpenMP threads` above the second. Without those
arguments the plan remains a structural, portable dry-run; it does not guess
scheduler capacity from the login host. A real `validation` run or resume
requires both declarations, so an unknown or insufficient allocation cannot
silently enter the correctness matrix.
Unknown values and an empty result fail explicitly. Validation includes
`p_test > p_trial`, trial degree at least three for `temporal-polynomial` and
at least one for `spatial-cosine`,
maximum degree nine, positive dimensions, `NP=proc-x*proc-y*proc-z`, a
distributable process grid, and positive runtime controls.

The isolated build lives only under `benchmarking/build/<profile>/` and uses
the selected repository `CONFIG`, compiler, and libraries. No MPI or MUMPS
path is hardcoded. `make clean-benchmark-build` removes only marker-owned
benchmark build/cache content, never benchmark results.

## Frozen runner, provenance, resume, and failures

A new temporal run builds the release adapters, exclusively creates its run
directory, writes the complete manifest, reads that manifest back through the
strict decoder, and executes only the decoded frozen cases:

```bash
make benchmark-convergence \
  RUN_ID=temporal-validation-local \
  BENCHMARK_CONVERGENCE_PROFILE=temporal-validation
```

An existing run ID is never adopted or overwritten by a new run. To continue
an interrupted run, request the same profile and filters explicitly:

```bash
make benchmark-resume \
  RUN_ID=temporal-validation-local \
  BENCHMARK_RESUME_PROFILE=temporal-validation
```

Transient retries are opt-in. The same public target can, for example, make up
to two additional attempts after each retryable failure:

```bash
make benchmark-resume \
  RUN_ID=temporal-validation-local \
  BENCHMARK_RESUME_PROFILE=temporal-validation \
  BENCHMARK_PLAN_ARGS='--max-retries 2'
```

`--max-retries` is also accepted by new-run and shard-run workflows through
`BENCHMARK_PLAN_ARGS`; its default is zero.

The same neutral resume variable applies to spatial runs, for example
`BENCHMARK_RESUME_PROFILE=h-convergence-full` or
`BENCHMARK_RESUME_PROFILE=p-convergence-full`. For compatibility,
`BENCHMARK_RESUME_PROFILE` defaults to `BENCHMARK_CONVERGENCE_PROFILE` when it
is not set explicitly.

Resume requires the current profile expansion, filters, commit SHA, dirty
state, and content fingerprint of the nonignored worktree to match the frozen
manifest. The expanded `config_hash` must match too. Thus two different dirty
source trees are not treated as compatible.
It holds an exclusive execution lock and skips a case only when `status.json`,
`result.json`, both logs, the complete normalized configuration, and the
adapter's domain result all revalidate. The tagged stdout is parsed again by
the registered adapter and must reproduce the saved domain result exactly.
Missing, failed, timed-out, incomplete, or tampered cases are retried; a
different configuration or source state is refused.

Every newly executed run also has a schema-versioned `execution.json`. Its
`manifest_hash` binds it transitively to the manifest's full Git SHA, dirty
flag, nonignored-worktree fingerprint, expanded configuration, and case set.
The execution record adds:

- compiler command and version, compile/link flags, debug/release profile, and
  hashes of the executable and build stamps;
- configured and linked MUMPS, BLAS, LAPACK, ScaLAPACK, METIS, GKlib, and
  other libraries, including the MUMPS version where discoverable;
- MPI or scheduler launcher kind, executable, resolved path, argv, version,
  and declared MPI/CPU/thread capacity;
- hostname, OS, kernel, architecture, CPU model, physical/logical core counts,
  physical memory, timestamp, timezone name, and UTC offset;
- every planned rank grid and OpenMP binding tuple, plus relevant
  `OMP_*`, `GOMP_*`, and `KMP_*` environment values.

Each unavailable observation is represented explicitly as
`{"value": null, "reason": "..."}`; the collector does not invent a version
or hardware fact. Present observations use `{"value": ..., "reason": null}`.
The record carries separate SHA-256 hashes for the build identity, execution
compatibility identity, and whole record.

A resume now needs both an exactly compatible frozen manifest and a valid
`execution.json` whose build/launcher compatibility hash equals the current
one. It still checks the exact launcher argv of every reusable passed case.
Historical runs without execution provenance remain readable by offline
analysis, but they cannot be resumed into a mixed build.

Before any retry, every reusable completed case must also have the same full
launcher command as the current resume request. A changed `MPIEXEC`, rank
flag, or other launcher prefix therefore refuses the resume before creating or
rewriting any case; offline analysis remains independent of the currently
installed launcher.

The lower-level equivalents are `make -C benchmarking convergence ...` and
`make -C benchmarking resume ...`. Cases move explicitly through
`planned`, `running`, `passed`, `failed`, or `timeout`. A successful result is
written atomically only after process exit, tagged-record parsing, and domain
validation, so interruption cannot create an apparently completed case.
Manifest, execution, status, result, log, and analysis publication is contained
by `ResultStore` and uses atomic/exclusive writes as appropriate. A timeout
terminates the whole process group and escalates to `SIGKILL` if necessary.

Every unsuccessful attempt has exactly one `failure_kind`: `numerical`,
`mpi`, `timeout`, `resource`, or `configuration`. Only `mpi`, `timeout`, and
`resource` are retried; numerical and configuration failures stop immediately.
Because MPI launchers propagate payload exit codes, an ordinary nonzero exit is
conservatively `numerical` even when the command is wrapped by `mpiexec` or a
scheduler. It is never retried. The `mpi` category is reserved for an explicit,
typed launcher/control-plane signal; it is not inferred from the launcher name,
exit code, or diagnostic text. A custom launcher supplies that signal through
the optional `classify_process_failure` hook, which may return only `"mpi"` or
`None`. Invalid values and hook errors are configuration failures rather than
retryable MPI failures. The built-in `mpiexec` and command-template launchers
do not infer such a signal, so an untyped launcher failure remains
conservatively `numerical`; installations that can expose a reliable
control-plane side channel must register a launcher implementing the hook.
When retries are enabled, final `status.json` retains a compact `attempts`
history. `stdout.log` and `stderr.log` deliberately contain only the last
attempt; for a repeated timing case that means the last launched process of
that attempt. There is no per-attempt or per-sample log archive, so use the
status history for earlier classifications.
Before a new run directory is created, execution preflight validates every
registered adapter/launcher command and the availability of both payload and
launcher executables. Invalid decompositions are rejected during ordinary
case validation. These errors therefore cannot leave a nominal run behind.

Process completion is not equivalent to a passing convergence analysis; the
selected analyzer must still accept the scientific series.
