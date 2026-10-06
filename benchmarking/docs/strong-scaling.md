# Strong scaling

[Benchmarking index](../README.md) · [Repository README](../../README.md)

The full matrix is deliberately manual. A structural dry-run performs no
build or result writes:

```bash
make benchmark-plan BENCHMARK_PROFILE=cluster-scaling
```

Real runs require explicit allocation limits. The public target builds the
isolated `release` tree selected by `CONFIG`, executes the frozen plan, runs
field-gated analysis, and requests the three-panel scaling plot:

```bash
make benchmark-strong RUN_ID=strong-full \
  BENCHMARK_PLAN_ARGS='--available-mpi-slots 8 --available-cpu-slots 64'
```

`MPIEXEC`, `MPI_NP_FLAG`, and shell-parsed `MPIEXEC_FLAGS` select the local or
cluster launcher; no scheduler or absolute launcher path is embedded in a
profile. `openmp.dynamic` is fixed to false for strong scaling, while
`openmp.proc_bind` and `openmp.places` are profile data and therefore frozen
in every case identity and result.

A practical local verification slice contains two resource levels, DG/PR,
and two degree pairs:

```bash
make benchmark-strong RUN_ID=strong-local-check \
  BENCHMARK_STRONG_PROFILE=strong-scaling-smoke \
  BENCHMARK_PLAN_ARGS='--available-mpi-slots 2 --available-cpu-slots 2'
```

Resume uses the same profile, filters, allocation, launcher, source
fingerprint, and run ID through `make benchmark-resume`. A changed launcher
prefix or normalized configuration is rejected instead of being mixed into
one run. Any filtered scaling slice must retain `MPI=1,OMP=1` for every
selected problem/scheme/degree/mesh combination, because that case supplies
the mandatory full-field reference.

```bash
make benchmark-resume RUN_ID=strong-local-check \
  BENCHMARK_RESUME_PROFILE=strong-scaling-smoke \
  BENCHMARK_PLAN_ARGS='--available-mpi-slots 2 --available-cpu-slots 2'
```

Each strong case starts the identical executable nine times: two validated
warmups followed by seven measured repetitions. The primary sample is the
Fortran `physical_step_wall_seconds`, never process wall time. The interval
starts after initialization, initial projection, and the pre-step MPI barrier;
it ends immediately after the fixed physical-step package. Exact-field
evaluation, regular-grid sampling, CSV output, analysis, and cleanup stay
outside it. `MPI_Wtime` is reduced with `MPI_MAX`, so every sample represents
the slowest rank. The safety status reduction performed by each physical step
is part of this package.

All raw warmup, measured, and process-wall samples are retained. Analysis
reports median, minimum, maximum, median absolute deviation (MAD), and relative
MAD. A sample below the profile's `minimum_sample_seconds` marks the
configuration unreliable; it remains visible in JSON/CSV but is excluded from
speedup and efficiency. No configuration silently receives a different number
of physical steps.

Before a timing is eligible, its complete numerical field must match the
otherwise identical `MPI=1, OMP=1` field with absolute tolerance `1e-11` and
relative tolerance `1e-10`. Speedup uses the median of the smallest eligible
resource configuration for the identical global problem. With
`R = MPI ranks x OpenMP threads`, efficiency is
`speedup / (R/R_reference)`. Equal-resource candidates are resolved
deterministically by MPI rank count, process-grid vector, OpenMP thread count,
and case ID; the selected case is recorded in every series. The analyzer emits
`analysis.json`, a flat CSV
with every sample and MPI/OpenMP coordinate, and `strong-scaling.png` with
time, speedup, and efficiency panels. Plot labels contain the full problem,
scheme, degree vectors, global mesh, MPI grid, OpenMP values, and stable series
ID. To prevent an unreadable plot from silently hiding thousands of identities,
plotting is explicitly omitted above 12 series; use the normal plan filters to
select an interpretable slice. JSON and CSV are always complete.

The full profile represents 12,960 configurations and 116,640 solver
invocations before retries. Final `17^3` field CSV files can consume roughly
10 GB in aggregate. Planning it is not a claim that the matrix was executed.
