# Weak and hybrid scaling

[Benchmarking index](../README.md) · [Repository README](../../README.md)

Strong and weak scaling answer different questions and are never pooled in
one series. A strong series keeps the global mesh fixed while resources grow.
A weak series in this repository fixes a three-dimensional workload per MPI
rank and derives the global mesh componentwise:

```text
global_elements[d] = local_elements[d] * process_grid[d]
```

The workload basis is explicitly `per-rank`. It is not per core: for fixed
local elements, raising the OpenMP thread count reduces work per core. Each
reported weak series therefore fixes problem, scheme, time package, degree
vectors, local element vector, OpenMP thread count, and the complete OpenMP
binding policy, then varies only the MPI process grid. The matrix still
contains all three useful views:

- MPI-only points use OMP1 across the MPI layouts;
- OpenMP-only points use MPI1 with OMP `1,2,4,8`;
- hybrid points use more than one rank and OMP `2,4,8`.

The OMP1, OMP2, OMP4, and OMP8 data remain distinct series. In particular,
the analyzer does not combine the MPI1/OMP1 and MPI1/OMP8 points into a
constant-work weak-scaling ladder, because their per-core workloads differ.

For a passing series the mandatory MPI1 member is the timing baseline. If its
median measured physical-step time is `t_1` and another level takes `t_r`, the
reported metric is

```text
weak_scaling_efficiency = t_1 / t_r
```

It is named only weak-scaling efficiency. It is not a strong-scaling speedup,
and no speedup field is emitted by the weak analyzer. A failed or unreliable
timing remains in JSON/CSV but cannot make the series pass.

Weak cases use the same repeated timing method: two validated warmups and seven
measured executable launches, a barrier before the fixed physical-step
package, `MPI_Wtime`, and `MPI_MAX` across ranks. Setup, initial projection,
regular-grid sampling, field output, analysis, and cleanup remain outside the
timed region. The report retains raw solver and process-wall samples plus
median, min, max, MAD, and relative MAD. The output is `analysis.json`, a flat
`analysis.csv`, and, when requested and readable, a two-panel
`weak-scaling.png` containing time and weak-scaling efficiency.

Every timed field is checked against exactly one MPI1/OMP1 field on the same
global mesh and with identical physics, degrees, time package, and sample
grid. A timing configuration can itself provide that reference. Otherwise the
planner adds a deduplicated `field-reference` helper case. Helpers are visible
in the frozen manifest and must pass, but are excluded from timing series and
efficiency calculations. This is necessary in weak scaling because different
MPI grids intentionally produce different global meshes; comparing all levels
to the small MPI1 baseline field would compare different discrete problems.

The profiles are:

- `weak-scaling-smoke`: local `4^3`, MPI1 and MPI2/`2x1x1`, OMP1; two
  measurements plus one helper;
- `local-weak-scaling`: local `8^3,12^3,16^3`, four process grids through
  `2x2x2`, and OMP `1,2,4,8`; 6,480 measurements plus 1,080 helpers;
- `cluster-weak-scaling`: the same local work and OpenMP levels over eight
  balanced and asymmetric grids through `4x4x4`; 12,960 measurements plus
  2,025 helpers.

The two complete profiles cover `igrm_l2`, `igrm_heat`, and
`pure_diffusion_igrm`, DG/PR/BE, and all 15 supported isotropic degree pairs.
They are release-only manual workflows. A real local smoke is:

```bash
make benchmark-weak RUN_ID=weak-local-check \
  BENCHMARK_WEAK_PROFILE=weak-scaling-smoke \
  BENCHMARK_PLAN_ARGS='--available-mpi-slots 2 --available-cpu-slots 2'
```

Resume uses the identical frozen selection, repository state, launcher, and
allocation declaration:

```bash
make benchmark-resume RUN_ID=weak-local-check \
  BENCHMARK_RESUME_PROFILE=weak-scaling-smoke \
  BENCHMARK_PLAN_ARGS='--available-mpi-slots 2 --available-cpu-slots 2'
```

Planning the full cluster profile validates all 14,985 cases without claiming
that any solver was launched:

```bash
make benchmark-plan BENCHMARK_PROFILE=cluster-weak-scaling
```

The local launcher remains the configured `MPIEXEC`, `MPIEXEC_FLAGS`, and
`MPI_NP_FLAG`. For a scheduler, `--launcher-template` (or
`BENCHMARK_LAUNCHER_TEMPLATE` on the Make workflow) supplies one shell-free
argv template. Supported placeholders are `{ranks}`, `{threads}`, `{procx}`,
`{procy}`, `{procz}`, and one mandatory final standalone `{payload}`. Hybrid
execution requires `{threads}`. For example:

```bash
cd benchmarking
python3 -m ads_benchmark.cli plan \
  --repository-root .. --config-dir configs \
  --profile weak-scaling-smoke \
  --available-mpi-slots 2 --available-cpu-slots 8 \
  --launcher-template \
    'srun --ntasks={ranks} --cpus-per-task={threads} {payload}' \
  --show-commands
```

`--show-commands` prints each final shell-quoted argv, including the complete
solver payload, without executing it. The template is tokenized once and is
never evaluated through a shell. Profiles therefore contain no account,
partition, node count, launcher path, or site-specific MPI installation. The
runner validates the requested rank count, `ranks * threads`, process-grid
product, allocation limits, and hybrid binding contract before execution.
SLURM's `SLURM_NTASKS` and `SLURM_CPUS_PER_TASK` can further constrain a
template launcher, but explicit `--available-mpi-slots` and
`--available-cpu-slots` keep an execution request auditable and portable.

A successful two-rank smoke proves that the real local executable, repeated
timing, full-field gate, and weak analyzer worked at those two resource
levels. It does not establish multi-node behavior, cover the full degree and
physics matrix, or prove that either complete profile was executed. Machine
CPU, memory, scheduler, and wall-time limits must be reviewed before selecting
a larger slice.
