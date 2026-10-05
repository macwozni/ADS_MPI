# ADS MPI benchmarking framework

This directory contains the tracked, reusable benchmark infrastructure.
Generated manifests, logs, fields, and measurements belong under the ignored
`benchmarks/<run-id>/` tree. The pre-existing
`benchmarks/igrm_strong_scaling/` prototype is user data and is never adopted,
overwritten, moved, or cleaned by this framework.

Stage 2 adds a real manufactured transient, a shared Fortran lifecycle, and
executable adapters for `igrm_l2`, `igrm_heat`, and
`pure_diffusion_igrm`. Stage 3 adds frozen temporal runs, verified resume, and
a registered convergence analyzer. Stage 4 adds an independent non-polynomial
manufactured case, mesh-size and degree-convergence workflows, and explicit
temporal-error qualification of every spatial point. Stage 5 adds a reusable
regular-grid field comparator, MPI/OpenMP correctness matrices, and a timing
eligibility gate based on the complete numerical field. Stage 6 adds real
strong scaling with repeated solver-side timings, field-gated statistics, and
portable MPI/OpenMP runtime policy. Expensive runs remain separate from
`make test`. Stage 7 adds per-rank weak and hybrid scaling, a neutral
scheduler launcher template, deterministic plan sharding, and verified shard
result merging. It uses the same timing and complete-field correctness engine;
weak-scaling ratios are not relabeled as strong-scaling speedup.
Stage 8 binds every executed run to explicit build and machine provenance,
classifies failures and retries only transient ones, and adds a numerics-first
A/B comparison. It can compare two completed run directories or execute two
Git refs from isolated detached worktrees without checking out or modifying the
caller's working tree.

## Architecture and extension contract

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

### How to add a benchmark

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

### How to add a problem adapter

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

## Manufactured transients

The registered exact case is `temporal-polynomial` on the unit cube:

```text
q(s)       = s^2 (3 - 2s)
Q(x,y,z)   = q(x) q(y) q(z)
u(x,y,z,t) = exp(-t) Q(x,y,z)
```

Because `q'(0)=q'(1)=0`, the normal flux is zero on every face. The
production weak form is `M u_t + K u = f`, with `K` representing
`-Delta`. Since `q''(s)=6-12s`, the source for
`u_t - Delta(u) = f` is

```text
f(x,y,z,t) = exp(-t) [
    -Q(x,y,z)
    + (12x-6) q(y)q(z)
    + (12y-6) q(x)q(z)
    + (12z-6) q(x)q(y)
]
```

The initial field is `Q`. It is exactly representable when every trial
degree is at least three. The harness performs one mass-only projection at
`t=0`, verifies that state, and then performs exactly `N` physical updates
with `dt=T/N`. It asserts both the computed final time and the ADS state time
equal `T`.

Source evaluation follows the production scheme tables:

| Scheme | source time in physical step `[t_n,t_n+dt]` |
| --- | --- |
| Douglas-Gunn (`dg`) | `t_n + dt/2` |
| split Backward Euler (`be`) | `t_n + dt` |
| cyclic Peaceman-Rachford (`pr`) | `t_n + dt/6`, `t_n + dt/2`, `t_n + 5dt/6` |

The callback time used inside OpenMP RHS assembly is thread-private.

The independent `spatial-cosine` case is used by the `h`, `p`, `validation`,
and `strong` families:

```text
R(x,y,z)   = cos(pi*x) cos(pi*y) cos(pi*z)
u(x,y,z,t) = exp(-t) R(x,y,z)
f(x,y,z,t) = (3*pi^2 - 1) exp(-t) R(x,y,z)
```

Here `Delta(R)=-3*pi^2*R`, so the source has the displayed sign for the
production convention `u_t-Delta(u)=f`. The normal derivative vanishes on all
six faces. Unlike `temporal-polynomial`, this field is not exactly
representable by any finite spline degree. Its initial projection error is
therefore a measured spatial error, not a failed temporal-case initialization
gate.

## Problem adapters

All adapters use the shared initialization, projection, DG/PR/BE wrappers,
measurement, MPI status handling, and cleanup. They are still distinct
benchmark entry points and select the appropriate production path:

| Adapter | benchmark-local responsibility |
| --- | --- |
| `igrm_l2` | Uses the generic production `ComputePointForRHS` path. The public problem performs one solve; this adapter deliberately reuses the initialized ADS state for `N` physical updates to `T`, so the run is a real temporal experiment. |
| `igrm_heat` | Uses the problem's production `heat_igrm_rhs_point` full-RHS callback and preserves its physical-time branch, while replacing the scalar source with the manufactured source. No VTI is written in timed steps. |
| `pure_diffusion_igrm` | Uses the generic production RHS path and accepts independent test/trial degrees in all three axes instead of the public driver's isotropic degree shortcut. |

No production source or public problem CLI is changed. The executables are:

```text
benchmarking/build/<debug|release>/EXEC/igrm_l2_manufactured
benchmarking/build/<debug|release>/EXEC/igrm_heat_manufactured
benchmarking/build/<debug|release>/EXEC/pure_diffusion_igrm_manufactured
```

Their common direct CLI is:

```text
<scheme> <T> <steps> <nx> <ny> <nz>
<ptest-x> <ptest-y> <ptest-z>
<ptrial-x> <ptrial-y> <ptrial-z>
<proc-x> <proc-y> <proc-z>
<sample-points> <write-samples:0|1> [exact-case]
```

Normally use the Python runner, which derives this command from the normalized
case and checks the returned record against it.

## Profiles and normalized cases

Profiles are strict JSON and require no third-party Python packages. Every
profile supplies the family, problem, scheme, exact case, `T`, step count or
`dt`, three-dimensional mesh and degrees, MPI layout, OpenMP threads,
sampling, measurement controls, build profile, and launcher.

Time values are normalized with exact rational arithmetic. A recurring value
such as `0.1/3` is stored as `1/30`, then deliberately converted to binary64
by the executable adapter. `T/dt` must be an exact positive integer.

Registered profiles are:

- `smoke`: the real 3 problems x 3 schemes matrix at `N=4`, `3x3x3`,
  `p_test=4`, `p_trial=3`, MPI 1, OMP 1;
- `smoke-refined`: the matching nine cases at `N=8`;
- `temporal-validation`: 72 cases on `3x3x3` at `(p_test,p_trial)=(4,3)`,
  spanning all three problems, all three schemes, and all eight time levels;
- `temporal-full`: the specified 792-case matrix;
- `h-convergence-smoke`: 54 cases over `2^3,4^3,8^3`;
- `h-convergence-full`: 90 cases over `2^3,4^3,8^3,16^3,32^3`, with
  fixed `(p_test,p_trial)=(3,2)` in every axis;
- `p-convergence-smoke`: 108 cases covering `p_trial=1,2,3` with both
  isotropic enrichments `+1` and `+2`;
- `p-convergence-full`: all 270 isotropic cases for `p_trial=1,...,8`,
  `p_test=p_trial+1`, and the admissible `p_test=p_trial+2 <= 9` cases, on
  the fixed `2^3` mesh; the required `(p_test,p_trial)=(2,1)` is covered by
  the smoke profile and by the post-fix validation recorded in
  [`reproducers/spatial-degree-solver-defects.md`](reproducers/spatial-degree-solver-defects.md);
- `p-anisotropic-smoke`: 108 cases over three cyclic rotations of `(3,4,5)`;
- `p-anisotropic-full`: 216 cases over all six rotations of `(3,4,5)`, each
  with componentwise test enrichment `+1` and `+2`;
- `strong-scaling-smoke`: eight executable strong-scaling cases covering
  DG/PR, two degree pairs, and the `1x1x1`/`2x1x1` resource levels;
- `local-scaling`: 2,160 strong-scaling cases at global mesh `16^3`, covering
  all problems, schemes, 15 isotropic degree pairs, OpenMP `1,2,4,8`, and the
  serial plus X/Y/Z two-rank layouts;
- `cluster-scaling`: the 12,960-case full strong-scaling matrix, adding global
  meshes `32^3,64^3`, the XY/XZ/YZ four-rank layouts, and `2x2x2`.
- `weak-scaling-smoke`: two timed `igrm_l2`/DG configurations at local `4^3`,
  MPI `1x1x1` and `2x1x1`, OMP1, plus one deduplicated serial full-field
  helper for the two-rank global mesh (three cases total);
- `local-weak-scaling`: 6,480 timed configurations over all three problems,
  DG/PR/BE, all 15 degree pairs, local `8^3,12^3,16^3`, process grids
  `1x1x1,2x1x1,2x2x1,2x2x2`, and OMP `1,2,4,8`, plus 1,080 serial
  full-field helpers (7,560 cases total);
- `cluster-weak-scaling`: the same physics, degrees, local work, and OpenMP
  levels over `1x1x1,2x1x1,2x2x1,2x2x2,4x2x2,3x3x3,4x4x2,4x4x4`;
  it has 12,960 timed configurations plus 2,025 deduplicated helpers, for
  14,985 cases total.

`temporal-full` expands

```text
3 problems x 3 schemes x 8 step counts x 11 degree pairs = 792
```

at `T=0.1`, `N=4,...,512`, and fixed `4x4x4` mesh. Planning it works,
but it is not currently a passing scientific qualification. The smaller
`temporal-validation` profile expands

```text
3 problems x 3 schemes x 8 step counts x 1 degree pair = 72
```

on `3x3x3`. It is a diagnostic matrix, not a workaround: wider validation
found deterministic outliers on odd and even meshes and across different
degree pairs. See
[`reproducers/temporal-convergence-instability.md`](reproducers/temporal-convergence-instability.md)
and the original
[`reproducers/even-mesh-transient.md`](reproducers/even-mesh-transient.md).
The Stage-2 smoke profiles retain their valid two-level refinement check; the
oracle was not weakened, and that check is not presented as a full temporal
qualification.

Every spatial profile uses the short observation window `T=0.01` and contains
both `N=128` and `N=256` for each otherwise identical point. The shorter window
keeps evolution error below the projection error over a useful spatial range;
it is not treated as proof of separation. The resulting `dt/dt2` pair is a
candidate qualification pair. Every run writes the common-grid field samples
needed to qualify L2 and Linf independently; the smoke and degree profiles use
`33^3` points and the full h profile uses `65^3`.
Each metric must then satisfy the spatial analyzer's fixed 10% temporal-share
rule. If a point fails, the report requires a finer pair and excludes that
point from spatial orders or degree-decrease calculations.

Each normalized case is encoded as sorted compact JSON. Its ID is
`<problem>-<scheme>-<20 hex digits>`, derived from SHA-256 of that semantic
document. Profile name, run ID, timestamps, output paths, Git state, and JSON
key order do not affect the ID. Duplicate semantic cases are rejected before
filtering.

## Plan, build, and run

From the repository root:

```bash
# No-write planning.
make benchmark-plan
make benchmark-plan BENCHMARK_PROFILE=temporal-full

# Isolated ADS library plus all three adapters.
make benchmark-build BUILD=debug

# Build and run N=4 and N=8, then require improvement for all nine pairs.
make benchmark-smoke BENCHMARK_RUN_ID=stage2-smoke

# Plan the complete and diagnostic temporal matrices without running them.
make benchmark-plan BENCHMARK_PROFILE=temporal-full
make benchmark-plan BENCHMARK_PROFILE=temporal-validation

# Plan Stage-4 spatial matrices without running them.
make benchmark-plan BENCHMARK_PROFILE=h-convergence-full
make benchmark-plan BENCHMARK_PROFILE=p-convergence-full
make benchmark-plan BENCHMARK_PROFILE=p-anisotropic-full

# Execute frozen full workflows; select the matching *-smoke profile for smoke.
make benchmark-h-convergence RUN_ID=stage4-h
make benchmark-p-convergence RUN_ID=stage4-p
make benchmark-p-convergence RUN_ID=stage4-p-anisotropic \
  BENCHMARK_P_CONVERGENCE_PROFILE=p-anisotropic-full

make benchmark-h-convergence RUN_ID=stage4-h-smoke \
  BENCHMARK_H_CONVERGENCE_PROFILE=h-convergence-smoke
make benchmark-p-convergence RUN_ID=stage4-p-smoke \
  BENCHMARK_P_CONVERGENCE_PROFILE=p-convergence-smoke

# Stage-7 weak matrices: planning writes nothing and needs no allocation.
make benchmark-plan BENCHMARK_PROFILE=local-weak-scaling
make benchmark-plan BENCHMARK_PROFILE=cluster-weak-scaling

# Real two-level weak-scaling check on a declared local allocation.
make benchmark-weak RUN_ID=stage7-weak-smoke \
  BENCHMARK_WEAK_PROFILE=weak-scaling-smoke \
  BENCHMARK_PLAN_ARGS='--available-mpi-slots 2 --available-cpu-slots 2'
```

`benchmark-smoke` requires a new base ID. It creates
`benchmarks/stage2-smoke-n4/` and `benchmarks/stage2-smoke-n8/`; existing
run directories are never overwritten. The lower-level equivalents include:

```bash
make -C benchmarking build BUILD=release
make -C benchmarking build-igrm_heat BUILD=debug
make -C benchmarking smoke BENCHMARK_RUN_ID=stage2-smoke
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
  ./benchmarking/build/debug/EXEC/igrm_l2_manufactured \
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
  RUN_ID=stage3-validation \
  BENCHMARK_CONVERGENCE_PROFILE=temporal-validation
```

An existing run ID is never adopted or overwritten by a new run. To continue
an interrupted run, request the same profile and filters explicitly:

```bash
make benchmark-resume \
  RUN_ID=stage3-validation \
  BENCHMARK_RESUME_PROFILE=temporal-validation
```

Transient retries are opt-in. The same public target can, for example, make up
to two additional attempts after each retryable failure:

```bash
make benchmark-resume \
  RUN_ID=stage3-validation \
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
- configured and linked MUMPS, BLAS, LAPACK, ScaLAPACK, ParMETIS, METIS,
  GKlib, and other libraries, including the MUMPS version where discoverable;
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

## Temporal convergence analysis

Analyze a complete frozen run and optionally request a log-log plot:

```bash
make benchmark-analyze RUN_ID=stage3-validation
make benchmark-analyze RUN_ID=stage3-validation \
  BENCHMARK_ANALYZE_ARGS=--plot
```

The analyzer first requires every result named by the frozen manifest and
rejects missing, duplicated, misnamed, or configuration-mismatched records.
It then groups cases that differ only in time resolution. An unfiltered full
run requires all eight levels `N=4,8,...,512`; an explicitly filtered run may
select a common, complete sequence of at least four levels. L2 and Linf are
analyzed independently using errors against the registered analytical
solution, not differences between two numerical solutions. A full report
contains all seven local estimates

```text
log(e_i/e_(i+1)) / log(dt_i/dt_(i+1))
```

and a centered log-log regression over the four finest usable consecutive
points. The acceptance contract is deliberately fixed:

| scheme | theoretical order | accepted regression order |
| --- | ---: | ---: |
| DG | 2 | `[1.70, 2.30]` |
| PR | 1 | `[0.80, 1.20]` |
| BE | 1 | `[0.80, 1.20]` |

Both metrics also require regression `R^2 >= 0.98` and at least a fourfold
reduction of analytical error. Non-finite data, a nonzero solver status,
wrong final time, missing levels, duplicate IDs, or nonmonotonic analytical
error before a credible plateau are hard failures.

A positive local order above four times the scheme's theoretical order is
also a hard metric failure. This deliberately retains catastrophic coarse
solves followed by an apparently clean tail, such as the measured PR defect;
the tail regression alone cannot certify the full series.

A plateau may remove only a finest-level suffix from regression. Its local
orders must stay below `max(0.20, 0.25*p_expected)` and its errors within a
`max/min <= 1.25` band. At least two consecutive plateau transitions are
required, except that one final transition is allowed at an actual roundoff
scale (`<= 1e-10` times the solution norm). The excluded points and their
local orders remain in the report. Isolated spikes and large coarse-grid
errors are therefore never relabeled as plateau.

Analysis writes `analysis/analysis.json` and `analysis/analysis.csv` below the
run directory. `--plot` additionally requests `analysis/convergence.png`; if
matplotlib is unavailable, only the plot is omitted and numerical analysis
still runs. Reanalysis removes any older PNG before deciding whether the new
report has a plot, so a stale image can never describe newer JSON/CSV output.

## Spatial and degree convergence analysis

Run `make benchmark-analyze RUN_ID=<name>` for either an `h` or a `p` run.
The analyzer first enforces the same frozen-manifest and result-integrity
checks as temporal analysis. For each spatial point and each metric it then
subtracts the two sampled error fields. Because both runs have the same exact
case and final time, this is exactly the sampled numerical-field difference:

```text
delta_i = (u_dt - u_exact)_i - (u_dt/2 - u_exact)_i
        = u_dt,i - u_dt/2,i

D_L2   = sqrt(composite_trapezoid_3d(delta_i^2))
D_Linf = max_i abs(delta_i)
c_t    = D_metric / E_metric(dt/2)
```

The CSV header, row count, X-fastest regular-grid order, coordinates, exact
values, pointwise errors, reported sampled Linf, and deterministic field
checksum are all rebound and checked before subtraction. Merely subtracting
the two scalar error norms is
not used: that would be only a lower bound and could falsely pass two very
different fields with equal norms. A point is spatially reliable only when
`c_t <= 0.10`. The accepted value is `E(dt/2)`. A larger or non-finite
indicator marks the metric
`time_dominated/unreliable`, requests a finer `dt` pair, and suppresses every
order or degree-drop estimate involving that point. This same-layout temporal
check is local to Stage 4; cross-layout MPI/OpenMP field equivalence remains a
separate workflow.

`D_L2` and `D_Linf` are discrete common-grid diagnostics. The reported L2
error in the denominator is the independent Gauss-integrated error, while the
reported Linf error and `D_Linf` are maxima on the declared regular grid, not
a proof of the continuous supremum. In particular, `33^3` is a smoke/degree
resolution and `65^3` gives only two sampling intervals per element at the
finest `32^3` h level. Reports preserve `sample_points_per_axis` and the
difference method so a denser study can be distinguished rather than silently
compared as the same Linf experiment.

During Stage 4 the spatial workflow exposed a structurally singular mixed
iGRM system at the coarser PR step and at required low degree. The solver fix
was kept in the separate core commit `5e1161b`; pre-fix commands, root cause,
and post-fix numerical evidence are retained in
[`reproducers/spatial-degree-solver-defects.md`](reproducers/spatial-degree-solver-defects.md).

For `h`, only cases that differ in the element counts are compared. The report
retains all qualified L2 and Linf values, every local
`log(E_i/E_(i+1))/log(h_i/h_(i+1))` estimate, and the regression slope versus
`h`. For isotropic `p`, the two enrichment families are kept separate and the
report records error decrease and adjacent ratios versus the actual trial
degree; it does not invent a fixed algebraic convergence order. A plateau is
accepted only after at least two strict, qualified error decreases, so a flat
sequence from the first degree cannot pass as convergence. Anisotropic degree
vectors remain vectors in configuration, identity, grouping, JSON, and CSV
reports. Their rotation spread is an audit result, not an ordered convergence
gate: the ADI axis order does not justify an arbitrary rotational-equality
tolerance. Solver/roundoff plateaus and time-dominated points or ranges are
retained as auditable data rather than silently discarded.

## Full-field MPI/OpenMP validation

The Stage 5 oracle is implemented in the shared `ads_benchmark.validation`
package and is also used by the spatial analyzer. It parses the canonical
physical-grid artifact, verifies finite values, row count, X-fastest ordering,
regular coordinates, the analytical component, `error=numerical-exact`, the
reported sampled Linf, and the solver checksum. Comparisons then inspect every
named component with

```text
abs(candidate-reference) <= 1e-11 + 1e-10 * max(abs(reference),abs(candidate))
```

for layout equivalence. The report additionally retains composite-trapezoid
L2 and pointwise Linf differences, both field norms, canonical SHA-256
checksums, mismatch counts, the largest absolute difference and its physical
coordinate, and (on failure) the point with the largest normalized tolerance
violation.
A serial `MPI=1, grid=1x1x1, OMP=1` case is the reference for each otherwise
identical problem, scheme, time level, mesh, degree, sample grid, and build.
Every non-reference timing has `timing_valid=false` unless its complete field
passes this gate. At the finest level, the serial and parallel timings are
also ineligible unless the scheme passes the analytical and cross-scheme
accuracy/trend gates. Neither an L2 norm nor a checksum can substitute for
the pointwise decision; the self-test swaps two local values while preserving
the sum and L2 norm and requires the comparator to fail.

The checked profiles use the same decomposition-independent `17^3` physical
grid, spatial-cosine solution, `2^3` element mesh, test/trial degrees `3/2`,
`T=0.01`, and time levels `N=128,256`:

- `validation-smoke`: one problem/scheme, 36 cases;
- `validation-full`: three problems, DG/PR/BE, 324 cases.

Each level includes OMP1 and OMP4 for `1x1x1`, `2x1x1`, `1x2x1`, `1x1x2`,
`2x2x1`, `2x1x2`, `1x2x2`, `2x2x2`, and the uneven `3x2x1` layout. The
ordinary case validator enforces `NP=procx*procy*procz` and at least one owned
trial-space DOF per axis. Uneven quotient/remainder partitions are supported,
so divisibility is deliberately not required; the library's arbitrary-width
coefficient-halo exchange is exercised by the matrix.

At the finest planned `dt`, each serial scheme must be pointwise within an
absolute `0.05` of the analytical solution (5% of this manufactured field's
unit amplitude, with relative term zero). DG/PR/BE pairs must be pointwise
within absolute `0.01` of one another (1% of unit amplitude, again with
relative term zero). Their L2 and Linf differences are reported at both time
levels; from coarse to fine neither may grow by more than 0.1%, and at least
one must contract by at least 0.1% (an already bitwise-identical final pair
also passes). Thus a coarse difference is recorded rather than mistaken for
a parallel failure, while final agreement and a measurable refinement trend
remain hard gates. These fixed limits are documented policy, not values
enlarged after observing a failed run.

Run and analyze the complete matrix with one public target, declaring an
allocation that can provide its largest 8-rank/OMP4 case:

```bash
make benchmark-validate RUN_ID=stage5-full \
  BENCHMARK_PLAN_ARGS='--available-mpi-slots 8 --available-cpu-slots 32'
```

`make -C benchmarking plan BENCHMARK_PROFILE=validation-full` remains the
portable full-matrix dry-run when that allocation is unavailable.

For a declared 24-CPU/6-rank allocation, a representative Y/Z/uneven subset
can be selected without changing the profile:

```bash
make benchmark-validate RUN_ID=stage5-subset \
  BENCHMARK_PLAN_ARGS='--problem igrm_l2 \
    --mpi-grid 1,1,1 --mpi-grid 1,2,1 --mpi-grid 1,1,2 --mpi-grid 3,2,1 \
    --available-mpi-slots 6 --available-cpu-slots 24'
```

The target writes `analysis.json` and a readable `analysis.csv` with
`analytic-level`, `parallel-variant`, and `scheme-pair` rows. Each difference
row names the problem, scheme or scheme pair, decomposition, OpenMP size,
coordinates, values, tolerances, and maximum deviation. Scheme-pair rows also
carry the aggregate agreement/trend result and the coarse/final L2 and Linf
differences, so a trend failure is visible without consulting JSON.
Validation has no convergence plot; `--plot` is explicitly omitted rather
than producing a misleading image. VTI remains optional and outside the
measured physical-step interval; validation uses normalized numerical samples
and never compares VTI bytes or metadata ordering.

## Strong scaling

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

## Weak and hybrid scaling

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

Weak cases use the same Stage-6 timing method: two validated warmups and seven
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

## Sharded benchmark runs

Large weak (or other) frozen plans can be split without editing their
configuration. `problem-scheme-degree` creates one deterministic shard per
problem/scheme/test-degree/trial-degree group. `index` assigns the sorted case
IDs round-robin to an explicit shard count:

```bash
# Public root-Make forms.
make benchmark-shard-plan RUN_ID=weak-parent-grouped \
  BENCHMARK_SHARD_PROFILE=cluster-weak-scaling \
  BENCHMARK_SHARD_STRATEGY=problem-scheme-degree

make benchmark-shard-plan RUN_ID=weak-parent \
  BENCHMARK_SHARD_PROFILE=cluster-weak-scaling \
  BENCHMARK_SHARD_STRATEGY=index BENCHMARK_SHARD_COUNT=3
```

The lower-level CLI equivalents make every frozen input explicit:

```bash
cd benchmarking

# Create the immutable parent manifest and grouped shard envelopes.
python3 -m ads_benchmark.cli shard \
  --repository-root .. --config-dir configs \
  --profile cluster-weak-scaling --run-id weak-parent-grouped \
  --strategy problem-scheme-degree

# Alternative used below: exactly three index shards.
python3 -m ads_benchmark.cli shard \
  --repository-root .. --config-dir configs \
  --profile cluster-weak-scaling --run-id weak-parent \
  --strategy index --shard-count 3
```

The parent owns `benchmarks/<parent>/manifest.json` and numbered immutable
envelopes below `benchmarks/<parent>/shards/`. Each envelope binds the complete
parent metadata, repository/source identity, full expected case digest, shard
index/count, and subset digest. Generation refuses an existing parent.

After building the release adapters, execute each shard into its own run ID.
The command below illustrates shard zero; allocation values must cover that
shard's largest case:

```bash
make benchmark-run-shard \
  PARENT_RUN_ID=weak-parent SHARD_INDEX=0 \
  RUN_ID=weak-parent-part-000 \
  BENCHMARK_LAUNCHER_TEMPLATE='srun --ntasks={ranks} --cpus-per-task={threads} {payload}' \
  BENCHMARK_PLAN_ARGS='--available-mpi-slots 64 --available-cpu-slots 512'
```

The public target builds the release adapters. Its lower-level equivalent,
when already inside `benchmarking/`, is:

```bash
make build BUILD=release

python3 -m ads_benchmark.cli run-shard \
  --repository-root .. \
  --parent-run-id weak-parent --shard-index 0 \
  --run-id weak-parent-part-000 \
  --available-mpi-slots 64 --available-cpu-slots 512 \
  --launcher-template \
    'srun --ntasks={ranks} --cpus-per-task={threads} {payload}'
```

`run-shard --resume` verifies and resumes that isolated shard run by the same
rules as an ordinary frozen run. Once every shard has completed, list every
child run explicitly and merge into the parent analysis:

```bash
make benchmark-merge-shards PARENT_RUN_ID=weak-parent \
  SHARD_RUNS='weak-parent-part-000 weak-parent-part-001 weak-parent-part-002' \
  BENCHMARK_ANALYZE_ARGS=--plot
```

The direct equivalent is:

```bash
python3 -m ads_benchmark.cli merge-shards \
  --repository-root .. --parent-run-id weak-parent \
  --shard-run weak-parent-part-000 \
  --shard-run weak-parent-part-001 \
  --shard-run weak-parent-part-002 \
  --plot
```

This merge matches the three-shard index parent created above. Merge refuses a
missing shard/index or case, a duplicated case or child contribution,
conflicting case content, incompatible config hashes, and incompatible
repository/source identities. It analyzes only after the contributed child
manifests reconstruct the exact parent plan; it does not silently accept a
partial report. Generated results and shard manifests remain ignored runtime
data and are not committed with the framework.

## A/B numerical and performance comparison

Compare two already completed run directories through the public root target:

```bash
make benchmark-compare \
  BENCHMARK_BASELINE_RESULTS=benchmarks/baseline-run \
  BENCHMARK_CANDIDATE_RESULTS=benchmarks/candidate-run
```

The loader applies the normal manifest, status, command, log, result, and
adapter verification to both runs. It then requires identical configuration
hashes, case-ID sets, and per-case normalized configurations, plus compatible
machine, toolchain, external-library, launcher, topology, binding, and OpenMP
execution provenance. Repository/build-root paths are normalized so isolated
worktrees can be compared; the source refs themselves may differ.

The decision is deliberately numerics-first. Every candidate field must pass
the complete-field comparison against its baseline counterpart before any
timing is classified. One numerical mismatch blocks timing for the whole
comparison. For every eligible case the report retains each side's sample
count, minimum, maximum, median, MAD, relative MAD, and range, then computes

```text
median_ratio = candidate_median / baseline_median
```

The default regression threshold is 5%, so a regression requires
`median_ratio > 1.05`. The default minimum is five measured samples per side;
the policy rejects values below three. Default complete-field tolerances are
absolute `1e-11` and relative `1e-10`. Too few samples, a timing already marked
unreliable, malformed timing, or missing legacy execution provenance produces
`inconclusive`, never a green result or a regression. Legacy runs can therefore
still be inspected numerically without pretending their timings are
reproducible. Incompatible configuration or present-but-different execution
provenance is reported as `incompatible`. Only
`no-regression-detected` returns success; regression, numerical mismatch,
incompatibility, and inconclusive evidence return a nonzero verdict.

The default command prints a summary. Optional versioned JSON and flat CSV,
as well as a different threshold, use the same public Make target:

```bash
make benchmark-compare \
  BENCHMARK_BASELINE_RESULTS=benchmarks/baseline-run \
  BENCHMARK_CANDIDATE_RESULTS=benchmarks/candidate-run \
  BENCHMARK_COMPARE_ARGS='--regression-threshold 0.03 --minimum-samples 5 --output-json /tmp/ads-ab.json --output-csv /tmp/ads-ab.csv'
```

To build and run two refs on the same host without touching a dirty main
worktree, supply refs instead of result directories:

```bash
make benchmark-compare \
  RUN_ID=stage8-ab \
  BENCHMARK_BASELINE_REF=HEAD~1 \
  BENCHMARK_CANDIDATE_REF=HEAD \
  BENCHMARK_COMPARE_PROFILE=strong-scaling-smoke \
  BENCHMARK_PLAN_ARGS='--available-mpi-slots 2 --available-cpu-slots 2' \
  BENCHMARK_COMPARE_ARGS='--max-retries 1'
```

The ref orchestrator resolves both refs once to full commit SHAs, creates two
distinct detached marker-owned worktrees and separate isolated build roots,
and never checks out the caller's tree. Matched cases use a deterministic
case-level `AB, BA, AB, ...` order to reduce monotonic machine drift. This is
not sample-level alternation: each side of a case still performs all of its own
warmups followed by all of its measured samples. The two exclusive result runs
are `benchmarks/stage8-ab-baseline/` and
`benchmarks/stage8-ab-candidate/` for the example above. Both receive the
versioned schedule; the candidate run receives the JSON/CSV comparison report.

After both sides execute and the report is published, the orchestrator removes
only clean, marker-verified worktrees through `git worktree remove`, even when
the comparison verdict is regression or inconclusive. An execution error or
other exception before report publication, and any modified worktree, is
retained for diagnosis rather than force-deleted. Set
`BENCHMARK_WORKSPACE_PARENT` to choose the temporary parent; otherwise the
system temporary directory is used. The parent must be outside the benchmark
repository so a retained diagnostic workspace cannot silently become an
untracked change in the caller's tree.

For a retained failure, first inspect the printed workspace and both worktrees.
When they are clean and no longer needed, remove the two registered paths with
`git worktree remove <workspace>/baseline` and
`git worktree remove <workspace>/candidate`; then inspect the ownership marker,
unlink that one marker, and remove the now-empty workspace with `rmdir`. Never
use `git worktree remove --force`, `git worktree prune`, or a recursive delete
as a substitute for resolving modified/unowned contents.

The comparator reports MAD and range but deliberately applies a transparent
median-ratio threshold rather than a statistical significance test. Very short
or noisy measurements can therefore cross a tight threshold even for identical
builds; choose a duration and sample count appropriate to the machine and treat
the recorded dispersion as part of the verdict review.

Result-directory inputs and ref inputs are mutually exclusive, and both sides
of the selected mode are required. Use normal profile filters in
`BENCHMARK_PLAN_ARGS` for a small local slice; a full cluster profile remains
a manual allocation-aware workflow.

## Measurement and machine-readable results

Legacy convergence and validation cases retain the exact single-invocation
result schema (`warmups=0`, `samples=1`). Strong and weak cases use the same
executor, adapter, launcher, and result store but add the repeated timing arrays,
reliability threshold/flag, and effective OpenMP environment under `timing`.
Any failed, timed-out, unparsable, or domain-invalid repetition fails the
whole case; resume accepts only a complete, strictly revalidated aggregate.
The local and cluster profiles use a 0.05-second minimum solver-side sample;
the small executable smoke profile uses 0.001 seconds and still reports any
shorter sample as unreliable.

The benchmark harness computes the L2 error and solution norm with an
independent Gauss rule of `min(10,p_trial+3)` points per axis. It deliberately
does not reuse the trial-space assembly rule: doing so can alias the
non-polynomial cosine field and report a false zero error at low degree. Ten
points is the library's supported limit and exactly integrates a squared
degree-nine spline. For Linf and field fingerprinting, rank zero reconstructs
the global coefficients and evaluates an endpoint-inclusive regular grid of
`sample_points_per_axis^3` points, with X varying fastest. Compensated sums
produce a deterministic `field_checksum`; Linf and checksum are then
broadcast. Sampling and optional CSV output occur outside the measured
physical-step interval.

With `write_samples=1`, the case directory also receives
`field_samples.csv` with columns:

```text
x,y,z,numerical,exact,error
```

Spatial analysis requires this artifact for both members of every `dt/dt2`
pair. It is opened relative to the verified case directory without following
symlinks, must be a regular UTF-8 file, and is bounded to 64 MiB before it is
parsed. Profile validation therefore caps `points_per_axis` at 76 whenever
`write_samples=true`; the fixed-width CSV for `77^3` samples cannot fit under
that storage contract. Unwritten sampling retains the executable's general
limit of 257 points per axis. Run loading retains only a lazy safe reference,
and parsed errors use a
compact binary64 array; a full profile therefore never keeps all CSV texts or
per-sample Python objects in memory. The artifact is not copied into
`result.json` or the final report.

Each successful executable writes exactly one stdout line beginning with
`ADS_BENCHMARK_RESULT `, followed by strict JSON. Its domain fields are:

```text
schema_version, kind, exact_case, problem, scheme,
requested_final_time, actual_final_time, time_step, steps,
initial_l2_error, initial_linf_error,
l2_error, linf_error, solution_l2_norm, field_checksum,
sample_points_per_axis, field_samples_written,
physical_step_wall_seconds, solver_status
```

The parser rejects duplicate/missing/extra keys, non-finite values, a nonzero
solver status, multiple tagged records, and any mismatch with the planned
case. The stored outer result adds the stable `case_id`, complete normalized
`configuration`, process `timing.wall_seconds`, and the parsed
`domain_result`.

## Numerical gates

Before timing, every executable requires finite initial metrics. The exactly
representable `temporal-polynomial` additionally requires
`initial_l2_error <= 1e-10`, `initial_linf_error <= 1e-10`, and agreement of
the projected solution norm with `(13/35)^(3/2)` within `1e-10`. The
non-polynomial spatial case retains its finite initial projection errors
instead of applying that exact-representation gate. Failures, MPI/MUMPS
errors, NaN, or Inf produce a nonzero process result.

The smoke oracle loads exactly the nine passed `N=4` results and their nine
`N=8` counterparts, rechecks the initial errors, and requires a strict L2
decrease for every problem/scheme pair. This is a refinement sanity check, not
yet a formal order estimate.

## Result layout and safety

A written or executed plan owns:

```text
benchmarks/<run-id>/manifest.json
benchmarks/<run-id>/execution.json                    # executed runs
benchmarks/<run-id>/cases/<case-id>/status.json
benchmarks/<run-id>/cases/<case-id>/result.json
benchmarks/<run-id>/cases/<case-id>/stdout.log
benchmarks/<run-id>/cases/<case-id>/stderr.log
benchmarks/<run-id>/cases/<case-id>/field_samples.csv  # only when requested
benchmarks/<run-id>/analysis/analysis.json
benchmarks/<run-id>/analysis/analysis.csv
benchmarks/<run-id>/analysis/convergence.png           # only with --plot
benchmarks/<run-id>/analysis/strong-scaling.png         # strong --plot
benchmarks/<run-id>/analysis/weak-scaling.png           # weak --plot
benchmarks/<ab-id>-baseline/analysis/ab-schedule.json    # ref A/B workflow
benchmarks/<ab-id>-candidate/analysis/ab-schedule.json   # ref A/B workflow
benchmarks/<ab-id>-candidate/analysis/comparison.json    # ref A/B report
benchmarks/<ab-id>-candidate/analysis/comparison.csv     # ref A/B report
benchmarks/<parent-run-id>/shards/shard-000000.json     # sharded plans
```

The manifest contains the schema version, full expanded configuration and case
IDs, configuration hash, Git commit, dirty-tree flag, and a SHA-256 content
fingerprint of tracked and nonignored untracked worktree state. Ignored
benchmark/build output is excluded. Manifest, execution, result, and comparison
documents carry a `schema_version` and a `kind`; status documents carry their
own `schema_version` plus the exact case identity and state. Strict input readers
for manifests, execution records, statuses, and results reject unknown, missing,
or inconsistent identity fields. Directory-mode A/B outputs are
written only when explicit `--output-json`/`--output-csv` paths are passed
through `BENCHMARK_COMPARE_ARGS`.

Writes are atomic and exclusive.
Descriptor-relative operations reject traversal, symlink swaps,
sibling-prefix tricks, and adoption of unrelated directories, including the
legacy prototype.

The entire `benchmarks/` tree is ignored runtime data, not versioned source.
Do not stage manifests, fields, logs, reports, or comparison runs. Neither
`make clean` nor `make clean-benchmark-build` deletes a run; the latter removes
only marker-owned framework caches and isolated builds. There is deliberately
no broad results-clean target and no `rm -rf benchmarks` workflow. Remove an
individual run manually only after identifying that exact run directory and
deciding its data is no longer needed.

## Benchmark framework tests

```bash
make benchmark-self-test
```

All test sources and helpers live under `tests/benchmarking/`; the
`benchmarking/` tree contains only benchmark implementation, profiles,
adapters, documentation, and reproducers. The suite is also registered in the
ordinary test hierarchy, so `make test` runs these contract tests without
executing a numerical benchmark profile.

The dependency-free suite covers the 792-case temporal plan, all six spatial
profiles, both validation matrices, both strong-scaling matrices and their
exact case counts, all three weak profiles and their measurement/helper
counts, the final 792/12,960/14,985 configuration audit across all three
problems, DG/PR/BE, eight temporal levels, 11 temporal and 15 scaling degree
pairs, strong/weak, X/Y/Z decomposition, and OMP `1,2,4,8`, plus case-ID
collision and controlled-run-path checks, execution provenance and explicit
unavailable observations,
resume build identity, failure classification and bounded retry histories,
numerics-first A/B policy, detached-worktree ownership and alternating
case scheduling, weak-efficiency and same-global-mesh field gates, launcher-template
expansion and resource rejection, deterministic index/group sharding, and
missing/duplicate/conflicting shard rejection, stable
vector-preserving identifiers,
filters and validation, exact-time conversion, safe storage, process timeout
and termination, frozen-manifest decoding, resume compatibility and verified
skip/retry behavior, strict tagged-result parsing, planned/result binding, the
registered manufactured adapters, smoke-refinement verification, synthetic
temporal, spatial, full-field, strong-scaling, and weak-scaling series, raw
repeated samples, MAD/speedup/efficiency, short-region rejection, local perturbations, separation
failures, plateau and corrupted convergence inputs, fake-adapter extensibility,
and marker-guarded Make cleanup. These tests exercise framework contracts
without launching the costly full numerical matrices. The ordinary
`make test` target does not execute full benchmark profiles; real benchmark
runs remain explicit, and the separate short `make test-performance` gate is
not a substitute for them.
