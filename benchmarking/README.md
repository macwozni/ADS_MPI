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
temporal-error qualification of every spatial point. Expensive runs remain
separate from `make test`.

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

The Fortran side follows the same separation. A registered
`BenchmarkAdapter` supplies a manufactured-case descriptor and procedures
for initialize, initial projection, one physical step, measurement, and
cleanup. `benchmark_harness.F90` owns their ordering and the `1..N` time
loop. The three main programs only register one adapter and invoke that shared
harness.

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

The independent `spatial-cosine` case is used only by the `h` and `p`
families:

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
- `local-scaling` and `cluster-scaling`: planning presets whose complete
  scientific scaling workflows belong to later stages.

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
`--steps`.
Unknown values and an empty result fail explicitly. Validation includes
`p_test > p_trial`, trial degree at least three for `temporal-polynomial` and
at least one for `spatial-cosine`,
maximum degree nine, positive dimensions, `NP=proc-x*proc-y*proc-z`, a
distributable process grid, and positive runtime controls.

The isolated build lives only under `benchmarking/build/<profile>/` and uses
the selected repository `CONFIG`, compiler, and libraries. No MPI or MUMPS
path is hardcoded. `make clean-benchmark-build` removes only marker-owned
benchmark build/cache content, never benchmark results.

## Frozen convergence runner and resume

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

The same neutral resume variable applies to spatial runs, for example
`BENCHMARK_RESUME_PROFILE=h-convergence-full` or
`BENCHMARK_RESUME_PROFILE=p-convergence-full`. For compatibility,
`BENCHMARK_RESUME_PROFILE` defaults to `BENCHMARK_CONVERGENCE_PROFILE` when it
is not set explicitly.

Resume requires the current profile expansion, filters, commit SHA, dirty
state, and content fingerprint of the nonignored worktree to match the frozen
manifest. Thus two different dirty source trees are not treated as compatible.
It holds an exclusive execution lock and skips a case only when `status.json`,
`result.json`, both logs, the complete normalized configuration, and the
adapter's domain result all revalidate. The tagged stdout is parsed again by
the registered adapter and must reproduce the saved domain result exactly.
Missing, failed, timed-out, incomplete, or tampered cases are retried; a
different configuration or source state is refused.

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

## Measurement and machine-readable results

Stage 3 execution supports exactly one measured invocation per case
(`warmups=0`, `samples=1`). Profiles with other values remain valid for
planning and dry-run, but `run` and `resume` reject them during execution
preflight until repetition support is added with the scaling workflow. This
prevents a manifest from claiming measurements the runner did not perform.

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
benchmarks/<run-id>/cases/<case-id>/status.json
benchmarks/<run-id>/cases/<case-id>/result.json
benchmarks/<run-id>/cases/<case-id>/stdout.log
benchmarks/<run-id>/cases/<case-id>/stderr.log
benchmarks/<run-id>/cases/<case-id>/field_samples.csv  # only when requested
benchmarks/<run-id>/analysis/analysis.json
benchmarks/<run-id>/analysis/analysis.csv
benchmarks/<run-id>/analysis/convergence.png           # only with --plot
```

The manifest contains the schema version, full expanded configuration and case
IDs, configuration hash, Git commit, dirty-tree flag, and a SHA-256 content
fingerprint of tracked and nonignored untracked worktree state. Ignored
benchmark/build output is excluded. Writes are atomic and exclusive.
Descriptor-relative operations reject traversal, symlink swaps,
sibling-prefix tricks, and adoption of unrelated directories, including the
legacy prototype.

## Self-tests

```bash
make benchmark-self-test
```

The dependency-free suite covers the 792-case temporal plan, all six spatial
profiles and their exact case counts, stable vector-preserving identifiers,
filters and validation, exact-time conversion, safe storage, process timeout
and termination, frozen-manifest decoding, resume compatibility and verified
skip/retry behavior, strict tagged-result parsing, planned/result binding, the
registered manufactured adapters, smoke-refinement verification, synthetic
temporal and spatial series, separation failures, plateau and corrupted
convergence inputs, fake-adapter extensibility, and marker-guarded Make
cleanup. These tests exercise framework contracts without adding costly
numerical benchmark runs to `make test`.
