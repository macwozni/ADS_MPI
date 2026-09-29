# ADS MPI benchmarking framework

This directory contains the tracked, reusable benchmark infrastructure.
Generated manifests, logs, fields, and measurements belong under the ignored
`benchmarks/<run-id>/` tree. The pre-existing
`benchmarks/igrm_strong_scaling/` prototype is user data and is never adopted,
overwritten, moved, or cleaned by this framework.

Stage 2 adds a real manufactured transient, a shared Fortran lifecycle, and
executable adapters for `igrm_l2`, `igrm_heat`, and
`pure_diffusion_igrm`. Expensive runs remain separate from `make test`.

## Architecture and extension contract

```text
configs/*.json
       |
       v
catalog.py --> neutral planner/executor/storage --> benchmarks/<run-id>/
       |                         |
       |                         +--> MPI launcher
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
- launcher.

The generic executor owns the case working directory, environment, OpenMP
settings, process group, timeout, logs, status transitions, and common result
envelope. A Python problem adapter validates one planned case, constructs the
payload argv, parses the tagged solver record, and binds it back to the exact
planned configuration.

The Fortran side follows the same separation. A registered
`BenchmarkAdapter` supplies a manufactured-case descriptor and procedures
for initialize, initial projection, one physical step, measurement, and
cleanup. `benchmark_harness.F90` owns their ordering and the `1..N` time
loop. The three main programs only register one adapter and invoke that shared
harness.

## Manufactured transient

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
<sample-points> <write-samples:0|1>
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
- `temporal-full`: the specified 792-case matrix;
- `local-scaling` and `cluster-scaling`: planning presets whose complete
  scientific scaling workflows belong to later stages.

`temporal-full` expands

```text
3 problems x 3 schemes x 8 step counts x 11 degree pairs = 792
```

at `T=0.1`, `N=4,...,512`, and fixed `4x4x4` mesh. Planning it works,
but executing it as a scientific convergence run is currently blocked by the
known even-mesh transient defect documented in
[`reproducers/even-mesh-transient.md`](reproducers/even-mesh-transient.md).
The Stage-2 smoke profiles intentionally use the verified odd `3x3x3` mesh;
the oracle was not weakened.

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
  dg 0.1 4 3 3 3 4 4 4 3 3 3 1 1 1 17 0
```

Repeatable planner/runner filters include `--problem`, `--scheme`,
`--degree-pair`, `--mesh`, `--mpi-grid`, `--mpi-ranks`, and `--omp`.
Unknown values and an empty result fail explicitly. Validation includes
`p_test > p_trial`, trial degree at least three for this exact case,
maximum degree nine, positive dimensions, `NP=proc-x*proc-y*proc-z`, a
distributable process grid, and positive runtime controls.

The isolated build lives only under `benchmarking/build/<profile>/` and uses
the selected repository `CONFIG`, compiler, and libraries. No MPI or MUMPS
path is hardcoded. `make clean-benchmark-build` removes only marker-owned
benchmark build/cache content, never benchmark results.

## Measurement and machine-readable results

`NormL2` computes the final L2 error against the exact solution and the
solution L2 norm using the production quadrature. For Linf and field
fingerprinting, rank zero reconstructs the global coefficients and evaluates
an endpoint-inclusive regular grid of
`sample_points_per_axis^3` points, with X varying fastest. Compensated sums
produce a deterministic `field_checksum`; Linf and checksum are then
broadcast. Sampling and optional CSV output occur outside the measured
physical-step interval.

With `write_samples=1`, the case directory also receives
`field_samples.csv` with columns:

```text
x,y,z,numerical,exact,error
```

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

Before timing, every executable requires finite initial metrics,
`initial_l2_error <= 1e-10`, `initial_linf_error <= 1e-10`, and agreement
of the projected solution norm with `(13/35)^(3/2)` within `1e-10`.
Failures, MPI/MUMPS errors, NaN, or Inf produce a nonzero process result.

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
```

The manifest contains the schema version, full expanded configuration and case
IDs, configuration hash, Git commit, and dirty-tree flag. Writes are atomic
and exclusive. Descriptor-relative operations reject traversal, symlink
swaps, sibling-prefix tricks, and adoption of unrelated directories, including
the legacy prototype.

## Self-tests

```bash
make benchmark-self-test
```

The dependency-free suite covers the 792-case plan, stable identifiers,
filters and validation, exact-time conversion, safe storage, process timeout
and termination, strict tagged-result parsing, planned/result binding, the
registered manufactured adapters, smoke-refinement verification, fake-adapter
extensibility, and marker-guarded Make cleanup. These tests exercise framework
contracts without adding costly numerical benchmark runs to `make test`.
