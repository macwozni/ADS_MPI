# Temporal-convergence instability across meshes and degree pairs

Status: historical pre-fix evidence recorded on 2026-09-29. Core commit
`5e1161b` removed several failures listed below, so their exact values must
not be read as current-HEAD results. Post-fix spot checks still expose a
degree-dependent non-monotone series; full temporal qualification remains
open and the oracle has not been weakened.

The Stage-2 investigation first exposed pathological transient results on a
`2x2x2` mesh. Wider Stage-3 runs show that "even mesh" is not a sufficient
diagnosis: deterministic outliers also occur on `3x3x3`, while some
`4x4x4` series have a correct asymptotic tail. The failure depends on the
combination of mesh, time step, scheme, and test/trial degrees.

## Reproducer

Build the isolated `igrm_l2` adapter from the repository root. Release is used
below only to shorten the run; the first pathological `4x4x4` PR case was also
run with the debug/ASan executable and returned exactly the same errors.

```bash
make -C benchmarking build-igrm_l2 BUILD=release
```

Every command uses `T=0.1`, one MPI rank, one OpenMP thread, a `1x1x1`
process grid, 17 sampling points per axis, and no sample-file output. Use the
MPI launcher associated with the compiler selected by `CONFIG`.

An isolated nonmonotone DG point on the odd mesh is:

```bash
env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close \
  mpiexec -n 1 \
  ./benchmarking/build/release/EXEC/igrm_l2_manufactured \
  dg 0.1 32 3 3 3 4 4 4 3 3 3 1 1 1 17 0
```

It reports `l2_error=3.5498768978062012e-1` and
`linf_error=5.4807048290814064`, although the adjacent `N=16` and `N=64`
runs report L2 errors `5.1522031505851026e-6` and
`3.1575064843927544e-7`. Repeating `N=32` returns the same result.

A coarse-grid failure on the official `4x4x4` temporal mesh is:

```bash
env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close \
  mpiexec -n 1 \
  ./benchmarking/build/release/EXEC/igrm_l2_manufactured \
  pr 0.1 4 4 4 4 4 4 4 3 3 3 1 1 1 17 0
```

It reports `l2_error=2.5014012371112409e8` and
`linf_error=8.6492818523908844e9`. The debug/ASan executable returns the same
values, so this is not an optimization artifact.

The degree dependence can be reproduced by changing only the test degree:

```bash
env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close \
  mpiexec -n 1 \
  ./benchmarking/build/release/EXEC/igrm_l2_manufactured \
  dg 0.1 4 4 4 4 5 5 5 3 3 3 1 1 1 17 0

env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close \
  mpiexec -n 1 \
  ./benchmarking/build/release/EXEC/igrm_l2_manufactured \
  dg 0.1 4 4 4 4 5 5 5 4 4 4 1 1 1 17 0
```

The `(p_test,p_trial)=(5,3)` command returns L2 error about `3.060e2`,
whereas `(5,4)` returns the expected `8.251e-5`. Conversely, PR at `N=4`
returns about `1.110e-1` for `(5,3)` and `1.813e8` for `(5,4)`.

## Observed temporal series

All records below were finite, had `solver_status=0`, ended at the requested
`T`, and passed the exact initial-projection oracle. Thus process success or a
correct initial projection does not certify the evolved field.

For `3x3x3`, `(p_test,p_trial)=(4,3)`:

| scheme | N | L2 error | local L2 order | Linf error | local Linf order |
| --- | ---: | ---: | ---: | ---: | ---: |
| DG | 4 | `8.251576e-5` | - | `3.772257e-4` | - |
| DG | 8 | `2.039529e-5` | `2.016` | `9.326909e-5` | `2.016` |
| DG | 16 | `5.152203e-6` | `1.985` | `3.272828e-5` | `1.511` |
| DG | 32 | `3.549877e-1` | `-16.072` | `5.480705e0` | `-17.353` |
| DG | 64 | `3.157506e-7` | `20.101` | `1.443726e-6` | `21.856` |
| PR | 4 | `2.176948e-4` | - | `1.302970e-3` | - |
| PR | 8 | `9.372403e-5` | `1.216` | `4.628776e-4` | `1.493` |
| PR | 16 | `4.575522e-5` | `1.034` | `2.279649e-4` | `1.022` |
| PR | 32 | `2.261230e-5` | `1.017` | `1.124437e-4` | `1.020` |
| PR | 64 | `1.124119e-5` | `1.008` | `5.580797e-5` | `1.011` |
| BE | 4 | `4.169852e-4` | - | `2.202145e-3` | - |
| BE | 8 | `1.435181e-4` | `1.539` | `7.294717e-4` | `1.594` |
| BE | 16 | `5.778504e-5` | `1.312` | `2.658487e-4` | `1.456` |
| BE | 32 | `2.620820e-5` | `1.141` | `1.077232e-4` | `1.303` |
| BE | 64 | `1.872249e-2` | `-9.481` | `3.670692e-1` | `-11.735` |

For the official `4x4x4` mesh and the same degrees:

| scheme | N | L2 error | local L2 order | Linf error | local Linf order |
| --- | ---: | ---: | ---: | ---: | ---: |
| DG | 4 | `8.251611e-5` | - | `4.067121e-4` | - |
| DG | 8 | `2.561321e-5` | `1.688` | `4.750393e-4` | `-0.224` |
| DG | 16 | `5.071582e-6` | `2.336` | `2.320140e-5` | `4.356` |
| DG | 32 | `1.265771e-6` | `2.002` | `5.787106e-6` | `2.003` |
| DG | 64 | `3.159916e-7` | `2.002` | `1.445023e-6` | `2.002` |
| DG | 128 | `7.869264e-8` | `2.006` | `3.582485e-7` | `2.012` |
| PR | 4 | `2.501401e8` | - | `8.649282e9` | - |
| PR | 8 | `3.737280e-2` | `32.640` | `2.030622e0` | `31.988` |
| PR | 16 | `2.320375e-4` | `7.331` | `6.827116e-3` | `8.216` |
| PR | 32 | `2.257949e-5` | `3.361` | `1.125050e-4` | `5.923` |
| PR | 64 | `1.124023e-5` | `1.006` | `5.585438e-5` | `1.010` |
| PR | 128 | `5.605335e-6` | `1.004` | `2.779671e-5` | `1.007` |
| BE | 4 | `5.477451e-1` | - | `2.074543e1` | - |
| BE | 8 | `1.435762e-4` | `11.897` | `7.325840e-4` | `14.789` |
| BE | 16 | `5.777958e-5` | `1.313` | `2.641202e-4` | `1.472` |
| BE | 32 | `2.620821e-5` | `1.141` | `1.077406e-4` | `1.294` |
| BE | 64 | `1.258296e-5` | `1.059` | `4.752084e-5` | `1.181` |
| BE | 128 | `6.185226e-6` | `1.025` | `2.217765e-5` | `1.099` |

Clean four-level tails demonstrate the intended formal behavior when the
defect is not triggered: DG regression is about `2.00`, while PR and BE
regressions approach `1.00`. The isolated spikes and enormous coarse errors
are not roundoff or solver plateaus: they are not a narrow, flat suffix and
must not be excluded from the oracle as such.

## Consequence for Stage 3

The original `2x2x2` observations remain reproducible and are retained in
[`even-mesh-transient.md`](even-mesh-transient.md), but its original working
diagnosis is superseded by this broader evidence.

Both temporal profiles remain useful for planning and diagnosis:

- `temporal-full` is the required 792-case `4x4x4` matrix;
- `temporal-validation` is the 72-case `3x3x3`, `(4,3)` matrix.

Neither historical run is a passing scientific qualification. Core fix
`5e1161b` invalidated the specific DG/BE failures quoted above, but it did not
by itself qualify either complete profile. The analyzer must retain its
analytic errors, strict monotonicity checks, and order thresholds. Both
profiles need a fresh post-fix run, and any remaining outlier needs its own
production diagnosis before the profile can be claimed green.
