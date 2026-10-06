# Temporal-convergence instability across meshes and degree pairs

Status: fixed in two production-core steps. Commit `5e1161b` removed the
mixed-system singularity behind the original DG/BE/PR outliers. A remaining
PR-only series was then traced to the conditional stability of cyclic PR in
three dimensions and removed by the stabilized `theta=2/3` correction. The
historical tables remain below as evidence; they are not current-HEAD results.
The oracle and its monotonicity/order thresholds were not weakened.

The Stage-2 investigation first exposed pathological transient results on a
`2x2x2` mesh. Wider Stage-3 runs show that "even mesh" is not a sufficient
diagnosis: deterministic outliers also occur on `3x3x3`, while some
`4x4x4` series have a correct asymptotic tail. The failure depends on the
combination of mesh, time step, scheme, and test/trial degrees.

## Residual PR defect after `5e1161b`

A raw spot check on base commit `6d06e82` isolated the remaining defect to
`igrm_l2`, mesh `4x4x4`, `(p_test,p_trial)=(5,4)`, and PR. Every process
returned zero, the final time was exactly `0.1`, and the initial L2 error was
`3.5043802835924535e-15`, but the temporal error increased from `N=4` to
`N=8` before collapsing at `N=16`:

| N | L2 error | Linf error |
| ---: | ---: | ---: |
| 4 | `1.3680492566742378e-2` | `1.4888076539493222e-1` |
| 8 | `1.9124019292174381e-2` | `2.1878176151757733e-1` |
| 16 | `4.6118629766672209e-5` | `2.6224943194563810e-4` |

These values came from direct executable invocations, not a frozen framework
run. The exact command, with `N` replaced by `4`, `8`, or `16`, was:

```bash
env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close \
  mpiexec -n 1 \
  ./benchmarking/build/release/EXEC/igrm_l2_manufactured \
  pr 0.1 N 4 4 4 5 5 5 4 4 4 1 1 1 17 0 temporal-polynomial
```

For mass-normalized diffusion eigenvalues, the former cyclic table had the
modal amplification

```text
product_i (1 - q_j - q_k) / (1 + q_i),  q_i = dt*lambda_i/3.
```

For equal stiff modes this tends to `-8`, so a successful directional solve
could still amplify its input. This explains both the degree dependence and
the counterintuitive failure at an intermediate time resolution.

## Verification after the PR correction

The public `pr` selector now uses an equilibrium-preserving Douglas correction
with `theta=2/3`. Its amplification is

```text
1 - sum(z_i) / product_i(1 + theta*z_i),  z_i = dt*lambda_i,
```

which is contractive for nonnegative diffusion eigenvalues. Repeating the same
case on the correction worktree on 2026-10-06 produced a strictly decreasing
four-level series:

| N | L2 error | Linf error | solver status |
| ---: | ---: | ---: | ---: |
| 4 | `1.7527725990903214e-4` | `9.1521757964363459e-4` | 0 |
| 8 | `5.4771014954474597e-5` | `2.8385157178889564e-4` | 0 |
| 16 | `2.0510879724187206e-5` | `9.8379305864204625e-5` | 0 |
| 32 | `8.9560854801929062e-6` | `3.8312769803261482e-5` | 0 |

These are direct spot checks of the original failing slice. The complete
post-fix gates were subsequently executed from frozen manifests as recorded
below.

## Frozen post-fix gates from 2026-10-06

Both authoritative runs used GCC/GFortran, MPICH, and an isolated MUMPS build
with the bundled PORD ordering only; neither Intel nor ParMETIS was used.

The validation run
`temporal-validation-post-pr-pord-20261006` finished with `72/72` passed
processes. Its analyzer accepted all nine problem/scheme series (`9/9`). This
confirms the corrected PR behavior on the complete diagnostic profile without
relaxing any oracle.

The full gate used parent run
`temporal-full-post-pr-pord-12c-20261006` and twelve deterministic index shards
named `...-part-000` through `...-part-011`. All 792 cases reached a terminal
state:

| state / failure kind | count |
| --- | ---: |
| passed | 653 |
| failed / `numerical` | 72 |
| `timeout` | 67 |
| failed / `mpi` | 0 |

Every numerical failure was the same independent initial-projection defect:
return code 1, test/trial degrees `(9,9,9)/(8,8,8)`, and stdout beginning
`initial manufactured state failed its oracle`. There was exactly one such
case for every problem, scheme, and `N=4,...,512`. Every timeout expired at the
configured 1800-second limit and was classified as `timeout`; no failed case
violated these invariants.

The full result is therefore evidence that the original PR reproducer is fixed,
but it is not a passing full scientific qualification. The official shard
merge deliberately requires every case to pass, so it produced no partial
analysis report. The `(9,8)` projection defect and the runtime envelope must be
handled separately before a fresh `temporal-full` run can be called green.

## Historical pre-fix reproducer

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

## Historical consequence for Stage 3

The original `2x2x2` observations remain reproducible and are retained in
[`even-mesh-transient.md`](even-mesh-transient.md), but its original working
diagnosis is superseded by this broader evidence.

Both temporal profiles remain useful for planning and diagnosis:

- `temporal-full` is the required 792-case `4x4x4` matrix;
- `temporal-validation` is the 72-case `3x3x3`, `(4,3)` matrix.

Neither historical run is a passing scientific qualification. Core commit
`5e1161b` invalidated the specific DG/BE failures quoted above, and the later
stabilized PR correction removed the residual raw outlier without changing the
analyzer. The fresh post-fix validation profile is green; the full profile has
now been executed but remains red for the separate projection-oracle and
timeout reasons recorded above.
