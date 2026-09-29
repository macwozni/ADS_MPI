# Even-mesh manufactured transient instability

Status: open core defect, observed during Stage-2 validation on 2026-09-29.
The benchmark oracle is unchanged, and the production core is not fixed in the
Stage-2 benchmark commit.

## Reproducer

Build the isolated debug adapters from the repository root:

```bash
make benchmark-build BUILD=debug
```

The following commands use one MPI rank, one OpenMP thread, the
`igrm_l2` benchmark adapter, `T=0.1`, a `2x2x2` element mesh,
`p_test=(4,4,4)`, `p_trial=(3,3,3)`, a `1x1x1` process grid, 17 sample
points per axis, and no CSV output:

```bash
env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close \
  mpiexec -n 1 \
  ./benchmarking/build/debug/EXEC/igrm_l2_manufactured \
  dg 0.1 4 2 2 2 4 4 4 3 3 3 1 1 1 17 0

env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close \
  mpiexec -n 1 \
  ./benchmarking/build/debug/EXEC/igrm_l2_manufactured \
  dg 0.1 8 2 2 2 4 4 4 3 3 3 1 1 1 17 0

env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close \
  mpiexec -n 1 \
  ./benchmarking/build/debug/EXEC/igrm_l2_manufactured \
  pr 0.1 4 2 2 2 4 4 4 3 3 3 1 1 1 17 0

env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close \
  mpiexec -n 1 \
  ./benchmarking/build/debug/EXEC/igrm_l2_manufactured \
  pr 0.1 8 2 2 2 4 4 4 3 3 3 1 1 1 17 0

env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close \
  mpiexec -n 1 \
  ./benchmarking/build/debug/EXEC/igrm_l2_manufactured \
  pr 0.1 16 2 2 2 4 4 4 3 3 3 1 1 1 17 0

env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close \
  mpiexec -n 1 \
  ./benchmarking/build/debug/EXEC/igrm_l2_manufactured \
  be 0.1 4 2 2 2 4 4 4 3 3 3 1 1 1 17 0

env OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close \
  mpiexec -n 1 \
  ./benchmarking/build/debug/EXEC/igrm_l2_manufactured \
  be 0.1 8 2 2 2 4 4 4 3 3 3 1 1 1 17 0
```

Use the MPI launcher associated with the compiler selected by `CONFIG`.

## Observed result

All records were syntactically valid and finite. The initial projection
remained accurate (L2 about `3.013e-15`, sampled Linf about `1.049e-13`),
but the final error depended pathologically on the time-step count:

| scheme | steps | `dt` | final L2 error |
| --- | ---: | ---: | ---: |
| DG | 4 | 0.025 | `5.2683403681e-4` |
| DG | 8 | 0.0125 | `1.2228819666e8` |
| PR | 4 | 0.025 | `9.2763596092e-1` |
| PR | 8 | 0.0125 | `4.3434098042e-1` |
| PR | 16 | 0.00625 | `1.0856729877e4` |
| BE | 4 | 0.025 | `7.1503905089e12` |
| BE | 8 | 0.0125 | `3.5171153024e-4` |

The DG `N=8` command was repeated and returned the same
`1.2228819666e8` L2 error, so this was not a transient process-launch
failure.

As a control, changing only the mesh to `3x3x3` produced the expected
decrease from `N=4` to `N=8`:

| scheme | `N=4` L2 | `N=8` L2 |
| --- | ---: | ---: |
| DG | about `8.25e-5` | about `2.04e-5` |
| PR | about `2.18e-4` | about `9.37e-5` |
| BE | about `4.17e-4` | about `1.44e-4` |

This isolates the failure from the exactly representable initial projection
and makes an even-mesh issue in the shared transient/core path the current
working diagnosis. Root-cause localization and a production fix require a
separate core-change review and commit.

## Consequence for profiles

The Stage-2 `smoke` and `smoke-refined` profiles use the verified
`3x3x3` control mesh and retain strict initial-state and refinement gates.
The declared `temporal-full` profile intentionally remains the specified
`4x4x4` matrix, so its scientific execution is blocked pending the separate
core fix. A successful dry-run of its 792-case plan is not evidence that those
numerical runs are valid.
