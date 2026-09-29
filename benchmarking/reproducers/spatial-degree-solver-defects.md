# Spatial/degree solver defects exposed and resolved during Stage 4

Stage 4 initially exposed deterministic MUMPS factorization failures. They
were solver defects, not failed convergence tolerances, so no benchmark oracle
or requested degree range was relaxed. The core repair is isolated in commit
`5e1161b9d8d0411c92b95bdcea7ade121d8a3e5c` (`fix: restore nonsingular mixed
iGRM saddle system`). This document retains the pre-fix reproducers and the
post-fix evidence.

Commands were run with the release benchmark adapters, one MPI rank, and:

```bash
ROOT=/home/maciekw/Dropbox/shared/Maciek_Paszynski/CPC_MPI_2017/code/mpi-ads
MPIEXEC=/opt/lib/mpich-5.0.0/bin/mpiexec
```

## Pre-fix failures

The coarser PR member of the h-smoke temporal pair failed for all three
adapters. A representative command was:

```bash
"$MPIEXEC" -n 1 \
  "$ROOT/benchmarking/build/release/EXEC/igrm_l2_manufactured" \
  pr 0.01 64 2 2 2 3 3 3 2 2 2 1 1 1 33 1 spatial-cosine
```

MUMPS reported `INFO(1:2)=-10,8` during factorization. The same case at
`N=128` happened to pass. The generated diagnostic run
`stage4-all9-spatial-smoke-20260929` therefore completed only 15 of 18 cases.

The lowest required degree in the p workflow also failed at both members of
its temporal pair:

```bash
"$MPIEXEC" -n 1 \
  "$ROOT/benchmarking/build/release/EXEC/igrm_l2_manufactured" \
  dg 0.01 128 2 2 2 2 2 2 1 1 1 1 1 1 17 0 spatial-cosine
```

Both `N=128` and `N=256` reported `INFO(1:2)=-10,6`. Higher-degree variants
could either fail or produce a catastrophic field, depending on the numerical
pivot selected for the same structurally singular matrix.

## Root cause and repair

The mixed one-dimensional assembler used a trial-first saddle system,

```text
[ A_trial  C   ]
[ C^T      0_test ]
```

When the test space is enriched, `dim(test) > dim(trial)`. The zero test block
therefore contains an unavoidable nullspace that the coupling cannot remove.
Whether MUMPS rejected the matrix or returned a large solution depended on
roundoff and pivoting.

The repaired canonical residual-minimization system is test-first:

```text
[ G_test  B ] [r] = [Ft]
[ B^T     0 ] [u]   [ F]
```

`G_test` is the test-space Gram block, the zero block belongs to the smaller
trial space, and the directional right-hand side is packed as `[Ft;F]`.
For nonsymmetric first-derivative forms, transposing the coupling also swaps
the two derivative mix coefficients. Regression tests now check the exact
block layout, the coupling transpose, right-hand-side packing, and full rank
for mesh two at test/trial degrees `2/1` and `5/3`.

## Post-fix validation

Fresh release builds of all three manufactured adapters passed. For PR on
mesh `2^3`, test/trial degree `3/2`, `T=0.01`, and `N=64/128`, every adapter
returned process status zero and solver status zero; corresponding sampled
fields were byte-identical across adapters. For `igrm_l2` the fine-step errors
were:

| metric | `N=128` error | sampled temporal fraction |
| --- | ---: | ---: |
| L2 | `1.1584975657e-2` | `2.6146e-4` |
| Linf | `3.3286399353e-2` | `3.7745e-4` |

DG degree probes at `N=128/256` also passed:

| trial/test family | fine-step L2 | fine-step Linf | max temporal fraction |
| --- | ---: | ---: | ---: |
| `p=1`, enrichment `+1` | `7.47244e-2` | `6.93207e-1` | `2.17e-4` |
| `p=2`, enrichment `+1` | `1.16374e-2` | `3.52903e-2` | `1.01e-4` |
| `p=3`, enrichment `+1` | `3.26069e-3` | `1.28930e-2` | `4.58e-3` |
| `p=3`, enrichment `+2` | `3.26159e-3` | `1.28979e-2` | `5.55e-3` |

The required low-degree p case now passes, its error decreases strongly over
the first three degrees, and every sampled temporal fraction is below the
fixed `0.10` Stage-4 limit. A BE `p=2` probe at `N=64/128` passed as well.
Generated reproducer runs remain ignored and are not part of either commit.
