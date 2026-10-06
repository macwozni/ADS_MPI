# Profiles and normalized cases

[Benchmarking index](../README.md) · [Repository README](../../README.md)

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
  [`reproducers/spatial-degree-solver-defects.md`](../reproducers/spatial-degree-solver-defects.md);
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
[`reproducers/temporal-convergence-instability.md`](../reproducers/temporal-convergence-instability.md)
and the original
[`reproducers/even-mesh-transient.md`](../reproducers/even-mesh-transient.md).
The two smoke profiles retain their valid two-level refinement check; the
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
