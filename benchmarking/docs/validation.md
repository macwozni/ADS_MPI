# Full-field MPI/OpenMP validation

[Benchmarking index](../README.md) · [Repository README](../../README.md)

The oracle is implemented in the shared `ads_benchmark.validation`
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
make benchmark-validate RUN_ID=validation-full-local \
  BENCHMARK_PLAN_ARGS='--available-mpi-slots 8 --available-cpu-slots 32'
```

`make -C benchmarking plan BENCHMARK_PROFILE=validation-full` remains the
portable full-matrix dry-run when that allocation is unavailable.

For a declared 24-CPU/6-rank allocation, a representative Y/Z/uneven subset
can be selected without changing the profile:

```bash
make benchmark-validate RUN_ID=validation-subset-local \
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
