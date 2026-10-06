# OpenMP performance regression

[Test suite index](README.md) · [Repository README](../README.md)

The real wall-clock regression is part of both `make test-driver` and the
complete `make test`; it can also be run by itself with:

```bash
make test-performance
```

It builds only `igrm_l2` with `BUILD=release` in the isolated
`build/openmp-performance` tree, then runs the same `32x32x32`, test-degree-three,
trial-degree-two DG case with one and four threads on one MPI rank. One warm-up
pair is followed by three measured pairs in alternating OMP1/OMP4 order.
`OMP_DYNAMIC=FALSE`, `OMP_PROC_BIND=close`, and `OMP_PLACES=cores` are fixed by
default. Every run, including warm-ups, must produce the same finite `31^3`
VTK field as one common reference to `1e-12`.
The gate fails unless the median paired speedup is at least `1.10` and OMP4 is
faster in a strict majority of pairs. Results are written to
`build/openmp-performance/openmp-performance.json`.

Run this wall-clock gate on an otherwise idle worker with at least four
exclusive CPU cores. MPI launcher flags must leave the single rank access to
all four cores; binding that rank to one core makes the OMP4 measurement invalid
and is expected to fail the speedup gate.

On a pinned, otherwise idle runner, the result from a known-good revision can
also gate absolute regressions of both OMP1 and OMP4 medians:

```bash
make test-performance \
  PERFORMANCE_BASELINE=/absolute/path/to/known-good.json \
  PERFORMANCE_MAX_REGRESSION=1.15
```

`PERFORMANCE_WARMUPS`, `PERFORMANCE_SAMPLES`, `PERFORMANCE_TIMEOUT`, and
`PERFORMANCE_MIN_SPEEDUP` are configurable. `PERFORMANCE_TIMEOUT` applies to
one process launch; the complete timed suite has the independent
`PERFORMANCE_SUITE_TIMEOUT=3600s`. If the number of warm-ups or samples, or the
per-launch timeout, is raised substantially, raise this outer timeout too. A
baseline is accepted only when its schema, workload, thread counts, MPI launch
command/options, binding, and placement match the current run, and is meaningful
only on the same pinned hardware and software stack.
`make test-performance-self-test` exercises configuration, VTK parsing,
full-field comparison, report generation, thresholds, and baseline validation
without launching MPI.

Relative `PERFORMANCE_BUILD_ROOT` values are resolved at the repository root.
The release tree carries an ownership marker; builds refuse a nonempty unowned
directory, and cleanup resolves symlinks before rejecting source, test, normal
build, or other unsafe in-repository destinations. Use `make
clean-performance` to remove only the owned release artifacts and JSON report;
unrelated files in that tree are retained.

Normal oil runs remain stochastic when `ADS_OIL_RANDOM_SEED` is unset. The
equivalent smoke commands are listed below for manual diagnostics. `make
problems` builds all required executables first.

The DG, PR, and BE coefficient tables satisfy the independent first-order
operator/source balance checked by the time-scheme unit tests. A scalar
transient oracle additionally executes the RHS state selectors and requires
at least first-order convergence for PR and BE. The current MPI/OpenMP runs
preserve a nonzero discrete equilibrium at roundoff for all three schemes. The
manufactured pure-diffusion oracle checks every VTK sample against both the
analytic equilibrium and the initial field, then compares the complete DG,
PR, and BE fields pairwise at every step. Together these checks guard against
inconsistent directional, state-selection, or forcing weights.

```bash
make problems
"${MPIEXEC:-/opt/lib/mpich-5.0.0/bin/mpiexec}" -n 1 \
  ./mymake/EXEC/l2 2 2 2 1 1 1 1

"${MPIEXEC:-/opt/lib/mpich-5.0.0/bin/mpiexec}" -n 1 \
  ./mymake/EXEC/heat 2 1 1 0.01 1 1 1
"${MPIEXEC:-/opt/lib/mpich-5.0.0/bin/mpiexec}" -n 1 \
  ./mymake/EXEC/eriksson 2 1 1 0.01 1 1 1
"${MPIEXEC:-/opt/lib/mpich-5.0.0/bin/mpiexec}" -n 1 \
  ./mymake/EXEC/pure_diffusion_igrm 3 1 1 1 1 1 0.1 dg
"${MPIEXEC:-/opt/lib/mpich-5.0.0/bin/mpiexec}" -n 1 \
  ./mymake/EXEC/oil \
  2 1 1 1 1 1 0.1 \
  1 0.5 0.5 0.5 \
  1 0.25 0.25 0.25
"${MPIEXEC:-/opt/lib/mpich-5.0.0/bin/mpiexec}" -n 1 \
  ./mymake/EXEC/igrm_l2 \
  2 2 2 3 3 3 1 1 1 \
  1 1 1 pr
"${MPIEXEC:-/opt/lib/mpich-5.0.0/bin/mpiexec}" -n 1 \
  ./mymake/EXEC/igrm_heat \
  2 2 2 3 3 3 2 2 2 \
  1 1 1 1 0.001 dg
"${MPIEXEC:-/opt/lib/mpich-5.0.0/bin/mpiexec}" -n 1 \
  ./mymake/EXEC/igrm_eirksson \
  4 4 4 3 3 3 2 2 2 \
  1 1 1
"${MPIEXEC:-/opt/lib/mpich-5.0.0/bin/mpiexec}" -n 1 \
  ./mymake/EXEC/igrm_stokes \
  2 2 2 2 2 2 2 2 2 \
  1 1 1
ADS_POLLUTION_OUTPUT_RESOLUTION=4 \
  "${MPIEXEC:-/opt/lib/mpich-5.0.0/bin/mpiexec}" -n 1 \
  ./mymake/EXEC/igrm_pollution \
  4 0 1 0 2 1 1 1 1 1
```
