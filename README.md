# ADS MPI IGA

This repository contains a Fortran/MPI implementation of the alternating
direction solver (ADS) for isogeometric analysis (IGA), including the current
iGRM-oriented workflow and several example problems.

The code is research software. The repository-root `Makefile` is the public
interface and only orchestrates builds owned by `src/`, `problems/`, and
`tests/`.

## Repository Layout

```text
Makefile              Public orchestration, configuration, and documentation
m_options             Active local compiler, dependency, and tool configuration
makeconfig/           Selectable GNU and Intel example configurations
src/GNUmakefile       Core-library build
src/sources.mk        Ordered core-source manifest
problems/GNUmakefile  Problem-build aggregator
problems/problems.mk  Public problem-name/directory manifest
problems/*/GNUmakefile
                      Local source list, arguments, build, run, and cleanup
tests/GNUmakefile     Test-group aggregator
tests/{src,problems,driver}/GNUmakefile
                      Test-suite manifests for each group
tests/*/GNUmakefile   Build, run, and cleanup for one concrete suite
mymake/               Compatibility entry point and generated artifacts
doxygen/              Generated/documentation support files
```

## Dependencies

The default repository-root `m_options` expects:

- an MPI Fortran compiler
- MUMPS
- BLAS
- LAPACK
- ScaLAPACK
- METIS
- ParMETIS
- GKlib

The complete test procedure additionally expects:

- pFUnit 4 (for the pFUnit-based test suites)
- Python 3.10 or newer (for layout checks and numerical VTI integration checks)
- Bash and GNU `timeout` (for CLI and MPI error-path tests)

Edit `m_options` for local library paths and compiler flags.

### Library Versions In The Current Local Configuration

The active `m_options` currently points to:

```text
MPI compiler: /opt/lib/mpich-5.0.0/bin/mpif90
MUMPS:        /opt/lib/MUMPS_5.8.2/
BLAS:         /opt/lib/lapack-3.12.1/lib64/libblas.a
LAPACK:       /opt/lib/lapack-3.12.1/lib64/liblapack.a
ScaLAPACK:    /opt/lib/scalapack-2.3.2/lib64/libscalapack.a
METIS:        /opt/lib/metis/lib/libmetis.a
ParMETIS:     /opt/lib/parmetis/lib/libparmetis.a
GKlib:        /opt/lib/GKlib/lib64/libGKlib.a
```

The code is still an MPI program. Directional linear solves are executed with
MUMPS on `MPI_COMM_SELF`, so each rank solves its local gathered systems while
the ADS workflow itself remains distributed through MPI.

The current debug-style options use bounds checking and AddressSanitizer:

```make
-O0 -g -fcheck=all -fbounds-check -fsanitize=address
```

## Current Problem API

Each migrated problem follows this structure:

```text
problems/<problem>/input_data.F90   Command-line parameters and problem data
problems/<problem>/RHS_fun.F90      Pointwise forcing callback
problems/<problem>/*_solver.F90     Optional problem-local solver
problems/<problem>/main.F90         Driver
problems/<problem>/GNUmakefile      Local build/run contract
```

The active callback has the signature:

```fortran
function forcing(un, du, X) result(ret)
   real(kind = 8), intent(in) :: un
   real(kind = 8), intent(in), dimension(3) :: du
   real(kind = 8), intent(in), dimension(3) :: X
   real(kind = 8) :: ret
end function forcing
```

Problem-local legacy `RHS_eq.F90` files may still exist as references, but the
active source lists in the problem-local GNUmakefiles use `RHS_fun.F90`. The
shared quadrature-level RHS assembly lives in `src/RHS_eq.F90`.

## Time Integration Schemes

There are two ADS time-advancement paths:

- `Step` is the single-step Forward Euler path. Its `alpha_step = 1` is
  intentional and should not be treated as a Douglas-Gunn or ADI scheme.
- `MultiStep` is the three-substep directional path used by iGRM time schemes.

The implemented time-scheme wrappers are:

```text
ForwardEuler3DStep          3D Forward Euler wrapper
DouglasGunn3DStep           3D Douglas-Gunn wrapper
PeacemanRachford3DStep      3D cyclic Peaceman-Rachford wrapper
BackwardEuler3DStep         3D split Backward Euler wrapper
```

`ForwardEuler3DStep` delegates to `Step`. The iGRM/ADI path is configured
once before the time loop and then reused through a persistent `TimeScheme3D`
object. Use the persistent configurators in `src/time_scheme.F90`:

- `ConfigureDouglasGunn3DTimeScheme`
- `ConfigurePeacemanRachford3DTimeScheme`
- `ConfigureBackwardEuler3DTimeScheme`
- `ConfigureMassOnly3DTimeScheme`

The `*TimeScheme` configurators fill `TimeScheme3D` with the RHS coefficient
table, the RHS derivative-state selector, and the directional LHS mixing table.
Time loops should call the named DG/PR/BE wrappers with that existing object,
not rebuild scheme coefficients every step. The lower level
`ConfigureDouglasGunn3D`, `ConfigurePeacemanRachford3D`, and
`ConfigureBackwardEuler3D` routines remain available for callers that need raw
coefficient tables.

The active iGRM direction keeps the residual-minimization gram block mass-only,
while the scheme operator is applied through the coupling blocks. This avoids
singular mixed iGRM stiffness blocks.

iGRM space compatibility is checked once during problem setup with
`ValidateIGRMTimeSchemeSpaces`; it is not repeated inside each time-step
wrapper.

The iGRM/ADI configurators are diffusion-oriented by default. They accept an
optional `include_transport` flag for first-derivative transport terms; the
`pure_diffusion_igrm` driver configures its scheme with
`include_transport=.false.`.

## Building

The repository-root `Makefile` is the public build interface. Run `make help`
to list every target and `make show-config` to display the effective tool and
library paths.

The tracked repository-root `m_options` is the single active configuration
file for the compiler, build profile, MPI launcher, numerical libraries,
pFUnit, Python, Doxygen, and test timeouts. The old `mymake/m_options` path is
retained as a compatibility wrapper which includes that same root file, so
there are no duplicate settings to keep synchronized. Command-line
assignments still override values for one invocation:

```bash
make show-config
make BUILD=release
make MPIEXEC=/path/to/mpiexec MPIFC=/path/to/mpif90 test
```

Ready-to-select examples live in `makeconfig/`:

```bash
make list-configs
make CONFIG=makeconfig/gnu-debug.mk show-config
make CONFIG=makeconfig/gnu-release.mk all
make CONFIG=makeconfig/intel-oneapi-debug.mk all
make CONFIG=makeconfig/intel-oneapi-release.mk all
```

The GNU examples use `mpif90`/`mpiexec` from `PATH`. The Intel examples use
the oneAPI `mpiifx` wrapper and Intel-specific module, checking, and OpenMP
flags. They cover library and problem builds. The individual test-suite
Makefiles still contain GNU-only diagnostic flags, so a complete Intel test
run additionally requires parameterizing those flags and using an Intel-built
pFUnit. Numerical libraries must likewise match the selected compiler; their
local paths remain centralized in root `m_options`.

The default `BUILD=debug` preserves the existing bounds checks and
AddressSanitizer flags. `BUILD=release` selects `RELEASE_OPTS` from the same
configuration file.

The make hierarchy is executable at every level:

```text
Makefile
|-- src/GNUmakefile -- src/sources.mk
|-- problems/GNUmakefile -- problems/problems.mk
|   `-- problems/<problem>/GNUmakefile
`-- tests/GNUmakefile
    |-- tests/src/GNUmakefile ------ tests/<source-suite>/GNUmakefile
    |-- tests/problems/GNUmakefile - tests/<problem-suite>/GNUmakefile
    `-- tests/driver/GNUmakefile --- tests/driver_cli/GNUmakefile
```

The root does not own core or problem source lists, per-problem defaults, or
test-suite registration. Consequently the corresponding subtree can also be
used directly, for example:

```bash
make -C src library
make -C problems build PROBLEM=heat
make -C problems/heat show-run
make -C tests/src run-suite TEST_SUITE=rhs_assembly
```

### Build targets

The default command builds the static ADS library and all ten problem
drivers. Problem builds are deliberately serialized because problem-local
Fortran modules reuse the names `input_data`, `RHS_fun`, and `main`.

```bash
# Static library plus every problem.
make
make all
make build-all

# Individual layers.
make library
make problems

# One selected problem: equivalent forms.
make build PROBLEM=heat
make build-heat
make heat

make list-problems
```

Generated files are placed in:

```text
mymake/LIB/libads.a
mymake/EXEC/l2
mymake/EXEC/heat
mymake/EXEC/eriksson
mymake/EXEC/pure_diffusion_igrm
mymake/EXEC/oil
mymake/EXEC/igrm_l2
mymake/EXEC/igrm_heat
mymake/EXEC/igrm_eirksson
mymake/EXEC/igrm_stokes
mymake/EXEC/igrm_pollution
```

The lower-level `mymake/makefile` remains usable as a compatibility entry
point. It delegates `library` to `src/` and problem targets to `problems/`; it
no longer owns a duplicate source list or compilation implementation.
Historical named `EXEC=<problem>` builds remain supported. Explicit
`SOURCE_ALL` builds use a deliberately isolated compatibility adapter: it
compiles only the caller-supplied ordered list, puts its objects and modules in
`<BUILD_ROOT>/_LEGACY_OBJ/<EXEC>` (by default below `mymake/`), and never
reuses the official core or problem object directories. Thus legacy custom
programs remain buildable without moving source ownership back into `mymake`.

Both the core and each problem keep a stamp of their effective compiler and
flags. A debug/release or GNU/Intel switch therefore rebuilds incompatible
objects, while an unchanged second invocation remains incremental.

### Cleanup targets

```bash
make clean-build
make clean-problems
make clean-library
make clean-legacy
make clean-tests
make clean-docs
make clean          # all generated build, test, and documentation artifacts
make distclean      # same as clean; m_options and makeconfig/ are preserved
```

The CMake files are still present, but this README documents the current
make-based workflow.

## Running

Every problem has a root-level `run-<problem>` target. The generic form is:

```bash
make run PROBLEM=heat
make run-heat
make run-heat ARGS='4 2 3 0.01 2 1 1' NP=2 OMP_NUM_THREADS=4
```

Use `make run-help` to print the complete argument syntax and defaults for all
ten problems. `make show-run PROBLEM=heat` prints the effective executable,
arguments, MPI/OpenMP settings, environment, and output directory without
building or launching it.

If `ARGS` is omitted, each target uses a small one-rank example with
`steps=1` for transient problems. Output is written to `output/<problem>` by
default; set `RUN_DIR` to choose another directory.
The historical per-problem form, for example `heat_ARGS='...'`, is also
accepted when the common `ARGS` variable is empty.

All runtime controls can be supplied on the make command line:

```text
ARGS                 raw problem command-line arguments
NP                   MPI rank count, default 1
MPIEXEC              MPI launcher from m_options
MPIEXEC_FLAGS        additional launcher flags
MPI_NP_FLAG          rank-count flag, default -n
OMP_NUM_THREADS      OpenMP thread count, default 1
OMP_DYNAMIC          OpenMP dynamic teams, default FALSE
OMP_PROC_BIND        OpenMP binding policy, default close
RUN_DIR              output working directory, default output/<problem>
RUN_ENV              additional environment assignments
OIL_SEED             shortcut for ADS_OIL_RANDOM_SEED on run-oil
```

For example:

```bash
make run-oil OIL_SEED=20260811 OMP_NUM_THREADS=4
make run-igrm_heat ARGS='2 2 2 3 3 3 2 2 2 2 1 1 1 0.001 be' NP=2
make run-igrm_eirksson NP=1
make run-igrm_stokes NP=1
make run-igrm_pollution NP=1
```

The number of MPI ranks must match:

```text
procx * procy * procz
```

The executables may also be invoked directly:

```bash
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/l2 2 2 2 1 1 1 1
```

### L2 Projection

Arguments:

```text
<isizex> <isizey> <isizez> <order> <procx> <procy> <procz>
```

Example:

```bash
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/l2 2 2 2 1 1 1 1
```

### Heat

Arguments:

```text
<size> <order> <steps> <dt> <procx> <procy> <procz>
```

Example:

```bash
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/heat 2 1 1 0.01 1 1 1
```

#### `marcinlos/iga-ads` `heat_3d` compatibility

The heat driver mirrors the problem definition from
[`examples/heat/heat_3d`](https://github.com/marcinlos/iga-ads/blob/fa6e64b50dba44709039bdb0c37971de1bda9af3/examples/heat/heat_3d.hpp)
at revision `fa6e64b50dba44709039bdb0c37971de1bda9af3`. In particular, both use
the compactly supported initial state

```text
r2 = min(8 * ((x - 0.5)^2 + (y - 0.5)^2 + (z - 0.5)^2), 1)
u0 = (r2 - 1)^2 * (r2 + 1)^2
```

followed by an L2 projection and the zero-source Forward Euler update
`M u(n+1) = M u(n) - dt K u(n)`, with natural homogeneous Neumann boundary
conditions. The matching invocation for the upstream defaults is:

```bash
OMP_NUM_THREADS=1 /opt/lib/mpich-5.0.0/bin/mpiexec -n 1 \
  ./mymake/EXEC/heat 12 2 5000 1e-7 1 1 1
```

Here `step0.vti` is the projected initial state and `step1.vti` through
`step5000.vti` are exactly 5000 physical updates. The automated compatibility
test uses the same mesh, degree, and time step but stops after one update. It
compares all `31^3` sampled values at VTI precision against a numerical oracle
produced by a probe reproducing the C++ `heat_3d` assembly with the unmodified
upstream basis, quadrature, and LAPACK implementation at the pinned revision.
Results agree numerically, rather than bit for bit, because this implementation
uses sparse MUMPS solves and a different summation order. Unlike the silent
upstream example, the full local command writes 5001 VTI files; use `steps=1`
for a quick compatibility check.

### Eriksson

Arguments:

```text
<size> <order> <steps> <dt> <procx> <procy> <procz>
```

Example:

```bash
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/eriksson 2 1 1 0.01 1 1 1
```

### Pure Diffusion iGRM

Arguments:

```text
<size> <order> <procx> <procy> <procz> <steps> <dt> [scheme]
```

The optional `scheme` argument selects the iGRM time scheme:

```text
dg    Douglas-Gunn, default
pr    Peaceman-Rachford
be    Backward Euler
```

Example:

```bash
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/pure_diffusion_igrm 3 1 1 1 1 1 0.1 dg
```

Scheme-selection examples:

```bash
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/pure_diffusion_igrm 2 1 1 1 1 1 0.1 dg
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/pure_diffusion_igrm 2 1 1 1 1 1 0.1 pr
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/pure_diffusion_igrm 2 1 1 1 1 1 0.1 be
```

### Oil

Arguments:

```text
<size> <order> <procx> <procy> <procz> <steps> <dt> \
<npumps> <pump_x> <pump_y> <pump_z> ... \
<ndrains> <drain_x> <drain_y> <drain_z> ...
```

Example with one pump and one drain:

```bash
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/oil \
  2 1 1 1 1 1 0.1 \
  1 0.5 0.5 0.5 \
  1 0.25 0.25 0.25
```

### iGRM L2 Projection

Arguments expected by the parser:

```text
<nelem_x> <nelem_y> <nelem_z> <ptest_x> <ptest_y> <ptest_z> <ptrial_x> <ptrial_y> <ptrial_z> <procx> <procy> <procz>
<nelem_x> <nelem_y> <nelem_z> <ptest_x> <ptest_y> <ptest_z> <ptrial_x> <ptrial_y> <ptrial_z> <procx> <procy> <procz> <scheme>
<nelem_x> <nelem_y> <nelem_z> <ptest_x> <ptest_y> <ptest_z> <ptrial_x> <ptrial_y> <ptrial_z> <procx> <procy> <procz> <tau> <scheme>
```

The three `ptest_*` values configure the enriched iGRM test-space degrees.
The three `ptrial_*` values configure the trial-space degrees.

The optional `scheme` argument selects the iGRM time scheme. Forward Euler is
not accepted here because this driver exercises the `MultiStep` iGRM schemes:

```text
dg    Douglas-Gunn, default
pr    Peaceman-Rachford
be    Backward Euler
```

Example:

```bash
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/igrm_l2 2 2 2 3 3 3 1 1 1 1 1 1 pr
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/igrm_l2 2 2 2 3 3 3 1 1 1 1 1 1 1.0 be
```

### iGRM Heat

Arguments:

```text
<nelem_x> <nelem_y> <nelem_z> \
<ptest_x> <ptest_y> <ptest_z> \
<ptrial_x> <ptrial_y> <ptrial_z> \
<procx> <procy> <procz> <steps> <dt> [scheme]
```

The test-space degree must be greater than the trial-space degree in every
direction. `steps` is the number of physical heat steps after the initial
mass-only projection. The optional scheme is `dg` (the default), `pr`, or
`be`; all three physical schemes are configured without transport terms.

Example with one initial projection and one physical time step:

```bash
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/igrm_heat \
  2 2 2 3 3 3 2 2 2 \
  1 1 1 1 0.1 dg
```

### iGRM Eriksson (MUMPS)

The executable is named `igrm_eirksson` to match the repository target name.
It is a stationary 3D tensor-product extension of the upstream 2D Eriksson
iGRM problem on the unit cube:

```text
-0.01 Laplace(u) + (1, 1, 1) dot grad(u) = f,
u = 0 on the boundary.
```

The manufactured solution is `u(x,y,z) = g(x)g(y)g(z)`, where
`g(x) = x - (exp(-(1-x)/0.01) - exp(-100))/(1-exp(-100))`. This retains the
Eriksson outflow layer while evaluating its exponential without overflow.
The driver assembles the complete saddle-point system

```text
[ G   -B ] [ r ] = [ -f ]
[ B^T  0 ] [ u ]   [  0 ]
```

and performs one MUMPS factorization and solve. This is a direct iGRM-MUMPS
driver; it has no time-step or iterative defect-correction loop.
The test and trial spaces share one uniform mesh and differ by polynomial
enrichment. Upstream subdivision/adaptation modes are not exposed by this
baseline driver.

The saddle formulation and one-shot direct solve follow
[`examples/erikkson/erikkson_mumps.hpp`](https://github.com/marcinlos/iga-ads/blob/959ac6e03b6b0e7332c1e253c1907bc502269d74/examples/erikkson/erikkson_mumps.hpp)
at upstream revision `959ac6e03b6b0e7332c1e253c1907bc502269d74`; this repository extends
that two-dimensional reference problem tensorially to three dimensions.

Arguments:

```text
<nelem_x> <nelem_y> <nelem_z> \
<ptest_x> <ptest_y> <ptest_z> \
<ptrial_x> <ptrial_y> <ptrial_z> \
<procx> <procy> <procz>
```

Trial-space degrees must be positive, and each test-space degree must be
strictly greater than the corresponding trial-space degree.

The reference implementation assembles and solves the complete system on MPI
rank zero with MUMPS on `MPI_COMM_SELF`, then broadcasts the solution to the
other ranks. It is intended as a correctness baseline and is not scalable to
large distributed 3D systems. A scalable version requires distributed matrix
assembly and a collective MUMPS solve.

Example:

```bash
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/igrm_eirksson \
  4 4 4 3 3 3 2 2 2 1 1 1
```

### 3D DG-iGRM Stokes

`igrm_stokes` implements the manufactured, stationary three-dimensional
Stokes problem from the upstream `DGiGRM_stokes_3D` example:

```text
-Laplace(v) + grad(p) = f,
div(v) = 0
```

The manufactured velocity boundary values are imposed weakly with Nitsche
terms. The residual/test fields are discontinuous between elements, while the
trial velocity and pressure use conforming tensor-product B-spline spaces. The
problem-local assembler includes the volume velocity/pressure coupling and the
DG facet jump, average, and penalty terms. It forms the complete mixed iGRM
saddle system, augments it with a Lagrange multiplier enforcing zero mean
trial pressure, and solves it once with MUMPS.

The formulation follows
[`DGiGRM_stokes_3D`](https://github.com/marcinlos/iga-ads/blob/959ac6e03b6b0e7332c1e253c1907bc502269d74/examples/dg/laplace.cpp#L1552)
at upstream revision `959ac6e03b6b0e7332c1e253c1907bc502269d74`.

Arguments:

```text
<nelem_x> <nelem_y> <nelem_z> \
<ptest_x> <ptest_y> <ptest_z> \
<ptrial_x> <ptrial_y> <ptrial_z> \
<procx> <procy> <procz>
```

All polynomial degrees must be positive. Each test-space degree must be
greater than or equal to its trial-space counterpart. The test space is
discontinuous and the trial space uses the upstream `C1` continuity (`C0` for
linear splines, where `C1` is impossible). The default of four elements and
equal degree two in every direction exactly selects the upstream spaces.

The root rank assembles and solves this correctness baseline, then broadcasts
the four trial fields. `result.vti` contains a three-component `Velocity`
array and a scalar `Pressure` array. The driver also reports algebraic RMS and
relative residuals, velocity and pressure L2 errors, and the L2 divergence.

```bash
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/igrm_stokes \
  4 4 4 2 2 2 2 2 2 1 1 1
```

### 3D DPG/iGRM Pollution

`igrm_pollution` implements the transient three-dimensional
advection-diffusion source problem from the upstream `pollution_dpg_3d`
example on the cube `[0, 5000]^3`:

```text
du/dt - div(D grad(u)) + wind(t) dot grad(u) = emission,
u(0) = 0,
D = (50, 50, 0.5).
```

The compact source is centered at `(3000, 2000, 2000)` with radius `25`.
For `r2 = min(sum(((x - center) / 25)^2), 1)`, its value is
`(r2 - 1)^2 (r2 + 1)^2`. The time step is the upstream value `1.8`; the
three directional substeps use `dt/3`. The time-dependent upstream wind law
is evaluated consistently at the beginning of every physical step, including
`t=0`, removing the prototype's discontinuous first update.

The trial and enriched test spaces have independently selectable degrees and
conforming continuities (`0 <= C <= p-1`). The upstream parser also permits
`C=-1`, but that broken space is deliberately rejected here: the prototype
does not contain the facet flux and trace terms required for a DG diffusion
formulation. `adapt=1` enables the upstream nonuniform knot map in the x
direction; `adapt=0` uses a uniform mesh. The implementation uses only the
existing ADS interfaces and keeps all pollution-specific assembly and model
data inside `problems/igrm_pollution`.

As in the upstream prototype, diffusion uses the natural zero-flux boundary
condition. Advection remains in the strong volume form and has no separately
imposed inflow value. These conditions are part of the comparison model.

The formulation follows
[`pollution_dpg_3d.hpp`](https://github.com/marcinlos/iga-ads/blob/959ac6e03b6b0e7332c1e253c1907bc502269d74/examples/pollution/pollution_dpg_3d.hpp)
and
[`dpg3d.cpp`](https://github.com/marcinlos/iga-ads/blob/959ac6e03b6b0e7332c1e253c1907bc502269d74/examples/pollution/dpg3d.cpp)
at upstream revision `959ac6e03b6b0e7332c1e253c1907bc502269d74`.
Obvious prototype copy-and-paste defects in the z-direction indexing and
temporary-buffer selection are corrected rather than reproduced.

Arguments:

```text
<N> <adapt:0|1> <p_trial> <C_trial> <p_test> <C_test> \
<steps> <procx> <procy> <procz>
```

`N` is the number of elements in every direction. `steps` is the number of
physical updates and must be at least one for an actual simulation. The
process-grid product must equal the number of MPI ranks.

The initial field is written to `out_0.vti`, followed by `out_1.vti` through
`out_<steps>.vti`. The default output resolution is 100 intervals per
direction, matching upstream. Set `ADS_POLLUTION_OUTPUT_RESOLUTION` to a
positive integer for smaller diagnostic files; automated tests use `4`.

The current problem-local direct implementation computes the dense
directional solves on MPI rank zero and broadcasts the final coefficients.
Multiple ranks therefore verify deterministic replication but do not yet
accelerate this problem. The reported `maximum coefficient abs` diagnostic is
the largest trial coefficient magnitude, not a separately sampled field
maximum.

```bash
ADS_POLLUTION_OUTPUT_RESOLUTION=4 \
  /opt/lib/mpich-5.0.0/bin/mpiexec -n 1 \
  ./mymake/EXEC/igrm_pollution 4 0 1 0 2 1 1 1 1 1
```

## iGRM Mesh Assumptions

The shared mixed iGRM matrix path used by the scalar problems assumes:

- the test and trial spaces use the same geometric mesh,
- repeated knots are allowed and do not change the geometric mesh,
- the test degree is greater than the trial degree.

In other words, the distinct knot locations must match, while knot
multiplicities may differ. `igrm_stokes` uses its own DG facet assembler and
therefore permits equal test and trial degrees while retaining different
continuities. `igrm_pollution` likewise owns its DPG assembly locally and
accepts independently configured conforming trial/test continuities.

## Testing

The root Makefile delegates to `tests/GNUmakefile`, which delegates to the
`src`, `problems`, and `driver` group GNUmakefiles; those in turn delegate to
the individual suites. The active configuration is forwarded through every
level. The defaults below live in root `m_options` and may be overridden on
the command line:

```text
PFUNIT_ROOT=/opt/lib/pfunit/PFUNIT-4.16
MPIEXEC=/opt/lib/mpich-5.0.0/bin/mpiexec
MPIFC=/opt/lib/mpich-5.0.0/bin/mpif90
MUMPS_DIR=/opt/lib/MUMPS_5.8.2
SUITE_TIMEOUT=600s
DRIVER_CLI_TIMEOUT=20s
DRIVER_SMOKE_TIMEOUT=60s
DRIVER_INTEGRATION_TIMEOUT=90
SKIP_MPI_CASES=0
```

### One source file, one primary test file

Every active library source in the `SRC_FILES` manifest in `src/sources.mk` has
exactly one primary, authored test file. The tab-separated mapping is stored
in `tests/test-map.tsv`. A suite may still use fixtures, stubs, generated
pFUnit sources, or link other production modules; those support files are not
additional primary tests for the mapped source.

Validate this invariant before changing or running the suites:

```bash
make test-layout
```

`check-layout` asks `src/GNUmakefile` for its active source list and rejects
missing mappings, inactive sources, duplicate sources or test files, and
paths that do not exist. It also keeps the `tests/src` suite manifest
synchronized with the map, validates the three group runners, and verifies
that every library suite references its production source and primary test.
All registered library, problem, and driver suites must have `all`, `run`, and
`clean` targets; unregistered suite directories are rejected. Problem modules
are kept in separate suites because several drivers deliberately use the same
Fortran module names (`input_data` and `RHS_fun`).

### Test targets

The runner separates library tests, problem-specific callback tests, and
full-driver tests:

```bash
# Check only the one-to-one layout.
make test-layout

# Run the 28 suites mapped to src/*.F90.
make test-src

# Run all problem-specific input, RHS, and solver suites.
make test-problems

# Build all ten problem executables and run CLI, smoke, and numerical
# integration tests.
make test-driver

# Individual driver layers.
make test-cli
make test-smoke
make test-integration

# Run the complete regression above in one command.
make test
make check

# Clean-build every test without executing it, list suites, or clean-run one.
make test-build
make test-list
make test-suite TEST_SUITE=rhs_assembly
```

Every test level can be invoked directly. Driver executables built in
`mymake/EXEC` are deliberately retained by `clean-tests`; use `clean-build` or
`clean` to remove them.

```bash
make -j1 -C tests run-src
make -j1 -C tests/problems run-suite TEST_SUITE=heat_rhs_fun
make -j1 -C tests/rhs_assembly run
```

The aggregate `run-src`, `run-problems`, and `run-driver` targets clean each
suite before running it, so compiler or flag changes cannot silently reuse a
stale test executable. Driver executables are rebuilt unconditionally.

Pass non-default tool locations once at the root; they are forwarded to every
suite:

```bash
make test \
  PFUNIT_ROOT=/path/to/pfunit \
  MPIEXEC=/path/to/mpiexec \
  MPIFC=/path/to/mpif90 \
  MUMPS_DIR=/path/to/mumps
```

The runner is deliberately serialized to keep diagnostics deterministic and
avoid oversubscribing MPI/OpenMP test jobs. Each problem build nevertheless
has a private `_OBJ` directory, so identically named problem modules cannot be
reused accidentally. Each suite is also protected by `SUITE_TIMEOUT`. MPI
suites exercise up to eight ranks, and relevant OpenMP tests compare thread
counts 1, 2, 4, and 8. The top-level runner requires a POSIX environment with
Bash and the coreutils `timeout` command; selected error-path probes
additionally use POSIX process primitives.

### Positive smoke and numerical integration tests

`test` and `test-driver` execute positive one-rank smoke tests for every real
driver and the numerical integration matrix automatically. The integration
matrix does more than check process status:

- transient drivers use `steps=1`, so iteration 1 executes with `t>0`;
- L2 runs with degree-four support and a `2x2x2` MPI grid, validates local
  coefficients on every rank through a global reduction, and treats `not OK`
  as a test failure even though the driver itself currently exits successfully;
- heat and Eriksson VTI files must be valid `Float64` XML with exactly
  `31^3` finite values, change after initialization, and agree between serial
  and hybrid MPI/OpenMP runs;
- every sampled value of the standard 12-element, degree-two heat case must
  match, at VTI precision, the projected initial state and first time step
  generated by the pinned `iga-ads` `heat_3d` implementation;
- DG, PR, and BE run through the real MUMPS path for every iGRM driver, with
  a global finite zero-solution oracle for the pure-diffusion example;
- the iGRM-heat DG case must be dissipative for the tested time step, and all
  three schemes' serial and hybrid VTI results must agree;
- iGRM Eriksson requires finite, small algebraic residuals, a nonzero interior,
  six homogeneous faces, decreasing L2 error after refinement, and identical
  serial and hybrid VTI output;
- iGRM Stokes requires finite, small algebraic residuals, finite velocity,
  pressure, and divergence errors, improving refined errors, valid coupled
  velocity/pressure VTI output, and identical serial and hybrid results;
- iGRM pollution executes a real physical source step, validates the zero
  initial state and nonzero finite concentration in `out_1.vti`, and compares
  the serial and hybrid MPI/OpenMP fields on a reduced `5^3` output grid;
- oil uses the opt-in `ADS_OIL_RANDOM_SEED` test seed and compares a positive
  drained result across one/four OpenMP threads and a hybrid MPI run.

Normal oil runs remain stochastic when `ADS_OIL_RANDOM_SEED` is unset. The
equivalent smoke commands are listed below for manual diagnostics. `make
problems` builds all required executables first.

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

## Documentation

Documentation uses the tracked `doxygen.conf`. Doxygen must run from the
repository root because its input and output paths are relative.

```bash
make docs-check   # validate/expand the Doxygen configuration
make docs-html    # doxygen/html/index.html
make docs-pdf     # HTML plus doxygen/latex/refman.pdf
make docs         # same complete HTML+PDF documentation build
make clean-docs
```

`DOXYGEN` can be overridden in `m_options` or on the command line. PDF
generation additionally requires a LaTeX installation with `pdflatex` and
`makeindex`.

## Notes

- A warning about `/opt/lib/parmetis/lib/include` may appear if that include
  directory does not exist locally. The builds used during recent smoke tests
  still completed with this warning.
- `mymake/EXEC/`, `mymake/LIB/`, and `mymake/_OBJ/` contain the public
  generated executables, library, and core objects; each problem keeps its
  own generated `_OBJ/` directory. Explicit legacy `SOURCE_ALL` builds use
  marker-protected `<BUILD_ROOT>/_LEGACY_OBJ/<EXEC>` directories.
- Doxygen-style comments are used throughout `src` and the migrated problem
  drivers.
