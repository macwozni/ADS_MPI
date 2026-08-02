# ADS MPI IGA

This repository contains a Fortran/MPI implementation of the alternating
direction solver (ADS) for isogeometric analysis (IGA), including the current
iGRM-oriented workflow and several example problems.

The code is research software. The actively maintained build path in the
current tree is the `mymake` build directory.

## Repository Layout

```text
src/        Core ADS, IGA basis, MPI communication, sparse assembly, output
problems/   Problem drivers and problem-specific data/callbacks
mymake/     Active make-based build configuration
tests/      Unit, integration, MPI, and driver-CLI regression suites
doxygen/    Generated/documentation support files
```

## Dependencies

The default `mymake/m_options` expects:

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
- Bash and GNU `timeout` (for CLI and MPI error-path tests)

Edit `mymake/m_options` for local library paths and compiler flags.

### Library Versions In The Current Local Configuration

The active `mymake/m_options` currently points to:

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
problems/<problem>/main.F90         Driver
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
current `mymake/m_files` problem groups use `RHS_fun.F90`. The shared
quadrature-level RHS assembly lives in `src/RHS_eq.F90`.

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

Build the default L2 projection driver from the repository root:

```bash
make -C mymake clean
make -C mymake
```

The default build uses:

```make
SOURCE_ALL = $(SOURCES) $(L2)
EXEC = l2
```

so the default executable is:

```text
mymake/EXEC/l2
```

When building another problem, pass both `SOURCE_ALL` and `EXEC` to `make`.
Use the same variables for `make clean`, because the object list depends on the
selected source group.

Examples:

```bash
make -C mymake clean SOURCE_ALL='$(SOURCES) $(L2)' EXEC=l2
make -C mymake       SOURCE_ALL='$(SOURCES) $(L2)' EXEC=l2

make -C mymake clean SOURCE_ALL='$(SOURCES) $(HEAT)' EXEC=heat
make -C mymake       SOURCE_ALL='$(SOURCES) $(HEAT)' EXEC=heat

make -C mymake clean SOURCE_ALL='$(SOURCES) $(ERIKSSON)' EXEC=eriksson
make -C mymake       SOURCE_ALL='$(SOURCES) $(ERIKSSON)' EXEC=eriksson

make -C mymake clean SOURCE_ALL='$(SOURCES) $(PURE_DIFFUSION_IGRM)' EXEC=pure_diffusion_igrm
make -C mymake       SOURCE_ALL='$(SOURCES) $(PURE_DIFFUSION_IGRM)' EXEC=pure_diffusion_igrm

make -C mymake clean SOURCE_ALL='$(SOURCES) $(OIL)' EXEC=oil
make -C mymake       SOURCE_ALL='$(SOURCES) $(OIL)' EXEC=oil

make -C mymake clean SOURCE_ALL='$(SOURCES) $(IGRM_L2)' EXEC=igrm_l2
make -C mymake       SOURCE_ALL='$(SOURCES) $(IGRM_L2)' EXEC=igrm_l2

make -C mymake clean SOURCE_ALL='$(SOURCES) $(IGRM_HEAT)' EXEC=igrm_heat
make -C mymake       SOURCE_ALL='$(SOURCES) $(IGRM_HEAT)' EXEC=igrm_heat
```

The CMake files are still present, but this README documents the current
make-based workflow.

## Running

Run executables with MPI. The number of MPI ranks must match:

```text
procx * procy * procz
```

For one-rank smoke tests:

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
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/heat 2 1 0 1 1 1 1
```

### Eriksson

Arguments:

```text
<size> <order> <steps> <dt> <procx> <procy> <procz>
```

Example:

```bash
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/eriksson 2 1 0 1 1 1 1
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
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/pure_diffusion_igrm 2 1 1 1 1 0 1
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

## iGRM Mesh Assumptions

The mixed iGRM matrix path assumes:

- the test and trial spaces use the same geometric mesh,
- repeated knots are allowed and do not change the geometric mesh,
- the test degree is greater than the trial degree.

In other words, the distinct knot locations must match, while knot
multiplicities may differ.

## Testing

Run tests from the repository root. The pFUnit suites use
`PFUNIT_ROOT=/opt/lib/pfunit/PFUNIT-4.16` by default, and the MPI suites use
`MPIEXEC=/opt/lib/mpich-5.0.0/bin/mpiexec`. Override either variable in the
shell block below when using another installation.

Every maintained suite has its own `tests/<suite>/GNUmakefile`; the legacy
top-level `tests/GNUmakefile` is not a complete test runner. The following
procedure cleans and runs all 18 library/unit/MPI suites, then forcibly
rebuilds all seven problem drivers and runs the driver-CLI suite:

```bash
set -eu

MPIEXEC=${MPIEXEC:-/opt/lib/mpich-5.0.0/bin/mpiexec}
PFUNIT_ROOT=${PFUNIT_ROOT:-/opt/lib/pfunit/PFUNIT-4.16}

for suite in \
  ads_directional_solve \
  ads_error_propagation \
  ads_lifecycle \
  argument_parser \
  basis \
  communicators \
  igrm_space \
  int2str \
  mumps_solver \
  my_mpi \
  norm_l2 \
  operator_assembly \
  parallelism \
  reorder_rhs \
  rhs_assembly \
  rhs_eq \
  solution_reconstruction \
  time_scheme
do
  make -j1 -C "tests/${suite}" clean
  make -j1 -C "tests/${suite}" \
    MPIEXEC="${MPIEXEC}" PFUNIT_ROOT="${PFUNIT_ROOT}" run
done

make -j1 -C tests/driver_cli -B MPIEXEC="${MPIEXEC}" run
```

`-B` is intentional for `driver_cli`: that suite's `clean` target is a no-op,
and a forced build verifies every current problem group, including
`IGRM_HEAT`. Several MPI suites exercise up to eight ranks and OpenMP thread
counts 1, 2, and 4.

Finish the full regression by rebuilding the corrected default target and
running a positive one-rank smoke test for every driver. The preceding
`driver_cli` command has already built the six non-default executables.

```bash
make -j1 -C mymake clean
make -j1 -C mymake
"${MPIEXEC:-/opt/lib/mpich-5.0.0/bin/mpiexec}" -n 1 \
  ./mymake/EXEC/l2 2 2 2 1 1 1 1

"${MPIEXEC:-/opt/lib/mpich-5.0.0/bin/mpiexec}" -n 1 \
  ./mymake/EXEC/heat 2 1 0 1 1 1 1
"${MPIEXEC:-/opt/lib/mpich-5.0.0/bin/mpiexec}" -n 1 \
  ./mymake/EXEC/eriksson 2 1 0 1 1 1 1
"${MPIEXEC:-/opt/lib/mpich-5.0.0/bin/mpiexec}" -n 1 \
  ./mymake/EXEC/pure_diffusion_igrm 2 1 1 1 1 0 1 dg
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
  1 1 1 1 0.1 dg
```

## Notes

- A warning about `/opt/lib/parmetis/lib/include` may appear if that include
  directory does not exist locally. The builds used during recent smoke tests
  still completed with this warning.
- `mymake/EXEC/` and `mymake/_OBJ/` contain generated build outputs.
- Doxygen-style comments are used throughout `src` and the migrated problem
  drivers.
