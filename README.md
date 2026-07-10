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
tests/      Local/unit-test experiments
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

Build from `mymake`:

```bash
cd mymake
make
```

The default build uses:

```make
SOURCE_ALL = $(SOURCES) $(IGRM_L2)
EXEC = l2
```

so the default executable is:

```text
mymake/EXEC/l2
```

To clean the default build:

```bash
cd mymake
make clean
```

When building another problem, pass both `SOURCE_ALL` and `EXEC` to `make`.
Use the same variables for `make clean`, because the object list depends on the
selected source group.

Examples:

```bash
cd mymake

make clean SOURCE_ALL='$(SOURCES) $(L2)' EXEC=l2_projection
make       SOURCE_ALL='$(SOURCES) $(L2)' EXEC=l2_projection

make clean SOURCE_ALL='$(SOURCES) $(HEAT)' EXEC=heat
make       SOURCE_ALL='$(SOURCES) $(HEAT)' EXEC=heat

make clean SOURCE_ALL='$(SOURCES) $(ERIKSSON)' EXEC=eriksson
make       SOURCE_ALL='$(SOURCES) $(ERIKSSON)' EXEC=eriksson

make clean SOURCE_ALL='$(SOURCES) $(PURE_DIFFUSION_IGRM)' EXEC=pure_diffusion_igrm
make       SOURCE_ALL='$(SOURCES) $(PURE_DIFFUSION_IGRM)' EXEC=pure_diffusion_igrm

make clean SOURCE_ALL='$(SOURCES) $(OIL)' EXEC=oil
make       SOURCE_ALL='$(SOURCES) $(OIL)' EXEC=oil

make clean SOURCE_ALL='$(SOURCES) $(IGRM_L2)' EXEC=igrm_l2
make       SOURCE_ALL='$(SOURCES) $(IGRM_L2)' EXEC=igrm_l2
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
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/l2_projection 2 2 2 1 1 1 1
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
<isizex> <isizey> <isizez> <order> <procx> <procy> <procz>
<isizex> <isizey> <isizez> <order> <procx> <procy> <procz> <scheme>
<isizex> <isizey> <isizez> <order> <procx> <procy> <procz> <tau> <scheme>
```

The current driver is a fixed small iGRM demonstration: it parses these
arguments but uses a hard-coded `2 x 2 x 2` mesh with `p_test = 3` and
`p_trial = 1`.

The optional `scheme` argument selects the iGRM time scheme. Forward Euler is
not accepted here because this driver exercises the `MultiStep` iGRM schemes:

```text
dg    Douglas-Gunn, default
pr    Peaceman-Rachford
be    Backward Euler
```

Example:

```bash
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/igrm_l2 2 2 2 1 1 1 1 pr
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/igrm_l2 2 2 2 1 1 1 1 1.0 be
```

## iGRM Mesh Assumptions

The mixed iGRM matrix path assumes:

- the test and trial spaces use the same geometric mesh,
- repeated knots are allowed and do not change the geometric mesh,
- the test degree is greater than the trial degree.

In other words, the distinct knot locations must match, while knot
multiplicities may differ.

## Notes

- A warning about `/opt/lib/parmetis/lib/include` may appear if that include
  directory does not exist locally. The builds used during recent smoke tests
  still completed with this warning.
- `mymake/EXEC/` and `mymake/_OBJ/` contain generated build outputs.
- Doxygen-style comments are used throughout `src` and the migrated problem
  drivers.
