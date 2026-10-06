# iGRM and DPG problem drivers

[Running guide](running.md) · [Core architecture](architecture.md) · [Documentation index](README.md) · [Repository README](../README.md)

All paths and commands are relative to the repository root.

## Pure Diffusion iGRM

Arguments:

```text
<size> <order> <procx> <procy> <procz> <steps> <dt> [scheme]
```

The optional `scheme` argument selects the iGRM time scheme:

```text
dg    Douglas-Gunn, default
pr    stabilized first-order 3D PR selector
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

## iGRM L2 Projection

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
pr    stabilized first-order 3D PR selector
be    Backward Euler
```

Example:

```bash
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/igrm_l2 2 2 2 3 3 3 1 1 1 1 1 1 pr
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/igrm_l2 2 2 2 3 3 3 1 1 1 1 1 1 1.0 be
```

## iGRM Heat

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

## iGRM Eriksson (MUMPS)

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

## 3D DG-iGRM Stokes

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

## 3D DPG/iGRM Pollution

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
