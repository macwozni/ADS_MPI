# Classic ADS problem drivers

[Running guide](running.md) · [Documentation index](README.md) · [Repository README](../README.md)

All paths and commands are relative to the repository root.

## L2 Projection

Arguments:

```text
<isizex> <isizey> <isizez> <order> <procx> <procy> <procz>
```

Example:

```bash
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/l2 2 2 2 1 1 1 1
```

## Heat

Arguments:

```text
<size> <order> <steps> <dt> <procx> <procy> <procz>
```

Example:

```bash
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/heat 2 1 1 0.01 1 1 1
```

### `marcinlos/iga-ads` `heat_3d` compatibility

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

## Eriksson

Arguments:

```text
<size> <order> <steps> <dt> <procx> <procy> <procz>
```

Example:

```bash
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/eriksson 2 1 1 0.01 1 1 1
```

## Oil

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
