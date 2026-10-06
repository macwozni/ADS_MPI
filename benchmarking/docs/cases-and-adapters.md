# Manufactured cases and adapters

[Benchmarking index](../README.md) · [Repository README](../../README.md)

## Manufactured transients

The registered exact case is `temporal-polynomial` on the unit cube:

```text
q(s)       = s^2 (3 - 2s)
Q(x,y,z)   = q(x) q(y) q(z)
u(x,y,z,t) = exp(-t) Q(x,y,z)
```

Because `q'(0)=q'(1)=0`, the normal flux is zero on every face. The
production weak form is `M u_t + K u = f`, with `K` representing
`-Delta`. Since `q''(s)=6-12s`, the source for
`u_t - Delta(u) = f` is

```text
f(x,y,z,t) = exp(-t) [
    -Q(x,y,z)
    + (12x-6) q(y)q(z)
    + (12y-6) q(x)q(z)
    + (12z-6) q(x)q(y)
]
```

The initial field is `Q`. It is exactly representable when every trial
degree is at least three. The harness performs one mass-only projection at
`t=0`, verifies that state, and then performs exactly `N` physical updates
with `dt=T/N`. It asserts both the computed final time and the ADS state time
equal `T`.

Source evaluation follows the production scheme tables:

| Scheme | source time in physical step `[t_n,t_n+dt]` |
| --- | --- |
| Douglas-Gunn (`dg`) | `t_n + dt/2` |
| split Backward Euler (`be`) | `t_n + dt` |
| cyclic Peaceman-Rachford (`pr`) | `t_n + dt/6`, `t_n + dt/2`, `t_n + 5dt/6` |

The callback time used inside OpenMP RHS assembly is thread-private.

The independent `spatial-cosine` case is used by the `h`, `p`, `validation`,
and `strong` families:

```text
R(x,y,z)   = cos(pi*x) cos(pi*y) cos(pi*z)
u(x,y,z,t) = exp(-t) R(x,y,z)
f(x,y,z,t) = (3*pi^2 - 1) exp(-t) R(x,y,z)
```

Here `Delta(R)=-3*pi^2*R`, so the source has the displayed sign for the
production convention `u_t-Delta(u)=f`. The normal derivative vanishes on all
six faces. Unlike `temporal-polynomial`, this field is not exactly
representable by any finite spline degree. Its initial projection error is
therefore a measured spatial error, not a failed temporal-case initialization
gate.

## Problem adapters

All adapters use the shared initialization, projection, DG/PR/BE wrappers,
measurement, MPI status handling, and cleanup. They are still distinct
benchmark entry points and select the appropriate production path:

| Adapter | benchmark-local responsibility |
| --- | --- |
| `igrm_l2` | Uses the generic production `ComputePointForRHS` path. The public problem performs one solve; this adapter deliberately reuses the initialized ADS state for `N` physical updates to `T`, so the run is a real temporal experiment. |
| `igrm_heat` | Uses the problem's production `heat_igrm_rhs_point` full-RHS callback and preserves its physical-time branch, while replacing the scalar source with the manufactured source. No VTI is written in timed steps. |
| `pure_diffusion_igrm` | Uses the generic production RHS path and accepts independent test/trial degrees in all three axes instead of the public driver's isotropic degree shortcut. |

No production source or public problem CLI is changed. The executables are:

```text
benchmarking/build/<debug|release>/EXEC/igrm_l2_manufactured
benchmarking/build/<debug|release>/EXEC/igrm_heat_manufactured
benchmarking/build/<debug|release>/EXEC/pure_diffusion_igrm_manufactured
```

Their common direct CLI is:

```text
<scheme> <T> <steps> <nx> <ny> <nz>
<ptest-x> <ptest-y> <ptest-z>
<ptrial-x> <ptrial-y> <ptrial-z>
<proc-x> <proc-y> <proc-z>
<sample-points> <write-samples:0|1> [exact-case]
```

Normally use the Python runner, which derives this command from the normalized
case and checks the returned record against it.
