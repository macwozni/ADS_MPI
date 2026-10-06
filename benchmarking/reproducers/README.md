# Benchmark reproducers

[Benchmarking index](../README.md) · [Convergence guide](../docs/convergence.md)

These files preserve the command lines and observations that motivated core or
benchmark changes. They are evidence records, not a statement that every
number still reproduces at the current `HEAD`.

| Reproducer | Status |
|---|---|
| [Spatial/degree solver defects](spatial-degree-solver-defects.md) | Fixed by core commit `5e1161b`; pre-fix commands and post-fix evidence retained. |
| [Temporal convergence instability](temporal-convergence-instability.md) | Historical pre-fix measurements; several listed failures were removed by `5e1161b`, but post-fix spot checks still contain a degree-dependent non-monotone series. Full profiles require a fresh qualification run. |
| [Original even-mesh transient](even-mesh-transient.md) | Superseded by the wider temporal investigation and retained only as the first observation. |

When a reproducer changes status, update this index first. Keep old measured
tables with their commit/date context instead of rewriting them as current
results.
