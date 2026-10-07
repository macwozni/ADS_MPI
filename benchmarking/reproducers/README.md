# Benchmark reproducers

[Benchmarking index](../README.md) · [Convergence guide](../docs/convergence.md)

These files preserve the command lines and observations that motivated core or
benchmark changes. They are evidence records, not a statement that every
number still reproduces at the current `HEAD`.

| Reproducer | Status |
|---|---|
| [Spatial/degree solver defects](spatial-degree-solver-defects.md) | Fixed by core commit `5e1161b`; pre-fix commands and post-fix evidence retained. |
| [Temporal convergence instability](temporal-convergence-instability.md) | Fixed in two core steps: `5e1161b` repaired the mixed system and the later stabilized PR correction removed the residual non-monotone series. The post-fix 72-case validation passed; the 792-case run completed but remains red for separate `(9,8)` projection-oracle failures and timeouts. |
| [Original even-mesh transient](even-mesh-transient.md) | Superseded by the wider temporal investigation and retained only as the first observation. |

When a reproducer changes status, update this index first. Keep old measured
tables with their commit/date context instead of rewriting them as current
results.
