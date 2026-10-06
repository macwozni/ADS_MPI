# Temporal, spatial, and degree convergence

[Benchmarking index](../README.md) · [Repository README](../../README.md)

## Temporal convergence analysis

Analyze a complete frozen run and optionally request a log-log plot:

```bash
make benchmark-analyze RUN_ID=temporal-validation-local
make benchmark-analyze RUN_ID=temporal-validation-local \
  BENCHMARK_ANALYZE_ARGS=--plot
```

The analyzer first requires every result named by the frozen manifest and
rejects missing, duplicated, misnamed, or configuration-mismatched records.
It then groups cases that differ only in time resolution. An unfiltered full
run requires all eight levels `N=4,8,...,512`; an explicitly filtered run may
select a common, complete sequence of at least four levels. L2 and Linf are
analyzed independently using errors against the registered analytical
solution, not differences between two numerical solutions. A full report
contains all seven local estimates

```text
log(e_i/e_(i+1)) / log(dt_i/dt_(i+1))
```

and a centered log-log regression over the four finest usable consecutive
points. The acceptance contract is deliberately fixed:

| scheme | theoretical order | accepted regression order |
| --- | ---: | ---: |
| DG | 2 | `[1.70, 2.30]` |
| PR | 1 | `[0.80, 1.20]` |
| BE | 1 | `[0.80, 1.20]` |

Both metrics also require regression `R^2 >= 0.98` and at least a fourfold
reduction of analytical error. Non-finite data, a nonzero solver status,
wrong final time, missing levels, duplicate IDs, or nonmonotonic analytical
error before a credible plateau are hard failures.

A positive local order above four times the scheme's theoretical order is
also a hard metric failure. This deliberately retains catastrophic coarse
solves followed by an apparently clean tail, such as the measured PR defect;
the tail regression alone cannot certify the full series.

A plateau may remove only a finest-level suffix from regression. Its local
orders must stay below `max(0.20, 0.25*p_expected)` and its errors within a
`max/min <= 1.25` band. At least two consecutive plateau transitions are
required, except that one final transition is allowed at an actual roundoff
scale (`<= 1e-10` times the solution norm). The excluded points and their
local orders remain in the report. Isolated spikes and large coarse-grid
errors are therefore never relabeled as plateau.

Analysis writes `analysis/analysis.json` and `analysis/analysis.csv` below the
run directory. `--plot` additionally requests `analysis/convergence.png`; if
matplotlib is unavailable, only the plot is omitted and numerical analysis
still runs. Reanalysis removes any older PNG before deciding whether the new
report has a plot, so a stale image can never describe newer JSON/CSV output.

## Spatial and degree convergence analysis

Run `make benchmark-analyze RUN_ID=<name>` for either an `h` or a `p` run.
The analyzer first enforces the same frozen-manifest and result-integrity
checks as temporal analysis. For each spatial point and each metric it then
subtracts the two sampled error fields. Because both runs have the same exact
case and final time, this is exactly the sampled numerical-field difference:

```text
delta_i = (u_dt - u_exact)_i - (u_dt/2 - u_exact)_i
        = u_dt,i - u_dt/2,i

D_L2   = sqrt(composite_trapezoid_3d(delta_i^2))
D_Linf = max_i abs(delta_i)
c_t    = D_metric / E_metric(dt/2)
```

The CSV header, row count, X-fastest regular-grid order, coordinates, exact
values, pointwise errors, reported sampled Linf, and deterministic field
checksum are all rebound and checked before subtraction. Merely subtracting
the two scalar error norms is
not used: that would be only a lower bound and could falsely pass two very
different fields with equal norms. A point is spatially reliable only when
`c_t <= 0.10`. The accepted value is `E(dt/2)`. A larger or non-finite
indicator marks the metric
`time_dominated/unreliable`, requests a finer `dt` pair, and suppresses every
order or degree-drop estimate involving that point. This same-layout temporal
check is local to spatial convergence; cross-layout MPI/OpenMP field equivalence remains a
separate workflow.

`D_L2` and `D_Linf` are discrete common-grid diagnostics. The reported L2
error in the denominator is the independent Gauss-integrated error, while the
reported Linf error and `D_Linf` are maxima on the declared regular grid, not
a proof of the continuous supremum. In particular, `33^3` is a smoke/degree
resolution and `65^3` gives only two sampling intervals per element at the
finest `32^3` h level. Reports preserve `sample_points_per_axis` and the
difference method so a denser study can be distinguished rather than silently
compared as the same Linf experiment.

The spatial workflow exposed a structurally singular mixed
iGRM system at the coarser PR step and at required low degree. The solver fix
was kept in the separate core commit `5e1161b`; pre-fix commands, root cause,
and post-fix numerical evidence are retained in
[`reproducers/spatial-degree-solver-defects.md`](../reproducers/spatial-degree-solver-defects.md).

For `h`, only cases that differ in the element counts are compared. The report
retains all qualified L2 and Linf values, every local
`log(E_i/E_(i+1))/log(h_i/h_(i+1))` estimate, and the regression slope versus
`h`. For isotropic `p`, the two enrichment families are kept separate and the
report records error decrease and adjacent ratios versus the actual trial
degree; it does not invent a fixed algebraic convergence order. A plateau is
accepted only after at least two strict, qualified error decreases, so a flat
sequence from the first degree cannot pass as convergence. Anisotropic degree
vectors remain vectors in configuration, identity, grouping, JSON, and CSV
reports. Their rotation spread is an audit result, not an ordered convergence
gate: the ADI axis order does not justify an arbitrary rotational-equality
tolerance. Solver/roundoff plateaus and time-dominated points or ranges are
retained as auditable data rather than silently discarded.
