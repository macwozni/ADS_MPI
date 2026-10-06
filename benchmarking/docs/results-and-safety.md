# Measurements, numerical gates, results, and safety

[Benchmarking index](../README.md) · [Repository README](../../README.md)

## Measurement and machine-readable results

Legacy convergence and validation cases retain the exact single-invocation
result schema (`warmups=0`, `samples=1`). Strong and weak cases use the same
executor, adapter, launcher, and result store but add the repeated timing arrays,
reliability threshold/flag, and effective OpenMP environment under `timing`.
Any failed, timed-out, unparsable, or domain-invalid repetition fails the
whole case; resume accepts only a complete, strictly revalidated aggregate.
The local and cluster profiles use a 0.05-second minimum solver-side sample;
the small executable smoke profile uses 0.001 seconds and still reports any
shorter sample as unreliable.

The benchmark harness computes the L2 error and solution norm with an
independent Gauss rule of `min(10,p_trial+3)` points per axis. It deliberately
does not reuse the trial-space assembly rule: doing so can alias the
non-polynomial cosine field and report a false zero error at low degree. Ten
points is the library's supported limit and exactly integrates a squared
degree-nine spline. For Linf and field fingerprinting, rank zero reconstructs
the global coefficients and evaluates an endpoint-inclusive regular grid of
`sample_points_per_axis^3` points, with X varying fastest. Compensated sums
produce a deterministic `field_checksum`; Linf and checksum are then
broadcast. Sampling and optional CSV output occur outside the measured
physical-step interval.

With `write_samples=1`, the case directory also receives
`field_samples.csv` with columns:

```text
x,y,z,numerical,exact,error
```

Spatial analysis requires this artifact for both members of every `dt/dt2`
pair. It is opened relative to the verified case directory without following
symlinks, must be a regular UTF-8 file, and is bounded to 64 MiB before it is
parsed. Profile validation therefore caps `points_per_axis` at 76 whenever
`write_samples=true`; the fixed-width CSV for `77^3` samples cannot fit under
that storage contract. Unwritten sampling retains the executable's general
limit of 257 points per axis. Run loading retains only a lazy safe reference,
and parsed errors use a
compact binary64 array; a full profile therefore never keeps all CSV texts or
per-sample Python objects in memory. The artifact is not copied into
`result.json` or the final report.

Each successful executable writes exactly one stdout line beginning with
`ADS_BENCHMARK_RESULT `, followed by strict JSON. Its domain fields are:

```text
schema_version, kind, exact_case, problem, scheme,
requested_final_time, actual_final_time, time_step, steps,
initial_l2_error, initial_linf_error,
l2_error, linf_error, solution_l2_norm, field_checksum,
sample_points_per_axis, field_samples_written,
physical_step_wall_seconds, solver_status
```

The parser rejects duplicate/missing/extra keys, non-finite values, a nonzero
solver status, multiple tagged records, and any mismatch with the planned
case. The stored outer result adds the stable `case_id`, complete normalized
`configuration`, process `timing.wall_seconds`, and the parsed
`domain_result`.

## Numerical gates

Before timing, every executable requires finite initial metrics. The exactly
representable `temporal-polynomial` additionally requires
`initial_l2_error <= 1e-10`, `initial_linf_error <= 1e-10`, and agreement of
the projected solution norm with `(13/35)^(3/2)` within `1e-10`. The
non-polynomial spatial case retains its finite initial projection errors
instead of applying that exact-representation gate. Failures, MPI/MUMPS
errors, NaN, or Inf produce a nonzero process result.

The smoke oracle loads exactly the nine passed `N=4` results and their nine
`N=8` counterparts, rechecks the initial errors, and requires a strict L2
decrease for every problem/scheme pair. This is a refinement sanity check, not
yet a formal order estimate.

## Result layout and safety

A written or executed plan owns:

```text
benchmarks/<run-id>/manifest.json
benchmarks/<run-id>/execution.json                    # executed runs
benchmarks/<run-id>/cases/<case-id>/status.json
benchmarks/<run-id>/cases/<case-id>/result.json
benchmarks/<run-id>/cases/<case-id>/stdout.log
benchmarks/<run-id>/cases/<case-id>/stderr.log
benchmarks/<run-id>/cases/<case-id>/field_samples.csv  # only when requested
benchmarks/<run-id>/analysis/analysis.json
benchmarks/<run-id>/analysis/analysis.csv
benchmarks/<run-id>/analysis/convergence.png           # only with --plot
benchmarks/<run-id>/analysis/strong-scaling.png         # strong --plot
benchmarks/<run-id>/analysis/weak-scaling.png           # weak --plot
benchmarks/<ab-id>-baseline/analysis/ab-schedule.json    # ref A/B workflow
benchmarks/<ab-id>-candidate/analysis/ab-schedule.json   # ref A/B workflow
benchmarks/<ab-id>-candidate/analysis/comparison.json    # ref A/B report
benchmarks/<ab-id>-candidate/analysis/comparison.csv     # ref A/B report
benchmarks/<parent-run-id>/shards/shard-000000.json     # sharded plans
```

The manifest contains the schema version, full expanded configuration and case
IDs, configuration hash, Git commit, dirty-tree flag, and a SHA-256 content
fingerprint of tracked and nonignored untracked worktree state. Ignored
benchmark/build output is excluded. Manifest, execution, result, and comparison
documents carry a `schema_version` and a `kind`; status documents carry their
own `schema_version` plus the exact case identity and state. Strict input readers
for manifests, execution records, statuses, and results reject unknown, missing,
or inconsistent identity fields. Directory-mode A/B outputs are
written only when explicit `--output-json`/`--output-csv` paths are passed
through `BENCHMARK_COMPARE_ARGS`.

Writes are atomic and exclusive.
Descriptor-relative operations reject traversal, symlink swaps,
sibling-prefix tricks, and adoption of unrelated directories, including the
legacy prototype.

The entire `benchmarks/` tree is ignored runtime data, not versioned source.
Do not stage manifests, fields, logs, reports, or comparison runs. Neither
`make clean` nor `make clean-benchmark-build` deletes a run; the latter removes
only marker-owned framework caches and isolated builds. There is deliberately
no broad results-clean target and no `rm -rf benchmarks` workflow. Remove an
individual run manually only after identifying that exact run directory and
deciding its data is no longer needed.
