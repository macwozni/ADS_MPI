# A/B numerical and performance comparison

[Benchmarking index](../README.md) · [Repository README](../../README.md)

Compare two already completed run directories through the public root target:

```bash
make benchmark-compare \
  BENCHMARK_BASELINE_RESULTS=benchmarks/baseline-run \
  BENCHMARK_CANDIDATE_RESULTS=benchmarks/candidate-run
```

The loader applies the normal manifest, status, command, log, result, and
adapter verification to both runs. It then requires identical configuration
hashes, case-ID sets, and per-case normalized configurations, plus compatible
machine, toolchain, external-library, launcher, topology, binding, and OpenMP
execution provenance. Repository/build-root paths are normalized so isolated
worktrees can be compared; the source refs themselves may differ.

The decision is deliberately numerics-first. Every candidate field must pass
the complete-field comparison against its baseline counterpart before any
timing is classified. One numerical mismatch blocks timing for the whole
comparison. For every eligible case the report retains each side's sample
count, minimum, maximum, median, MAD, relative MAD, and range, then computes

```text
median_ratio = candidate_median / baseline_median
```

The default regression threshold is 5%, so a regression requires
`median_ratio > 1.05`. The default minimum is five measured samples per side;
the policy rejects values below three. Default complete-field tolerances are
absolute `1e-11` and relative `1e-10`. Too few samples, a timing already marked
unreliable, malformed timing, or missing legacy execution provenance produces
`inconclusive`, never a green result or a regression. Legacy runs can therefore
still be inspected numerically without pretending their timings are
reproducible. Incompatible configuration or present-but-different execution
provenance is reported as `incompatible`. Only
`no-regression-detected` returns success; regression, numerical mismatch,
incompatibility, and inconclusive evidence return a nonzero verdict.

The default command prints a summary. Optional versioned JSON and flat CSV,
as well as a different threshold, use the same public Make target:

```bash
make benchmark-compare \
  BENCHMARK_BASELINE_RESULTS=benchmarks/baseline-run \
  BENCHMARK_CANDIDATE_RESULTS=benchmarks/candidate-run \
  BENCHMARK_COMPARE_ARGS='--regression-threshold 0.03 --minimum-samples 5 --output-json /tmp/ads-ab.json --output-csv /tmp/ads-ab.csv'
```

To build and run two refs on the same host without touching a dirty main
worktree, supply refs instead of result directories:

```bash
make benchmark-compare \
  RUN_ID=ab-local \
  BENCHMARK_BASELINE_REF=HEAD~1 \
  BENCHMARK_CANDIDATE_REF=HEAD \
  BENCHMARK_COMPARE_PROFILE=strong-scaling-smoke \
  BENCHMARK_PLAN_ARGS='--available-mpi-slots 2 --available-cpu-slots 2' \
  BENCHMARK_COMPARE_ARGS='--max-retries 1'
```

The ref orchestrator resolves both refs once to full commit SHAs, creates two
distinct detached marker-owned worktrees and separate isolated build roots,
and never checks out the caller's tree. Matched cases use a deterministic
case-level `AB, BA, AB, ...` order to reduce monotonic machine drift. This is
not sample-level alternation: each side of a case still performs all of its own
warmups followed by all of its measured samples. The two exclusive result runs
are `benchmarks/ab-local-baseline/` and
`benchmarks/ab-local-candidate/` for the example above. Both receive the
versioned schedule; the candidate run receives the JSON/CSV comparison report.

After both sides execute and the report is published, the orchestrator removes
only clean, marker-verified worktrees through `git worktree remove`, even when
the comparison verdict is regression or inconclusive. An execution error or
other exception before report publication, and any modified worktree, is
retained for diagnosis rather than force-deleted. Set
`BENCHMARK_WORKSPACE_PARENT` to choose the temporary parent; otherwise the
system temporary directory is used. The parent must be outside the benchmark
repository so a retained diagnostic workspace cannot silently become an
untracked change in the caller's tree.

For a retained failure, first inspect the printed workspace and both worktrees.
When they are clean and no longer needed, remove the two registered paths with
`git worktree remove <workspace>/baseline` and
`git worktree remove <workspace>/candidate`; then inspect the ownership marker,
unlink that one marker, and remove the now-empty workspace with `rmdir`. Never
use `git worktree remove --force`, `git worktree prune`, or a recursive delete
as a substitute for resolving modified/unowned contents.

The comparator reports MAD and range but deliberately applies a transparent
median-ratio threshold rather than a statistical significance test. Very short
or noisy measurements can therefore cross a tight threshold even for identical
builds; choose a duration and sample count appropriate to the machine and treat
the recorded dispersion as part of the verdict review.

Result-directory inputs and ref inputs are mutually exclusive, and both sides
of the selected mode are required. Use normal profile filters in
`BENCHMARK_PLAN_ARGS` for a small local slice; a full cluster profile remains
a manual allocation-aware workflow.
