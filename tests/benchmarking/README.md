# Benchmark framework tests

[Repository test guide](../README.md) · [Benchmarking index](../../benchmarking/README.md)

```bash
make benchmark-self-test
```

All test sources and helpers live under `tests/benchmarking/`; the
`benchmarking/` tree contains only benchmark implementation, profiles,
adapters, documentation, and reproducers. The suite is also registered in the
ordinary test hierarchy. `make test` builds the three release adapters and
runs a bounded real-MPI integration gate: DG at `N=4` and `N=8` on each
adapter, plus one two-rank full-field comparison. This is a correctness test
with no timing assertion and no published benchmark run. It does not execute a
benchmark profile; set `SKIP_MPI_CASES=1` only when MPI execution is
intentionally unavailable.

The remaining dependency-free contract tests cover the 792-case temporal plan,
all six spatial profiles, both validation matrices, both strong-scaling
matrices and their exact case counts, all three weak profiles and their
measurement/helper counts, the final 792/12,960/14,985 configuration audit
across all three problems, DG/PR/BE, eight temporal levels, 11 temporal and 15
scaling degree pairs, strong/weak, X/Y/Z decomposition, and OMP `1,2,4,8`, plus case-ID
collision and controlled-run-path checks, execution provenance and explicit
unavailable observations,
resume build identity, failure classification and bounded retry histories,
numerics-first A/B policy, detached-worktree ownership and alternating
case scheduling, weak-efficiency and same-global-mesh field gates, launcher-template
expansion and resource rejection, deterministic index/group sharding, and
missing/duplicate/conflicting shard rejection, stable
vector-preserving identifiers,
filters and validation, exact-time conversion, safe storage, process timeout
and termination, frozen-manifest decoding, resume compatibility and verified
skip/retry behavior, strict tagged-result parsing, planned/result binding, the
registered manufactured adapters, smoke-refinement verification, synthetic
temporal, spatial, full-field, strong-scaling, and weak-scaling series, raw
repeated samples, MAD/speedup/efficiency, short-region rejection, local
perturbations, separation failures, plateau and corrupted convergence inputs,
fake-adapter extensibility,
and marker-guarded Make cleanup. The suite does not launch the costly full
numerical matrices. The ordinary
`make test` target does not execute full benchmark profiles; real benchmark
runs remain explicit, and the separate short `make test-performance` gate is
not a substitute for them.
