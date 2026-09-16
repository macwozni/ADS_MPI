# ADS MPI benchmarking framework

This directory is the tracked, reusable benchmark infrastructure. Generated
plans, logs, fields, and measurements belong under the ignored
`benchmarks/<run-id>/` tree. The pre-existing
`benchmarks/igrm_strong_scaling/` prototype is user data and is never adopted,
overwritten, moved, or cleaned by this framework.

Stage 1 implements strict configuration, deterministic planning, provenance,
safe result storage, and the generic process contract. It deliberately does
not build or run an ADS solver. The three real adapters are registered as
planning-only until the manufactured transient harness is added in stage 2.

## Architecture and dependency direction

```text
configs/*.json
       |
       v
catalog.py ----> components/planning.py
       |                    |
       v                    v
framework/{config,planner,registry,protocols,executor,storage,...}
       |
       v
benchmarks/<run-id>/...
```

The modules in `ads_benchmark/framework/` never import concrete problems.
They receive all extension points through `Catalog` registries. Concrete
adapters and launchers import the framework contracts; `catalog.py` is the
single composition root that registers them. There is no problem-name or
family `if/elif` chain in the planner or executor.

The current extension axes are:

- problem adapter: validates problem details, builds a payload argv, and
  parses one domain result;
- experiment family: `temporal`, `h`, `p`, `validation`, `strong`, or `weak`;
- exact/manufactured case;
- build profile;
- launcher.

Resource profiles and analyzers/reporters are separate planned registries, not
problem-adapter branches. Stage 1 records measurement/build/launcher resources
directly in each normalized case and emits a raw common result envelope. The
dedicated resource-profile registry arrives with strong/weak/hybrid execution,
and analyzer/reporter registries with convergence and scaling reports. Those
later components will be injected into the same catalog; the neutral engine
will not import them directly.

The generic executor alone owns the working directory, environment, OpenMP
settings, process lifetime, timeout, stdout/stderr capture, status transitions,
and common result envelope. Each payload runs from its owned
`cases/<case-id>/` directory, so relative solver artifacts remain isolated. A
launcher wraps an adapter's payload argv. Adding a component requires an
implementation plus one registration in `catalog.py`, not an engine change.
`selftests/fake_adapter.py` proves this contract by planning, launching,
checking the working directory, logging, and parsing a real small subprocess
through a fresh registration.

## Profiles and normalized cases

Profiles are strict JSON without third-party Python dependencies. Every
profile explicitly supplies:

- family, problem, scheme, and exact case;
- final time and either a time-step count or exact decimal `time_step`;
- three-dimensional element grid;
- test/trial degrees in all three directions;
- MPI rank count and three-dimensional process grid;
- OpenMP thread count;
- warmups, measured samples, timeout, build profile, and launcher.

Times are parsed with exact rational arithmetic and serialized canonically: a
finite decimal such as `0.025` remains decimal, while a recurring value such as
`0.1/3` becomes the reduced string `1/30`. This avoids binary floating-point
drift and Decimal-context rounding in validation and identifiers.
`final_time / time_step` must be an exact positive integer. Future Fortran
adapters must convert this exact representation deliberately at their runtime
precision; they must not pass a fraction literal such as `1/30` to a solver.

Registered profiles are data in `configs/`:

- `smoke`: 3 problems x 3 schemes = 9 planning cases;
- `temporal-full`: the required 792-case temporal matrix;
- `local-scaling`: a stage-1 local planning preset;
- `cluster-scaling`: a stage-1 distributed planning preset.

The last two presets establish the configuration contract; the complete
scientific strong/weak matrices are intentionally deferred to their later
implementation stages.

`temporal-full` fixes the otherwise unspecified mesh at `4x4x4` and expands:

```text
3 problems x 3 schemes x 8 step counts x 11 degree pairs = 792
```

It uses `T=0.1`, `N=4,...,512`, and the exact degree pairs requested by the
benchmark specification. All other axes are explicit singletons.

Each case is normalized and encoded as sorted compact JSON. Its ID is
`<problem>-<scheme>-<20 hex digits>` from SHA-256 of that semantic document.
Profile name, JSON key/axis order, run ID, timestamps, output paths, Git SHA,
and dirty state do not affect the ID. Duplicate semantic cases are rejected
before filtering rather than silently deduplicated.

## Planning

From the repository root:

```bash
make benchmark-plan
make benchmark-plan BENCHMARK_PROFILE=temporal-full
```

Both commands are dry runs: they validate the complete source matrix, apply
filters, inspect Git provenance, print a summary, and write nothing. The
standalone equivalent is:

```bash
make -C benchmarking plan BENCHMARK_PROFILE=smoke
```

To inspect the complete JSON plan or write a new manifest:

```bash
cd benchmarking
python3 -m ads_benchmark plan --profile smoke --json
python3 -m ads_benchmark plan --profile smoke \
  --write-manifest --run-id smoke-20260916
```

A write is exclusive. An existing run directory, including the legacy
`igrm_strong_scaling` directory, is an error and is never overwritten.

Repeatable, conjunctive filters are available:

```text
--problem igrm_heat
--scheme be
--degree-pair 4:3
--degree-pair 4x5x6:3x4x5
--mesh 4x4x4
--mpi-grid 2x1x1
--mpi-ranks 2
--omp 4
```

Unknown values, invalid filter syntax, and a zero-case result fail explicitly.
The complete profile is validated and duplicate-checked before filters are
applied, so a filter cannot conceal a broken case.

Validation includes registered component names, positive finite times,
integral step counts, degree range `1..9`, `p_test > p_trial` per direction,
positive dimensions, `NP=procx*procy*procz`, and positive OpenMP/measurement
settings. A process-grid axis may not exceed the trial-space DOF count
`nelem_i + p_trial_i`, matching the library's quotient/remainder ownership
partition. Per-case timeouts are capped at 30 days to remain representable by
the subprocess runtime.

## Result layout and safety

A written plan initially creates only:

```text
benchmarks/<run-id>/manifest.json
```

It does not create successful-looking case records. Once execution is enabled,
the common engine owns:

```text
benchmarks/<run-id>/cases/<case-id>/status.json
benchmarks/<run-id>/cases/<case-id>/result.json
benchmarks/<run-id>/cases/<case-id>/stdout.log
benchmarks/<run-id>/cases/<case-id>/stderr.log
```

The manifest contains a schema version, the full expanded configuration and
case IDs, configuration hash, Git commit, and dirty-tree flag. Result writes
are atomic. The store accepts exactly the canonical repository `benchmarks/`
root, safe single-component IDs, and a new run directory. A matching plan
manifest owns only that run; executor APIs refuse to adopt any other directory,
including the legacy prototype. Descriptor-relative POSIX operations with
`O_NOFOLLOW` reject traversal, sibling-prefix tricks, symlink swaps, and
existing targets. No ownership marker ever claims the shared `benchmarks/`
root.

`make clean` and `make clean-benchmark-build` never access benchmark results.
The latter removes Python bytecode and only marker-owned
`benchmarking/build/` content.

## Self-tests

```bash
make benchmark-self-test
```

The dependency-free suite covers the 792-case count, stable identifiers,
duplicate detection, every filter, strict JSON parsing, validation failures,
Git dirty-state provenance, path containment, exclusive manifests, and
successful/nonzero/timeout/parser-failure execution through the fake adapter.
They also exercise an ignoring descendant process, extension-point failures,
legacy-run non-adoption, symlink swaps, and marker-guarded Make cleanup.
It is intentionally outside `tests/`, so costly future benchmark execution
cannot enter the ordinary `make test` regression tree.
