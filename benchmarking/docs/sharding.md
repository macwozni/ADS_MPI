# Sharded benchmark runs

[Benchmarking index](../README.md) · [Repository README](../../README.md)

Large weak (or other) frozen plans can be split without editing their
configuration. `problem-scheme-degree` creates one deterministic shard per
problem/scheme/test-degree/trial-degree group. `index` assigns the sorted case
IDs round-robin to an explicit shard count:

```bash
# Public root-Make forms.
make benchmark-shard-plan RUN_ID=weak-parent-grouped \
  BENCHMARK_SHARD_PROFILE=cluster-weak-scaling \
  BENCHMARK_SHARD_STRATEGY=problem-scheme-degree

make benchmark-shard-plan RUN_ID=weak-parent \
  BENCHMARK_SHARD_PROFILE=cluster-weak-scaling \
  BENCHMARK_SHARD_STRATEGY=index BENCHMARK_SHARD_COUNT=3
```

The lower-level CLI equivalents make every frozen input explicit:

```bash
cd benchmarking

# Create the immutable parent manifest and grouped shard envelopes.
python3 -m ads_benchmark.cli shard \
  --repository-root .. --config-dir configs \
  --profile cluster-weak-scaling --run-id weak-parent-grouped \
  --strategy problem-scheme-degree

# Alternative used below: exactly three index shards.
python3 -m ads_benchmark.cli shard \
  --repository-root .. --config-dir configs \
  --profile cluster-weak-scaling --run-id weak-parent \
  --strategy index --shard-count 3
```

The parent owns `benchmarks/<parent>/manifest.json` and numbered immutable
envelopes below `benchmarks/<parent>/shards/`. Each envelope binds the complete
parent metadata, repository/source identity, full expected case digest, shard
index/count, and subset digest. Generation refuses an existing parent.

After building the release adapters, execute each shard into its own run ID.
The command below illustrates shard zero; allocation values must cover that
shard's largest case:

```bash
make benchmark-run-shard \
  PARENT_RUN_ID=weak-parent SHARD_INDEX=0 \
  RUN_ID=weak-parent-part-000 \
  BENCHMARK_LAUNCHER_TEMPLATE='srun --ntasks={ranks} --cpus-per-task={threads} {payload}' \
  BENCHMARK_PLAN_ARGS='--available-mpi-slots 64 --available-cpu-slots 512'
```

The public target builds the release adapters. Its lower-level equivalent,
when already inside `benchmarking/`, is:

```bash
make build BUILD=release

python3 -m ads_benchmark.cli run-shard \
  --repository-root .. \
  --parent-run-id weak-parent --shard-index 0 \
  --run-id weak-parent-part-000 \
  --available-mpi-slots 64 --available-cpu-slots 512 \
  --launcher-template \
    'srun --ntasks={ranks} --cpus-per-task={threads} {payload}'
```

`run-shard --resume` verifies and resumes that isolated shard run by the same
rules as an ordinary frozen run. Once every shard has completed, list every
child run explicitly and merge into the parent analysis:

```bash
make benchmark-merge-shards PARENT_RUN_ID=weak-parent \
  SHARD_RUNS='weak-parent-part-000 weak-parent-part-001 weak-parent-part-002' \
  BENCHMARK_ANALYZE_ARGS=--plot
```

The direct equivalent is:

```bash
python3 -m ads_benchmark.cli merge-shards \
  --repository-root .. --parent-run-id weak-parent \
  --shard-run weak-parent-part-000 \
  --shard-run weak-parent-part-001 \
  --shard-run weak-parent-part-002 \
  --plot
```

This merge matches the three-shard index parent created above. Merge refuses a
missing shard/index or case, a duplicated case or child contribution,
conflicting case content, incompatible config hashes, and incompatible
repository/source identities. It analyzes only after the contributed child
manifests reconstruct the exact parent plan; it does not silently accept a
partial report. Generated results and shard manifests remain ignored runtime
data and are not committed with the framework.
