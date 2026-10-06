# Running problem drivers

[Documentation index](README.md) · [Repository README](../README.md)

All commands in this guide are run from the repository root.

Every problem has a root-level `run-<problem>` target. The generic form is:

```bash
make run PROBLEM=heat
make run-heat
make run-heat ARGS='4 2 3 0.01 2 1 1' NP=2 OMP_NUM_THREADS=4
```

Use `make run-help` to print the complete argument syntax and defaults for all
ten problems. `make show-run PROBLEM=heat` prints the effective executable,
arguments, MPI/OpenMP settings, environment, and output directory without
building or launching it.

If `ARGS` is omitted, each target uses a small one-rank example with
`steps=1` for transient problems. Output is written to `output/<problem>` by
default; set `RUN_DIR` to choose another directory.
The historical per-problem form, for example `heat_ARGS='...'`, is also
accepted when the common `ARGS` variable is empty.

All runtime controls can be supplied on the make command line:

```text
ARGS                 raw problem command-line arguments
NP                   MPI rank count, default 1
MPIEXEC              MPI launcher from m_options
MPIEXEC_FLAGS        additional launcher flags
MPI_NP_FLAG          rank-count flag, default -n
OMP_NUM_THREADS      OpenMP thread count, default 1
OMP_DYNAMIC          OpenMP dynamic teams, default FALSE
OMP_PROC_BIND        OpenMP binding policy, default close
RUN_DIR              output working directory, default output/<problem>
RUN_ENV              additional environment assignments
OIL_SEED             shortcut for ADS_OIL_RANDOM_SEED on run-oil
```

For example:

```bash
make run-oil OIL_SEED=20260811 OMP_NUM_THREADS=4
make run-igrm_heat ARGS='2 2 2 3 3 3 2 2 2 2 1 1 1 0.001 be' NP=2
make run-igrm_eirksson NP=1
make run-igrm_stokes NP=1
make run-igrm_pollution NP=1
```

The number of MPI ranks must match:

```text
procx * procy * procz
```

The executables may also be invoked directly:

```bash
/opt/lib/mpich-5.0.0/bin/mpiexec -n 1 ./mymake/EXEC/l2 2 2 2 1 1 1 1
```

## Driver reference

- [Classic ADS drivers](problems-classic.md)
- [iGRM and DPG drivers](problems-igrm.md)
