# ADS MPI documentation

[Repository README](../README.md)

This index is the stable entry point for detailed project documentation. All
commands in the linked guides assume the repository root unless a guide says
otherwise.

## Core project

- [Building, dependencies, and generated files](building.md)
- [Core architecture and time integration](architecture.md)
- [Shared problem-running interface](running.md)
- [Classic ADS problem drivers](problems-classic.md)
- [iGRM and DPG problem drivers](problems-igrm.md)

## Verification

- [Test hierarchy and commands](../tests/README.md)
- [Positive numerical integration matrix](../tests/integration.md)
- [OpenMP performance regression](../tests/performance.md)
- [Benchmark framework tests](../tests/benchmarking/README.md)

## Benchmarking

- [Benchmark framework index](../benchmarking/README.md)
- [Benchmark topic guides](../benchmarking/docs/)
- [Scientific reproducers](../benchmarking/reproducers/)

## Generated documentation

```bash
make docs-check
make docs-html
make docs-pdf
```

Generated Doxygen output is written below ignored `doxygen/`; it is not
source documentation and should not be committed.
