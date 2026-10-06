# Building, dependencies, and generated files

[Documentation index](README.md) · [Repository README](../README.md)

## Dependencies

The supported configuration for current work expects:

- GCC/GFortran through an MPI Fortran wrapper
- MUMPS
- BLAS
- LAPACK
- ScaLAPACK
- METIS
- GKlib

Intel compiler configurations and ParMETIS are outside the current supported
scope. Legacy files or variables referring to them are not qualified build
paths and should not be selected for new work.

The complete test procedure additionally expects:

- pFUnit 4 (for the pFUnit-based test suites)
- Python 3.10 or newer (for layout checks and numerical VTI integration checks)
- Bash and GNU `timeout` (for CLI and MPI error-path tests)

`make test-coverage` additionally requires GNU Fortran with `gcov` and the
`lcov`/`genhtml` tools. Coverage is intentionally a separate target: compiler
instrumentation would invalidate the wall-clock performance gate.

Edit `m_options` for local library paths and compiler flags.

### Supported local toolchain

The supported local GCC/MUMPS paths are:

```text
MPI compiler: /opt/lib/mpich-5.0.0/bin/mpif90
MUMPS:        /opt/lib/MUMPS_5.8.2/
BLAS:         /opt/lib/lapack-3.12.1/lib64/libblas.a
LAPACK:       /opt/lib/lapack-3.12.1/lib64/liblapack.a
ScaLAPACK:    /opt/lib/scalapack-2.3.2/lib64/libscalapack.a
METIS:        /opt/lib/metis/lib/libmetis.a
GKlib:        /opt/lib/GKlib/lib64/libGKlib.a
```

The code is still an MPI program. Directional linear solves are executed with
MUMPS on `MPI_COMM_SELF`, so each rank solves its local gathered systems while
the ADS workflow itself remains distributed through MPI.

The current debug-style options use bounds checking and AddressSanitizer:

```make
-O0 -g -fcheck=all -fbounds-check -fsanitize=address
```

## Building

The repository-root `Makefile` is the public build interface. Run `make help`
to list every target and `make show-config` to display the effective tool and
library paths.

The tracked repository-root `m_options` is the single active configuration
file for the compiler, build profile, MPI launcher, numerical libraries,
pFUnit, Python, Doxygen, and test timeouts. The old `mymake/m_options` path is
retained as a compatibility wrapper which includes that same root file, so
there are no duplicate settings to keep synchronized. Command-line
assignments still override values for one invocation:

```bash
make show-config
make BUILD=release
make MPIEXEC=/path/to/mpiexec MPIFC=/path/to/mpif90 test
```

Ready-to-select examples live in `makeconfig/`:

```bash
make list-configs
make CONFIG=makeconfig/gnu-debug.mk show-config
make CONFIG=makeconfig/gnu-release.mk all
```

The GNU examples use `mpif90`/`mpiexec` from `PATH`. Numerical-library paths
remain centralized in root `m_options`.

The default `BUILD=debug` preserves the existing bounds checks and
AddressSanitizer flags. `BUILD=release` selects `RELEASE_OPTS` from the same
configuration file.

The make hierarchy is executable at every level:

```text
Makefile
|-- src/GNUmakefile -- src/sources.mk
|-- problems/GNUmakefile -- problems/problems.mk
|   `-- problems/<problem>/GNUmakefile
`-- tests/GNUmakefile
    |-- tests/src/GNUmakefile ------ tests/<source-suite>/GNUmakefile
    |-- tests/problems/GNUmakefile - tests/<problem-suite>/GNUmakefile
    |-- tests/driver/GNUmakefile --- tests/driver_cli/GNUmakefile
    `-- tests/build/GNUmakefile ---- tests/make_hierarchy/GNUmakefile
```

The root does not own core or problem source lists, per-problem defaults, or
test-suite registration. Consequently the corresponding subtree can also be
used directly, for example:

```bash
make -C src library
make -C problems build PROBLEM=heat
make -C problems/heat show-run
make -C tests/src run-suite TEST_SUITE=rhs_assembly
```

### Build targets

The default command builds the static ADS library and all ten problem
drivers. Problem builds are deliberately serialized because problem-local
Fortran modules reuse the names `input_data`, `RHS_fun`, and `main`.

```bash
# Static library plus every problem.
make
make all
make build-all

# Individual layers.
make library
make problems

# One selected problem: equivalent forms.
make build PROBLEM=heat
make build-heat
make heat

make list-problems
```

Generated files are placed in:

```text
mymake/LIB/libads.a
mymake/EXEC/l2
mymake/EXEC/heat
mymake/EXEC/eriksson
mymake/EXEC/pure_diffusion_igrm
mymake/EXEC/oil
mymake/EXEC/igrm_l2
mymake/EXEC/igrm_heat
mymake/EXEC/igrm_eirksson
mymake/EXEC/igrm_stokes
mymake/EXEC/igrm_pollution
```

The lower-level `mymake/makefile` remains usable as a compatibility entry
point. It delegates `library` to `src/` and problem targets to `problems/`; it
no longer owns a duplicate source list or compilation implementation.
Historical named `EXEC=<problem>` builds remain supported. Explicit
`SOURCE_ALL` builds use a deliberately isolated compatibility adapter: it
compiles only the caller-supplied ordered list, puts its objects and modules in
`<BUILD_ROOT>/_LEGACY_OBJ/<EXEC>` (by default below `mymake/`), and never
reuses the official core or problem object directories. Thus legacy custom
programs remain buildable without moving source ownership back into `mymake`.

Both the core and each problem keep a stamp of their effective compiler and
flags. A debug/release or compiler/flag switch therefore rebuilds incompatible
objects, while an unchanged second invocation remains incremental.

### Cleanup targets

```bash
make clean-build
make clean-problems
make clean-library
make clean-legacy
make clean-tests
make clean-coverage
make clean-docs
make clean-benchmark-build
make clean          # all generated build, test, and documentation artifacts
make distclean      # same as clean; m_options and makeconfig/ are preserved
```

Benchmark results under `benchmarks/<run-id>/` are never removed by these
targets. `clean-benchmark-build` owns only framework bytecode and a guarded
`benchmarking/build/` directory. Run manifests, logs, fields, and reports are
ignored, unversioned runtime data; there is no broad results-clean target.

The CMake files are still present, but this README documents the current
make-based workflow.

## Documentation

Documentation uses the tracked `doxygen.conf`. Doxygen must run from the
repository root because its input and output paths are relative.

```bash
make docs-check   # validate/expand the Doxygen configuration
make docs-html    # doxygen/html/index.html
make docs-pdf     # HTML plus doxygen/latex/refman.pdf
make docs         # same complete HTML+PDF documentation build
make clean-docs
```

`DOXYGEN` can be overridden in `m_options` or on the command line. PDF
generation additionally requires a LaTeX installation with `pdflatex` and
`makeindex`.

## Notes

- `mymake/EXEC/`, `mymake/LIB/`, and `mymake/_OBJ/` contain the public
  generated executables, library, and core objects; each problem keeps its
  own generated `_OBJ/` directory. Explicit legacy `SOURCE_ALL` builds use
  marker-protected `<BUILD_ROOT>/_LEGACY_OBJ/<EXEC>` directories.
- Doxygen-style comments are used throughout `src` and the migrated problem
  drivers.
