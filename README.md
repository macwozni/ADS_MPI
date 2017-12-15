IGA-ADS
=======

IGA-ADS is a Fortran framework designed to facilitate creating parallel numerical simulations for time-dependent PDEs using isogeometric finite element method.


Requirements
------------

1. Dependencies
- BLAS
- LAPACK
- MPI

2. Tools
- compiler: reasonable Fortran2008 support is required (framework has been tested with ifort 14.0.2 and gfortran 4.8.5)
- build system: CMake 3.1 or higher

3. Compilation

To compile the code, create a directory for the build and execute following commands:

\$ cmake <options> \${PROJECT_DIR}
\$ make

where \${PROJECT_DIR} is the root directory of the project source. Options allow customizing which parts of the project are compiled. By default example applications are compiled.

- SKIP_PROBLEMS - if ON, the example problems are not compiled (default: OFF)
- COMPILE_TESTS - wether the unit tests are compiled. This should be disabled if pfUNIT is not available (default: OFF)

Options are specified as -Doption=value, e.g. 
\$ cmake -DSKIP_PROBLEMS=ON ..

---

Contents
--------

### Top-level structure:
  src/      - framework code
  problems/  - example problem implementations
  test/         - unit tests

### Detailed description of the package contents:

CMake build configuration:
  CMakeLists.txt

1. Simulation infrastructure:
 - src/knot_vector.F90
 - src/parallelism.F90
 -  src/communicators.F90
 - src/gauss.F90
 - src/math.F90
 - src/reorderRHS.F90
 - src/my_mpi.F90
 - src/projection_engine.F90
 - src/analysis.F90
 - src/int2str.F90
 - src/plot.F90

2. ADS solver implementation:
 - src/ADS.F90
 
3. B-spline functions implementation:
 - src/basis.F90

4. File output:
 - src/gnuplot.F90
 - src/vtk.F90

5. Problem implementations:

6. Unit tests:

