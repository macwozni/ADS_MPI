# GNU Fortran + MPI wrappers, optimized build.
BUILD ?= release
_ADS_EXAMPLE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
include $(_ADS_EXAMPLE_DIR)../m_options

COMPILER = mpif90
MPIFC = $(COMPILER)
MPIEXEC = mpiexec
