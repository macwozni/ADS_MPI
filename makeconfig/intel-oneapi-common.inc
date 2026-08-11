# Shared Intel oneAPI overrides. Included by the two public Intel examples.
_ADS_EXAMPLE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
include $(_ADS_EXAMPLE_DIR)../m_options

COMPILER = mpiifx
MPIFC = $(COMPILER)
MPIEXEC = mpiexec
DIALECT = -stand f03 -implicitnone -fpp -qopenmp -warn all
MODULE_OUTPUT = -module .
BOUNDS = -check bounds
DEBUG_OPTS = -O0 -g -traceback -check all $(BOUNDS)
RELEASE_OPTS = -O3 -unroll
FPP = fpp -O3
