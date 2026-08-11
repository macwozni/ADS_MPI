# Intel oneAPI MPI wrapper with optimization enabled.
BUILD ?= release
_ADS_INTEL_ENTRY_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
include $(_ADS_INTEL_ENTRY_DIR)intel-oneapi-common.inc
