# Intel oneAPI MPI wrapper with runtime checks.
BUILD ?= debug
_ADS_INTEL_ENTRY_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
include $(_ADS_INTEL_ENTRY_DIR)intel-oneapi-common.inc
