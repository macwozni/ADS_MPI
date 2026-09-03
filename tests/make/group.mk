# Shared implementation for the three test-group runners.  The including
# makefile supplies TEST_GROUP_NAME and SUITES; individual suites continue to
# own all compilation, execution, and cleanup details.

ifndef TEST_GROUP_NAME
$(error TEST_GROUP_NAME must be set before including tests/make/group.mk)
endif

ifndef SUITES
$(error SUITES must be set before including tests/make/group.mk)
endif

.DEFAULT_GOAL := all

.NOTPARALLEL:

.PHONY: all tests run run-suite list list-inline clean

TESTS_ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/..)
CONFIG     ?= $(TESTS_ROOT)/../m_options
CONFIG_PATH := $(abspath $(CONFIG))
include $(CONFIG_PATH)

BUILD_ROOT ?= $(TESTS_ROOT)/../mymake
BUILD_ROOT_PATH := $(abspath $(BUILD_ROOT))

PFUNIT_ROOT ?= /opt/lib/pfunit/PFUNIT-4.16
MPIEXEC     ?= /opt/lib/mpich-5.0.0/bin/mpiexec
MPIFC       ?= /opt/lib/mpich-5.0.0/bin/mpif90
MUMPS_DIR   ?= /opt/lib/MUMPS_5.8.2
PYTHON      ?= python3
SUITE_TIMEOUT ?= 600s
DRIVER_CLI_TIMEOUT ?= 20s
DRIVER_SMOKE_TIMEOUT ?= 60s
DRIVER_INTEGRATION_TIMEOUT ?= 90
SKIP_MPI_CASES ?= 0

SUITE_OPTIONS = \
	CONFIG="$(CONFIG_PATH)" \
	BUILD_ROOT="$(BUILD_ROOT_PATH)" \
	PYTHON="$(PYTHON)" \
	PFUNIT_ROOT="$(PFUNIT_ROOT)" \
	MPIEXEC="$(MPIEXEC)" \
	MPIFC="$(MPIFC)" \
	FC="$(MPIFC)" \
	DRIVER_CLI_TIMEOUT="$(DRIVER_CLI_TIMEOUT)" \
	DRIVER_SMOKE_TIMEOUT="$(DRIVER_SMOKE_TIMEOUT)" \
	DRIVER_INTEGRATION_TIMEOUT="$(DRIVER_INTEGRATION_TIMEOUT)" \
	SKIP_MPI_CASES="$(SKIP_MPI_CASES)" \
	MUMPS_DIR="$(MUMPS_DIR)"

define run_suites
	@set -eu; \
	for suite in $(1); do \
		printf '\n==> tests/%s [%s]: %s\n' \
			"$$suite" "$(TEST_GROUP_NAME)" "$(2)"; \
		timeout $(SUITE_TIMEOUT) $(MAKE) --no-print-directory -j1 \
			-C "$(TESTS_ROOT)/$$suite" $(SUITE_OPTIONS) "$(2)"; \
	done
endef

define rebuild_and_run_suites
	@set -eu; \
	for suite in $(1); do \
		printf '\n==> tests/%s [%s]: clean + %s\n' \
			"$$suite" "$(TEST_GROUP_NAME)" \
			"$(if $(strip $(2)),$(2),run)"; \
		timeout $(SUITE_TIMEOUT) $(MAKE) --no-print-directory -j1 \
			-C "$(TESTS_ROOT)/$$suite" $(SUITE_OPTIONS) clean; \
		timeout $(SUITE_TIMEOUT) $(MAKE) --no-print-directory -j1 \
			-C "$(TESTS_ROOT)/$$suite" $(SUITE_OPTIONS) \
			"$(if $(strip $(2)),$(2),run)"; \
	done
endef

all:
	$(call run_suites,$(SUITES),all)

tests: all

run:
	$(call rebuild_and_run_suites,$(SUITES))

run-suite:
	@if [ -z "$(TEST_SUITE)" ]; then \
		printf '%s\n' \
			'TEST_SUITE is required, for example: make run-suite TEST_SUITE=rhs_assembly'; \
		exit 2; \
	fi
	@case ' $(SUITES) ' in \
		*' $(TEST_SUITE) '*) ;; \
		*) printf 'Test suite %s is not in the %s group. Available: %s\n' \
			'$(TEST_SUITE)' '$(TEST_GROUP_NAME)' '$(SUITES)'; exit 2 ;; \
	esac
	$(call rebuild_and_run_suites,$(TEST_SUITE))

list:
	@printf '%s\n' $(SUITES)

list-inline:
	@printf '%s\n' '$(SUITES)'

clean:
	$(call run_suites,$(SUITES),clean)
