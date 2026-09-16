.DEFAULT_GOAL := all

.NOTPARALLEL:

.PHONY: all tests run run-src run-problems run-driver run-build-system \
	run-cli run-smoke run-integration run-performance \
	run-performance-self-test run-coverage clean-performance \
	clean-coverage-data run-suite check-layout list clean

CONFIG      ?= ../m_options
CONFIG_PATH := $(abspath $(CONFIG))
include $(CONFIG_PATH)

BUILD_ROOT  ?= ../mymake
BUILD_ROOT_PATH := $(abspath $(BUILD_ROOT))

PYTHON      ?= python3
PFUNIT_ROOT ?= /opt/lib/pfunit/PFUNIT-4.16
MPIEXEC     ?= /opt/lib/mpich-5.0.0/bin/mpiexec
MPIFC       ?= /opt/lib/mpich-5.0.0/bin/mpif90
MPIEXEC_FLAGS ?=
MPI_NP_FLAG ?= -n
OMP_PROC_BIND ?= close
OMP_PLACES ?= cores
MUMPS_DIR   ?= /opt/lib/MUMPS_5.8.2
SUITE_TIMEOUT ?= 600s
DRIVER_CLI_TIMEOUT ?= 20s
DRIVER_SMOKE_TIMEOUT ?= 60s
DRIVER_INTEGRATION_TIMEOUT ?= 90
SKIP_MPI_CASES ?= 0
PERFORMANCE_BUILD_ROOT ?= $(abspath ../build/openmp-performance)
PERFORMANCE_BUILD_ROOT_PATH := $(abspath $(PERFORMANCE_BUILD_ROOT))
PERFORMANCE_TIMEOUT ?= 300
PERFORMANCE_SUITE_TIMEOUT ?= 3600s
PERFORMANCE_WARMUPS ?= 1
PERFORMANCE_SAMPLES ?= 3
PERFORMANCE_MIN_SPEEDUP ?= 1.10
PERFORMANCE_MAX_REGRESSION ?= 1.15
PERFORMANCE_BASELINE ?=
PERFORMANCE_BASELINE_PATH := $(if $(strip $(PERFORMANCE_BASELINE)),$(abspath $(PERFORMANCE_BASELINE)),)
COVERAGE_FLAGS ?=
COVERAGE_COMPILER = $(MPIFC) --coverage
COVERAGE_MAKEFILE := $(abspath make/coverage.mk)

# This runner knows only the test groups.  Suite registration and execution
# belong to the group makefiles, while each suite owns its actual build rules.
GROUPS := src problems driver build

GROUP_OPTIONS = \
	CONFIG="$(CONFIG_PATH)" \
	BUILD_ROOT="$(BUILD_ROOT_PATH)" \
	PFUNIT_ROOT="$(PFUNIT_ROOT)" \
	MPIEXEC="$(MPIEXEC)" \
	MPIFC="$(MPIFC)" \
	FC="$(MPIFC)" \
	MUMPS_DIR="$(MUMPS_DIR)" \
	PYTHON="$(PYTHON)" \
	SUITE_TIMEOUT="$(SUITE_TIMEOUT)" \
	DRIVER_CLI_TIMEOUT="$(DRIVER_CLI_TIMEOUT)" \
	DRIVER_SMOKE_TIMEOUT="$(DRIVER_SMOKE_TIMEOUT)" \
	DRIVER_INTEGRATION_TIMEOUT="$(DRIVER_INTEGRATION_TIMEOUT)" \
	SKIP_MPI_CASES="$(SKIP_MPI_CASES)" \
	PERFORMANCE_BUILD_ROOT="$(PERFORMANCE_BUILD_ROOT_PATH)" \
	PERFORMANCE_TIMEOUT="$(PERFORMANCE_TIMEOUT)" \
	PERFORMANCE_SUITE_TIMEOUT="$(PERFORMANCE_SUITE_TIMEOUT)" \
	PERFORMANCE_WARMUPS="$(PERFORMANCE_WARMUPS)" \
	PERFORMANCE_SAMPLES="$(PERFORMANCE_SAMPLES)" \
	PERFORMANCE_MIN_SPEEDUP="$(PERFORMANCE_MIN_SPEEDUP)" \
	PERFORMANCE_MAX_REGRESSION="$(PERFORMANCE_MAX_REGRESSION)" \
	PERFORMANCE_BASELINE="$(PERFORMANCE_BASELINE_PATH)" \
	MPIEXEC_FLAGS="$(MPIEXEC_FLAGS)" \
	MPI_NP_FLAG="$(MPI_NP_FLAG)" \
	OMP_PROC_BIND="$(OMP_PROC_BIND)" \
	OMP_PLACES="$(OMP_PLACES)" \
	COVERAGE_FLAGS="$(COVERAGE_FLAGS)"

define delegate_groups
	@set -eu; \
	for group in $(1); do \
		printf '\n==> tests/%s: %s\n' "$$group" "$(2)"; \
		$(MAKE) --no-print-directory -j1 -C "$$group" \
			$(GROUP_OPTIONS) "$(2)"; \
	done
endef

check-layout:
	$(PYTHON) ./check_test_layout.py

all: check-layout
	$(call delegate_groups,$(GROUPS),all)

tests: all

run-src: check-layout
	+$(MAKE) --no-print-directory -j1 -C src $(GROUP_OPTIONS) run

run-problems: check-layout
	+$(MAKE) --no-print-directory -j1 -C problems $(GROUP_OPTIONS) run

run-driver:
	+$(MAKE) --no-print-directory -j1 -C driver $(GROUP_OPTIONS) run

run-build-system: check-layout
	+$(MAKE) --no-print-directory -j1 -C build $(GROUP_OPTIONS) run

run-cli:
	+$(MAKE) --no-print-directory -j1 -C driver $(GROUP_OPTIONS) run-cli

run-smoke:
	+$(MAKE) --no-print-directory -j1 -C driver $(GROUP_OPTIONS) run-smoke

run-integration: check-layout
	+$(MAKE) --no-print-directory -j1 -C driver \
		$(GROUP_OPTIONS) run-integration

run-performance:
	+$(MAKE) --no-print-directory -j1 -C driver \
		$(GROUP_OPTIONS) run-performance

run-performance-self-test:
	+$(MAKE) --no-print-directory -j1 -C driver \
		$(GROUP_OPTIONS) run-performance-self-test

run-coverage: check-layout clean-coverage-data
	@if test -z '$(strip $(COVERAGE_FLAGS))'; then \
		printf 'COVERAGE_FLAGS must be set for run-coverage\n' >&2; exit 2; \
	fi
	+$(MAKE) --no-print-directory -j1 -C build \
		$(GROUP_OPTIONS) run
	+env MAKEFILES='$(COVERAGE_MAKEFILE)' $(MAKE) --no-print-directory -j1 \
		-C src $(GROUP_OPTIONS) MPIFC='$(COVERAGE_COMPILER)' \
		FC='$(COVERAGE_COMPILER)' run
	+env MAKEFILES='$(COVERAGE_MAKEFILE)' $(MAKE) --no-print-directory -j1 \
		-C problems $(GROUP_OPTIONS) MPIFC='$(COVERAGE_COMPILER)' \
		FC='$(COVERAGE_COMPILER)' run
	+env MAKEFILES='$(COVERAGE_MAKEFILE)' $(MAKE) --no-print-directory -j1 \
		-C driver $(GROUP_OPTIONS) MPIFC='$(COVERAGE_COMPILER)' \
		FC='$(COVERAGE_COMPILER)' run-coverage

clean-performance:
	+$(MAKE) --no-print-directory -j1 -C driver \
		$(GROUP_OPTIONS) clean-performance

run-suite:
	@if [ -z "$(TEST_SUITE)" ]; then \
		printf '%s\n' \
			'TEST_SUITE is required, for example: make run-suite TEST_SUITE=rhs_assembly'; \
		exit 2; \
	fi
	@set -eu; \
	found=0; \
	available=''; \
	for group in $(GROUPS); do \
		suites="$$( $(MAKE) --no-print-directory -s -C "$$group" list-inline )"; \
		available="$$available $$suites"; \
		case " $$suites " in \
			*" $(TEST_SUITE) "*) \
				$(MAKE) --no-print-directory -j1 -C "$$group" \
					$(GROUP_OPTIONS) TEST_SUITE="$(TEST_SUITE)" run-suite; \
				found=1; \
				break ;; \
		esac; \
	done; \
	if [ "$$found" -eq 0 ]; then \
		printf 'Unknown test suite: %s\nAvailable: %s\n' \
			'$(TEST_SUITE)' "$${available# }"; \
		exit 2; \
	fi

run: run-src run-problems run-driver run-build-system

list:
	@set -eu; \
	for group in $(GROUPS); do \
		$(MAKE) --no-print-directory -s -C "$$group" list; \
	done

clean-coverage-data:
	@find '$(CURDIR)' -type f \( -name '*.gcda' -o -name '*.gcno' \) -delete

clean: clean-coverage-data
	$(call delegate_groups,$(GROUPS),clean)
