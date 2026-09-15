.DEFAULT_GOAL := all

.NOTPARALLEL:

.PHONY: all tests run run-cli run-smoke run-integration \
	performance-build run-performance run-performance-self-test \
	validate-performance-paths prepare-performance-root \
	validate-performance-ownership clean-performance clean

CONFIG     ?= ../../m_options
CONFIG_PATH := $(abspath $(CONFIG))
include $(CONFIG_PATH)

MYMAKE_DIR ?= ../../mymake
SRC_DIR    ?= ../../src
PROBLEMS_DIR ?= ../../problems
MYMAKE_PATH := $(abspath $(MYMAKE_DIR))
SRC_PATH := $(abspath $(SRC_DIR))
PROBLEMS_PATH := $(abspath $(PROBLEMS_DIR))
REPOSITORY_ROOT := $(abspath $(PROBLEMS_PATH)/..)
TESTS_PATH := $(abspath ..)
BUILD_ROOT ?= $(MYMAKE_PATH)
BUILD_ROOT_PATH := $(abspath $(BUILD_ROOT))
include $(PROBLEMS_PATH)/problems.mk

MPIEXEC    ?= /opt/lib/mpich-5.0.0/bin/mpiexec
MPIFC      ?= /opt/lib/mpich-5.0.0/bin/mpif90
MPIEXEC_FLAGS ?=
MPI_NP_FLAG ?= -n
OMP_PROC_BIND ?= close
OMP_PLACES ?= cores
PYTHON     ?= python3
DRIVER_CLI_TIMEOUT ?= 20s
DRIVER_SMOKE_TIMEOUT ?= 60s
DRIVER_INTEGRATION_TIMEOUT ?= 90
SKIP_MPI_CASES ?= 0
EXEC_DIR   := $(BUILD_ROOT_PATH)/EXEC
ORACLE_OBJ_DIR := $(BUILD_ROOT_PATH)/driver_cli_test_OBJ
ORACLE_SOURCE := pure_diffusion_nonzero_oracle.F90
ORACLE_OBJECT := $(ORACLE_OBJ_DIR)/pure_diffusion_nonzero_oracle.o
ORACLE_MODULE := $(ORACLE_OBJ_DIR)/pure_diffusion_oracle_rhs.mod
ORACLE_EXECUTABLE := $(EXEC_DIR)/pure_diffusion_nonzero_oracle
PERFORMANCE_BUILD_ROOT ?= $(REPOSITORY_ROOT)/build/openmp-performance
PERFORMANCE_BUILD_ROOT_PATH := $(abspath $(PERFORMANCE_BUILD_ROOT))
override PERFORMANCE_OBJ_DIR := $(PERFORMANCE_BUILD_ROOT_PATH)/igrm_l2_OBJ
override PERFORMANCE_EXECUTABLE := $(PERFORMANCE_BUILD_ROOT_PATH)/EXEC/igrm_l2
override PERFORMANCE_RESULT := $(PERFORMANCE_BUILD_ROOT_PATH)/openmp-performance.json
override PERFORMANCE_MARKER := $(PERFORMANCE_BUILD_ROOT_PATH)/.ads-openmp-performance-root
PERFORMANCE_TIMEOUT ?= 300
PERFORMANCE_WARMUPS ?= 1
PERFORMANCE_SAMPLES ?= 3
PERFORMANCE_MIN_SPEEDUP ?= 1.10
PERFORMANCE_MAX_REGRESSION ?= 1.15
PERFORMANCE_BASELINE ?=
PROBLEM_BUILD_OPTIONS := CONFIG='$(CONFIG_PATH)' \
	                         BUILD_ROOT='$(BUILD_ROOT_PATH)' \
	                         COMPILER='$(MPIFC)' MPIFC='$(MPIFC)'

MODULE_OUTPUT ?= -J.
ifneq ($(filter -J.%,$(MODULE_OUTPUT)),)
ORACLE_MODULE_FLAGS := -J$(ORACLE_OBJ_DIR)
else ifneq ($(filter -module,$(MODULE_OUTPUT)),)
ORACLE_MODULE_FLAGS := -module $(ORACLE_OBJ_DIR)
else
ORACLE_MODULE_FLAGS := $(MODULE_OUTPUT)
endif
ORACLE_FFLAGS := $(DIALECT) $(BUILD_OPTS) $(INCLUDE) \
	$(ORACLE_MODULE_FLAGS) -I$(ORACLE_OBJ_DIR) -I$(BUILD_ROOT_PATH)

COMMON_SOURCES := $(wildcard $(SRC_PATH)/*.F90)
BUILD_INPUTS   := GNUmakefile $(COMMON_SOURCES) $(SRC_PATH)/GNUmakefile \
	                  $(SRC_PATH)/sources.mk $(PROBLEMS_PATH)/GNUmakefile \
	                  $(PROBLEMS_PATH)/problem.mk $(PROBLEMS_PATH)/problems.mk \
	                  $(CONFIG_PATH)

DRIVERS := $(addprefix $(EXEC_DIR)/,$(PROBLEMS))

# A driver test is also a clean-build test.  Every executable is rebuilt via
# its owning problem makefile, with the test compiler/configuration forwarded.
.PHONY: $(DRIVERS)

define DEFINE_DRIVER_BUILD
$(EXEC_DIR)/$(1): $(BUILD_INPUTS) \
		$(wildcard $(PROBLEMS_PATH)/$($(1)_DIR)/*.F90) \
		$(PROBLEMS_PATH)/$($(1)_DIR)/GNUmakefile
	+$$(MAKE) --no-print-directory -j1 -C $$(PROBLEMS_PATH) \
		$$(PROBLEM_BUILD_OPTIONS) clean-$(1)
	+$$(MAKE) --no-print-directory -j1 -C $$(PROBLEMS_PATH) \
		$$(PROBLEM_BUILD_OPTIONS) build-$(1)
endef

$(foreach problem,$(PROBLEMS),$(eval $(call DEFINE_DRIVER_BUILD,$(problem))))

$(ORACLE_OBJ_DIR) $(EXEC_DIR):
	mkdir -p $@

# Depending on one real driver ensures that the core library and its module
# files were rebuilt with the exact compiler/configuration used by this suite.
$(ORACLE_OBJECT): $(ORACLE_SOURCE) $(firstword $(DRIVERS)) GNUmakefile \
		$(CONFIG_PATH) | $(ORACLE_OBJ_DIR)
	$(MPIFC) $(ORACLE_FFLAGS) -c $< -o $@

$(ORACLE_EXECUTABLE): $(ORACLE_OBJECT) | $(EXEC_DIR)
	$(MPIFC) $(ORACLE_FFLAGS) -o $@ $< \
		$(BUILD_ROOT_PATH)/LIB/libads.a $(USER_LIB)

all: $(DRIVERS) $(ORACLE_EXECUTABLE)

tests: all

run-cli: all
	MPIEXEC='$(MPIEXEC)' EXEC_DIR='$(abspath $(EXEC_DIR))' \
		DRIVER_CLI_TIMEOUT='$(DRIVER_CLI_TIMEOUT)' \
		SKIP_MPI_CASES='$(SKIP_MPI_CASES)' bash ./test_driver_cli.sh

run-smoke: all
	MPIEXEC='$(MPIEXEC)' EXEC_DIR='$(abspath $(EXEC_DIR))' \
		DRIVER_SMOKE_TIMEOUT='$(DRIVER_SMOKE_TIMEOUT)' \
		bash ./test_driver_smoke.sh

run-integration: all
	MPIEXEC='$(MPIEXEC)' EXEC_DIR='$(abspath $(EXEC_DIR))' \
		DRIVER_INTEGRATION_TIMEOUT='$(DRIVER_INTEGRATION_TIMEOUT)' \
		$(PYTHON) ./test_driver_integration.py

run-performance-self-test:
	$(PYTHON) ./test_openmp_performance.py --self-test

validate-performance-paths:
	@set -eu; \
	performance_root=$$(realpath -m -- '$(PERFORMANCE_BUILD_ROOT_PATH)'); \
	repository_root=$$(realpath -m -- '$(REPOSITORY_ROOT)'); \
	build_root=$$(realpath -m -- '$(BUILD_ROOT_PATH)'); \
	object_root=$$(realpath -m -- '$(PERFORMANCE_OBJ_DIR)'); \
	case "$$performance_root" in \
		/|"$$repository_root"|"$$repository_root"/build) \
			printf 'Unsafe PERFORMANCE_BUILD_ROOT: %s\n' \
				'$(PERFORMANCE_BUILD_ROOT_PATH)' >&2; exit 2 ;; \
		"$$repository_root"/build/*) ;; \
		"$$repository_root"/*) \
			printf 'Unsafe PERFORMANCE_BUILD_ROOT: %s\n' \
				'$(PERFORMANCE_BUILD_ROOT_PATH)' >&2; exit 2 ;; \
	esac; \
	case "$$performance_root" in \
		"$$build_root"|"$$build_root"/*) \
			printf 'Unsafe PERFORMANCE_BUILD_ROOT: %s\n' \
				'$(PERFORMANCE_BUILD_ROOT_PATH)' >&2; exit 2 ;; \
	esac; \
	case "$$object_root" in \
		"$$performance_root"/*) ;; \
		*) printf 'PERFORMANCE_OBJ_DIR must be inside PERFORMANCE_BUILD_ROOT: %s\n' \
			'$(PERFORMANCE_OBJ_DIR)' >&2; exit 2 ;; \
	esac

prepare-performance-root: validate-performance-paths
	@set -eu; \
	root='$(PERFORMANCE_BUILD_ROOT_PATH)'; \
	marker='$(PERFORMANCE_MARKER)'; \
	actual=$$(realpath -m -- "$$root"); \
	expected="ADS_MPI_OPENMP_PERFORMANCE_ROOT=$$actual"; \
	if test -e "$$root" && ! test -d "$$root"; then \
		printf 'PERFORMANCE_BUILD_ROOT is not a directory: %s\n' "$$root" >&2; \
		exit 2; \
	fi; \
	if test -d "$$root"; then \
		if test -L "$$marker"; then \
			printf 'Unsafe performance ownership marker: %s\n' "$$marker" >&2; \
			exit 2; \
		elif test -f "$$marker"; then \
			if test "$$(cat -- "$$marker")" != "$$expected"; then \
				printf 'Invalid performance ownership marker: %s\n' "$$marker" >&2; \
				exit 2; \
			fi; \
		elif test -n "$$(find "$$root" -mindepth 1 -maxdepth 1 -print -quit)"; then \
			printf 'Refusing non-empty unowned PERFORMANCE_BUILD_ROOT: %s\n' \
				"$$root" >&2; exit 2; \
		fi; \
	fi; \
	mkdir -p -- "$$root"; \
	if ! test -f "$$marker"; then \
		marker_tmp="$$marker.tmp"; \
		printf '%s\n' "$$expected" > "$$marker_tmp"; \
		mv -f -- "$$marker_tmp" "$$marker"; \
	fi

validate-performance-ownership: validate-performance-paths
	@set -eu; \
	root='$(PERFORMANCE_BUILD_ROOT_PATH)'; \
	marker='$(PERFORMANCE_MARKER)'; \
	if ! test -e "$$root"; then exit 0; fi; \
	if ! test -d "$$root"; then \
		printf 'PERFORMANCE_BUILD_ROOT is not a directory: %s\n' "$$root" >&2; \
		exit 2; \
	fi; \
	actual=$$(realpath -m -- "$$root"); \
	expected="ADS_MPI_OPENMP_PERFORMANCE_ROOT=$$actual"; \
	if test -L "$$marker" || ! test -f "$$marker" || \
		test "$$(cat -- "$$marker" 2>/dev/null || :)" != "$$expected"; then \
		printf 'Refusing unowned PERFORMANCE_BUILD_ROOT: %s\n' "$$root" >&2; \
		exit 2; \
	fi

performance-build: prepare-performance-root
	+$(MAKE) --no-print-directory -j1 \
		-C '$(PROBLEMS_PATH)/$(igrm_l2_DIR)' \
		CONFIG='$(CONFIG_PATH)' BUILD=release \
		BUILD_ROOT='$(PERFORMANCE_BUILD_ROOT_PATH)' \
		PROBLEM_OBJ_DIR='$(PERFORMANCE_OBJ_DIR)' \
		COMPILER='$(MPIFC)' MPIFC='$(MPIFC)' build

run-performance: run-performance-self-test performance-build
	MPIEXEC='$(MPIEXEC)' EXECUTABLE='$(PERFORMANCE_EXECUTABLE)' \
		PERFORMANCE_TIMEOUT='$(PERFORMANCE_TIMEOUT)' \
		PERFORMANCE_WARMUPS='$(PERFORMANCE_WARMUPS)' \
		PERFORMANCE_SAMPLES='$(PERFORMANCE_SAMPLES)' \
		PERFORMANCE_MIN_SPEEDUP='$(PERFORMANCE_MIN_SPEEDUP)' \
		PERFORMANCE_MAX_REGRESSION='$(PERFORMANCE_MAX_REGRESSION)' \
		PERFORMANCE_BASELINE='$(PERFORMANCE_BASELINE)' \
		PERFORMANCE_OUTPUT='$(PERFORMANCE_RESULT)' \
		MPIEXEC_FLAGS='$(MPIEXEC_FLAGS)' MPI_NP_FLAG='$(MPI_NP_FLAG)' \
		OMP_PROC_BIND='$(OMP_PROC_BIND)' OMP_PLACES='$(OMP_PLACES)' \
		$(PYTHON) ./test_openmp_performance.py

run: run-cli run-smoke run-integration run-performance

clean-performance: validate-performance-ownership
	+$(MAKE) --no-print-directory -j1 \
		-C '$(PROBLEMS_PATH)/$(igrm_l2_DIR)' \
		CONFIG='$(CONFIG_PATH)' BUILD=release \
		BUILD_ROOT='$(PERFORMANCE_BUILD_ROOT_PATH)' \
		PROBLEM_OBJ_DIR='$(PERFORMANCE_OBJ_DIR)' \
		COMPILER='$(MPIFC)' MPIFC='$(MPIFC)' clean
	+$(MAKE) --no-print-directory -j1 -C '$(SRC_PATH)' \
		CONFIG='$(CONFIG_PATH)' BUILD=release \
		BUILD_ROOT='$(PERFORMANCE_BUILD_ROOT_PATH)' \
		COMPILER='$(MPIFC)' MPIFC='$(MPIFC)' clean
	$(RM) -- '$(PERFORMANCE_RESULT)'
	@if test -d '$(PERFORMANCE_BUILD_ROOT_PATH)' && \
		test -z "$$(find '$(PERFORMANCE_BUILD_ROOT_PATH)' -mindepth 1 \
			-maxdepth 1 ! -name '.ads-openmp-performance-root' -print -quit)"; then \
		$(RM) -- '$(PERFORMANCE_MARKER)'; \
		rmdir -- '$(PERFORMANCE_BUILD_ROOT_PATH)' 2>/dev/null || :; \
	fi

clean: clean-performance
	$(RM) -- $(ORACLE_OBJECT) $(ORACLE_MODULE) $(ORACLE_EXECUTABLE)
	@rmdir -- '$(ORACLE_OBJ_DIR)' 2>/dev/null || :
