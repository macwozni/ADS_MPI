.DEFAULT_GOAL := all

# Public repository entry point.  Compilation belongs to src/, problem builds
# to problems/, and test registration/execution to tests/.
.NOTPARALLEL:

ROOT_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
SRC_DIR := $(ROOT_DIR)src
PROBLEMS_DIR := $(ROOT_DIR)problems
TESTS_DIR := $(ROOT_DIR)tests
MYMAKE_DIR := $(ROOT_DIR)mymake

CONFIG ?= m_options
CONFIG_PATH := $(abspath $(CONFIG))
include $(CONFIG_PATH)
BUILD_ROOT ?= mymake
BUILD_ROOT_PATH := $(abspath $(BUILD_ROOT))
include $(PROBLEMS_DIR)/problems.mk

PROBLEM ?= l2
DOXYFILE ?= doxygen.conf

NP ?= 1
MPIEXEC_FLAGS ?=
MPI_NP_FLAG ?= -n
OMP_NUM_THREADS ?= 1
OMP_DYNAMIC ?= FALSE
OMP_PROC_BIND ?= close
RUN_ENV ?=
RUN_DIR ?=
ARGS ?=
OIL_SEED ?=
TEST_SUITE ?=

PROBLEM_OPTIONS = \
	CONFIG="$(CONFIG_PATH)" \
	BUILD_ROOT="$(BUILD_ROOT_PATH)" \
	ARGS="$(ARGS)" \
	NP="$(NP)" \
	MPIEXEC="$(MPIEXEC)" \
	MPIEXEC_FLAGS="$(MPIEXEC_FLAGS)" \
	MPI_NP_FLAG="$(MPI_NP_FLAG)" \
	OMP_NUM_THREADS="$(OMP_NUM_THREADS)" \
	OMP_DYNAMIC="$(OMP_DYNAMIC)" \
	OMP_PROC_BIND="$(OMP_PROC_BIND)" \
	RUN_ENV="$(RUN_ENV)" \
	RUN_DIR="$(RUN_DIR)" \
	OIL_SEED="$(OIL_SEED)"

TEST_OPTIONS = \
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
	SKIP_MPI_CASES="$(SKIP_MPI_CASES)"

.PHONY: all build build-all library problems list-problems list-configs help targets \
	show-config config rebuild run run-help show-run \
	test check test-build test-layout test-src test-problems test-driver \
	test-cli test-smoke test-integration test-list test-suite \
	docs doc docs-html docs-pdf docs-check \
	clean clean-build clean-problems clean-library clean-legacy clean-tests clean-docs \
	distclean

all: library problems

build-all: all

library:
	+$(MAKE) --no-print-directory -j1 -C $(SRC_DIR) \
		CONFIG="$(CONFIG_PATH)" BUILD_ROOT="$(BUILD_ROOT_PATH)" library

problems:
	+$(MAKE) --no-print-directory -j1 -C $(PROBLEMS_DIR) \
		$(PROBLEM_OPTIONS) problems

build:
	+$(MAKE) --no-print-directory -j1 -C $(PROBLEMS_DIR) \
		$(PROBLEM_OPTIONS) PROBLEM="$(PROBLEM)" build

define DEFINE_PROBLEM_TARGETS
.PHONY: build-$(1) $(1) run-$(1)
build-$(1) $(1):
	+$$(MAKE) --no-print-directory -j1 -C $$(PROBLEMS_DIR) \
		$$(PROBLEM_OPTIONS) build-$(1)

run-$(1):
	+$$(MAKE) --no-print-directory -j1 -C $$(PROBLEMS_DIR) \
		$$(PROBLEM_OPTIONS) run-$(1)
endef
$(foreach problem,$(PROBLEMS),$(eval $(call DEFINE_PROBLEM_TARGETS,$(problem))))

rebuild: clean-build all

list-problems:
	+@$(MAKE) --no-print-directory -s -C $(PROBLEMS_DIR) list-problems

list-configs:
	@printf '%s\n' 'm_options (active local configuration)' \
		$(sort $(wildcard makeconfig/*.mk))

run:
	+$(MAKE) --no-print-directory -j1 -C $(PROBLEMS_DIR) \
		$(PROBLEM_OPTIONS) PROBLEM="$(PROBLEM)" run

run-help:
	+@$(MAKE) --no-print-directory -j1 -C $(PROBLEMS_DIR) \
		CONFIG="$(CONFIG_PATH)" BUILD_ROOT="$(BUILD_ROOT_PATH)" run-help

show-run:
	+@$(MAKE) --no-print-directory -j1 -C $(PROBLEMS_DIR) \
		$(PROBLEM_OPTIONS) PROBLEM="$(PROBLEM)" show-run

test check:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) $(TEST_OPTIONS) run

test-build:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) $(TEST_OPTIONS) clean
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) $(TEST_OPTIONS) all

test-layout:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) $(TEST_OPTIONS) check-layout

test-src:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) $(TEST_OPTIONS) run-src

test-problems:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) $(TEST_OPTIONS) run-problems

test-driver:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) $(TEST_OPTIONS) run-driver

test-cli:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) $(TEST_OPTIONS) run-cli

test-smoke:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) $(TEST_OPTIONS) run-smoke

test-integration:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) $(TEST_OPTIONS) run-integration

test-list:
	+@$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) $(TEST_OPTIONS) list

test-suite:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) \
		$(TEST_OPTIONS) TEST_SUITE="$(TEST_SUITE)" run-suite

docs doc: docs-pdf

docs-html:
	$(DOXYGEN) $(DOXYFILE)

docs-pdf: docs-html
	+$(MAKE) --no-print-directory -C doxygen/latex pdf

docs-check:
	$(DOXYGEN) -x $(DOXYFILE) >/dev/null

clean-build: clean-problems clean-legacy
	+$(MAKE) --no-print-directory -j1 -C $(SRC_DIR) \
		CONFIG="$(CONFIG_PATH)" BUILD_ROOT="$(BUILD_ROOT_PATH)" clean

clean-problems:
	+$(MAKE) --no-print-directory -j1 -C $(PROBLEMS_DIR) \
		CONFIG="$(CONFIG_PATH)" BUILD_ROOT="$(BUILD_ROOT_PATH)" clean

clean-library:
	+$(MAKE) --no-print-directory -j1 -C $(SRC_DIR) \
		CONFIG="$(CONFIG_PATH)" BUILD_ROOT="$(BUILD_ROOT_PATH)" clean-library

clean-legacy:
	+$(MAKE) --no-print-directory -j1 -C $(MYMAKE_DIR) \
		CONFIG="$(CONFIG_PATH)" BUILD_ROOT="$(BUILD_ROOT_PATH)" legacy-clean-all

clean-tests:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) $(TEST_OPTIONS) clean

clean-docs:
	$(RM) -r -- doxygen

clean: clean-build clean-tests clean-docs

# The selected configuration is user-owned and survives both cleanup targets.
distclean: clean

show-config config:
	@printf '%-28s %s\n' \
		'CONFIG' '$(CONFIG_PATH)' \
		'BUILD_ROOT' '$(BUILD_ROOT_PATH)' \
		'BUILD' '$(BUILD)' \
		'COMPILER' '$(COMPILER)' \
		'MPIFC' '$(MPIFC)' \
		'MODULE_OUTPUT' '$(MODULE_OUTPUT)' \
		'MPIEXEC' '$(MPIEXEC)' \
		'MUMPS_DIR' '$(MUMPS_DIR)' \
		'LAPACK_DIR' '$(LAPACK_DIR)' \
		'SCALAPACK_DIR' '$(SCALAPACK_DIR)' \
		'BLAS_DIR' '$(BLAS_DIR)' \
		'PARMETIS_DIR' '$(PARMETIS_DIR)' \
		'METIS_DIR' '$(METIS_DIR)' \
		'GKLIB_DIR' '$(GKLIB_DIR)' \
		'PFUNIT_ROOT' '$(PFUNIT_ROOT)' \
		'DOXYGEN' '$(DOXYGEN)' \
		'PYTHON' '$(PYTHON)' \
		'SUITE_TIMEOUT' '$(SUITE_TIMEOUT)' \
		'DRIVER_CLI_TIMEOUT' '$(DRIVER_CLI_TIMEOUT)' \
		'DRIVER_SMOKE_TIMEOUT' '$(DRIVER_SMOKE_TIMEOUT)' \
		'DRIVER_INTEGRATION_TIMEOUT' '$(DRIVER_INTEGRATION_TIMEOUT)' \
		'SKIP_MPI_CASES' '$(SKIP_MPI_CASES)'

targets help:
	@printf '%s\n' \
		'ADS MPI hierarchical make interface' \
		'' \
		'Build:' \
		'  make | make all | make build-all  static library + all ten problems' \
		'  make library                      delegate the core build to src/' \
		'  make problems                     delegate all drivers to problems/' \
		'  make build PROBLEM=heat           build one selected problem' \
		'  make build-heat | make heat        named problem build target' \
		'  make rebuild                      clean build artifacts and build everything' \
		'' \
		'Run:' \
		"  make run PROBLEM=heat [ARGS='...'] [NP=1]" \
		"  make run-heat [ARGS='...'] [NP=1]" \
		'  Named run targets exist for: $(PROBLEMS)' \
		'  make run-help                     print every problem argument syntax' \
		'  make show-run [PROBLEM=heat]       print the effective run settings' \
		'  Variables: ARGS NP MPIEXEC MPIEXEC_FLAGS MPI_NP_FLAG RUN_DIR RUN_ENV' \
		'             OMP_NUM_THREADS OMP_DYNAMIC OMP_PROC_BIND OIL_SEED' \
		'  NP must equal procx*procy*procz contained in ARGS.' \
		'' \
		'Tests:' \
		'  make test | make check             complete test run' \
		'  make test-build                    build all tests without running' \
		'  make test-layout                   verify one-to-one src/test layout' \
		'  make test-src | test-problems | test-driver' \
		'  make test-cli | test-smoke | test-integration' \
		'  make test-suite TEST_SUITE=rhs_assembly' \
		'  make test-list' \
		'' \
		'Documentation:' \
		'  make docs | make doc               HTML and PDF' \
		'  make docs-html | docs-pdf | docs-check' \
		'' \
		'Cleanup/configuration:' \
		'  make clean                         build + tests + generated documentation' \
		'  make clean-build | clean-problems | clean-library | clean-legacy' \
		'  make clean-tests | clean-docs' \
		'  make show-config                   print effective configuration' \
		'  make list-problems | list-configs' \
		'  BUILD=debug|release and all paths/tools are configured in root m_options.' \
		'  Select an example with CONFIG=makeconfig/<name>.mk.'
