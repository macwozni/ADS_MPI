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
PERFORMANCE_BUILD_ROOT ?= $(ROOT_DIR)build/openmp-performance
PERFORMANCE_BUILD_ROOT_PATH := $(abspath $(PERFORMANCE_BUILD_ROOT))
PERFORMANCE_TIMEOUT ?= 300
PERFORMANCE_SUITE_TIMEOUT ?= 3600s
PERFORMANCE_WARMUPS ?= 1
PERFORMANCE_SAMPLES ?= 3
PERFORMANCE_MIN_SPEEDUP ?= 1.10
PERFORMANCE_MAX_REGRESSION ?= 1.15
PERFORMANCE_BASELINE ?=
PERFORMANCE_BASELINE_PATH := $(if $(strip $(PERFORMANCE_BASELINE)),$(abspath $(PERFORMANCE_BASELINE)),)
COVERAGE_ROOT ?= $(ROOT_DIR)build/coverage
COVERAGE_ROOT_PATH := $(abspath $(COVERAGE_ROOT))
COVERAGE_BUILD_ROOT ?= $(COVERAGE_ROOT_PATH)/build
COVERAGE_BUILD_ROOT_PATH := $(abspath $(COVERAGE_BUILD_ROOT))
COVERAGE_TRACEFILE := $(COVERAGE_ROOT_PATH)/coverage.info
COVERAGE_SUMMARY := $(COVERAGE_ROOT_PATH)/coverage-summary.json
COVERAGE_HTML_DIR := $(COVERAGE_ROOT_PATH)/html
COVERAGE_MARKER := $(COVERAGE_ROOT_PATH)/.ads-coverage-root
COVERAGE_FLAGS ?= -O0 -g --coverage -fprofile-abs-path
COVERAGE_MIN_LINES ?= 90.0
COVERAGE_MIN_FUNCTIONS ?= 90.0
COVERAGE_MIN_BRANCHES ?= 50.0
GCOV ?= gcov
LCOV ?= lcov
GENHTML ?= genhtml
include $(PROBLEMS_DIR)/problems.mk

PROBLEM ?= l2
DOXYFILE ?= doxygen.conf

NP ?= 1
MPIEXEC_FLAGS ?=
MPI_NP_FLAG ?= -n
OMP_NUM_THREADS ?= 1
OMP_DYNAMIC ?= FALSE
OMP_PROC_BIND ?= close
OMP_PLACES ?= cores
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
	OMP_PLACES="$(OMP_PLACES)"

.PHONY: all build build-all library problems list-problems list-configs help targets \
	show-config config rebuild run run-help show-run \
	test check test-build test-layout test-src test-problems test-driver \
	test-build-system test-coverage _test-coverage-report validate-coverage-tools \
	validate-coverage-paths prepare-coverage-root validate-coverage-ownership \
	test-cli test-smoke test-integration test-performance \
	test-performance-self-test test-list test-suite \
	docs doc docs-html docs-pdf docs-check \
	clean clean-build clean-problems clean-library clean-legacy clean-tests clean-docs \
	clean-performance clean-coverage distclean

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

test-build-system:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) \
		$(TEST_OPTIONS) run-build-system

test-cli:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) $(TEST_OPTIONS) run-cli

test-smoke:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) $(TEST_OPTIONS) run-smoke

test-integration:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) $(TEST_OPTIONS) run-integration

test-performance:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) \
		$(TEST_OPTIONS) run-performance

test-performance-self-test:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) \
		$(TEST_OPTIONS) run-performance-self-test

validate-coverage-tools:
	@set -eu; \
	for tool in '$(GCOV)' '$(LCOV)' '$(GENHTML)'; do \
		if ! command -v "$$tool" >/dev/null 2>&1; then \
			printf 'Coverage tool not found: %s\n' "$$tool" >&2; \
			exit 2; \
		fi; \
	done; \
	compiler_version="$$( '$(MPIFC)' --version 2>&1 )"; \
	case "$$compiler_version" in \
		*'GNU Fortran'*) ;; \
		*) printf 'test-coverage requires GNU Fortran; MPIFC is %s\n' \
			'$(MPIFC)' >&2; exit 2 ;; \
	esac

validate-coverage-paths:
	@set -eu; \
	coverage_root=$$(realpath -m -- '$(COVERAGE_ROOT_PATH)'); \
	coverage_build=$$(realpath -m -- '$(COVERAGE_BUILD_ROOT_PATH)'); \
	repository_root=$$(realpath -m -- '$(ROOT_DIR)'); \
	build_root=$$(realpath -m -- '$(BUILD_ROOT_PATH)'); \
	performance_root=$$(realpath -m -- '$(PERFORMANCE_BUILD_ROOT_PATH)'); \
	case "$$coverage_root" in \
		"$$repository_root"/build/*) ;; \
		*) printf 'Unsafe COVERAGE_ROOT: %s\n' '$(COVERAGE_ROOT_PATH)' >&2; \
			exit 2 ;; \
	esac; \
	case "$$coverage_build" in \
		"$$coverage_root"/*) ;; \
		*) printf 'COVERAGE_BUILD_ROOT must be inside COVERAGE_ROOT: %s\n' \
			'$(COVERAGE_BUILD_ROOT_PATH)' >&2; exit 2 ;; \
	esac; \
	case "$$coverage_root" in \
		"$$build_root"|"$$build_root"/*|"$$performance_root"|"$$performance_root"/*) \
			printf 'COVERAGE_ROOT overlaps another build root: %s\n' \
				'$(COVERAGE_ROOT_PATH)' >&2; exit 2 ;; \
	esac; \
	case "$$build_root $$performance_root" in \
		*"$$coverage_root"/*) \
			printf 'COVERAGE_ROOT contains another build root: %s\n' \
				'$(COVERAGE_ROOT_PATH)' >&2; exit 2 ;; \
	esac; \
	if test -L '$(COVERAGE_ROOT_PATH)'; then \
		printf 'Unsafe COVERAGE_ROOT symlink: %s\n' '$(COVERAGE_ROOT_PATH)' >&2; \
		exit 2; \
	fi

prepare-coverage-root: validate-coverage-paths
	@set -eu; \
	root='$(COVERAGE_ROOT_PATH)'; \
	marker='$(COVERAGE_MARKER)'; \
	actual=$$(realpath -m -- "$$root"); \
	expected="ADS_MPI_COVERAGE_ROOT=$$actual"; \
	if test -e "$$root" && ! test -d "$$root"; then \
		printf 'COVERAGE_ROOT is not a directory: %s\n' "$$root" >&2; exit 2; \
	fi; \
	if test -d "$$root"; then \
		if test -L "$$marker"; then \
			printf 'Unsafe coverage ownership marker: %s\n' "$$marker" >&2; \
			exit 2; \
		elif test -f "$$marker"; then \
			if test "$$(cat -- "$$marker")" != "$$expected"; then \
				printf 'Invalid coverage ownership marker: %s\n' "$$marker" >&2; \
				exit 2; \
			fi; \
		elif test -n "$$(find "$$root" -mindepth 1 -maxdepth 1 -print -quit)"; then \
			printf 'Refusing non-empty unowned COVERAGE_ROOT: %s\n' "$$root" >&2; \
			exit 2; \
		fi; \
	fi; \
	$(RM) -r -- '$(COVERAGE_BUILD_ROOT_PATH)' '$(COVERAGE_HTML_DIR)'; \
	$(RM) -- '$(COVERAGE_TRACEFILE)' '$(COVERAGE_SUMMARY)' \
		'$(COVERAGE_SUMMARY).tmp'; \
	mkdir -p -- "$$root"; \
	if ! test -f "$$marker"; then \
		marker_tmp="$$marker.tmp"; \
		printf '%s\n' "$$expected" > "$$marker_tmp"; \
		mv -f -- "$$marker_tmp" "$$marker"; \
	fi

validate-coverage-ownership: validate-coverage-paths
	@set -eu; \
	root='$(COVERAGE_ROOT_PATH)'; \
	marker='$(COVERAGE_MARKER)'; \
	if ! test -e "$$root"; then exit 0; fi; \
	if ! test -d "$$root"; then \
		printf 'COVERAGE_ROOT is not a directory: %s\n' "$$root" >&2; exit 2; \
	fi; \
	actual=$$(realpath -m -- "$$root"); \
	expected="ADS_MPI_COVERAGE_ROOT=$$actual"; \
	if test -L "$$marker" || ! test -f "$$marker" || \
		test "$$(cat -- "$$marker" 2>/dev/null || :)" != "$$expected"; then \
		printf 'Refusing unowned COVERAGE_ROOT: %s\n' "$$root" >&2; exit 2; \
	fi

test-coverage:
	+@status=0; \
	$(MAKE) --no-print-directory -j1 _test-coverage-report || status=$$?; \
	cleanup_status=0; \
	$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) \
		$(TEST_OPTIONS) clean-coverage-data || cleanup_status=$$?; \
	if test "$$status" -ne 0; then exit "$$status"; fi; \
	if test "$$cleanup_status" -ne 0; then exit "$$cleanup_status"; fi; \
	printf 'LCOV tracefile: %s\nHTML report:    %s\nJSON summary:   %s\n' \
		'$(COVERAGE_TRACEFILE)' '$(COVERAGE_HTML_DIR)/index.html' \
		'$(COVERAGE_SUMMARY)'

_test-coverage-report: validate-coverage-tools prepare-coverage-root
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) \
		$(TEST_OPTIONS) BUILD_ROOT='$(COVERAGE_BUILD_ROOT_PATH)' \
		COVERAGE_FLAGS='$(COVERAGE_FLAGS)' run-coverage
	$(LCOV) --capture --directory '$(TESTS_DIR)' \
		--directory '$(COVERAGE_BUILD_ROOT_PATH)' \
		--base-directory '$(ROOT_DIR)' --gcov-tool '$(GCOV)' \
		--branch-coverage --function-coverage --filter branch \
		--rc geninfo_unexecuted_blocks=1 --ignore-errors inconsistent \
		--include '$(SRC_DIR)/*.F90' --output-file '$(COVERAGE_TRACEFILE)'
	$(GENHTML) --branch-coverage --function-coverage \
		--ignore-errors inconsistent \
		--title 'ADS MPI core coverage' --output-directory '$(COVERAGE_HTML_DIR)' \
		'$(COVERAGE_TRACEFILE)'
	$(PYTHON) '$(TESTS_DIR)/coverage_report/check_coverage.py' \
		--tracefile '$(COVERAGE_TRACEFILE)' \
		--source-manifest '$(SRC_DIR)/sources.mk' --source-root '$(SRC_DIR)' \
		--repository-root '$(ROOT_DIR)' \
		--object-root '$(COVERAGE_BUILD_ROOT_PATH)/_OBJ' --gcov '$(GCOV)' \
		--summary-json '$(COVERAGE_SUMMARY)' \
		--min-lines '$(COVERAGE_MIN_LINES)' \
		--min-functions '$(COVERAGE_MIN_FUNCTIONS)' \
		--min-branches '$(COVERAGE_MIN_BRANCHES)'

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

clean-performance:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) \
		$(TEST_OPTIONS) clean-performance

clean-coverage: validate-coverage-ownership
	@set -eu; \
	root='$(COVERAGE_ROOT_PATH)'; \
	marker='$(COVERAGE_MARKER)'; \
	if ! test -d "$$root"; then exit 0; fi; \
	$(RM) -r -- '$(COVERAGE_BUILD_ROOT_PATH)' '$(COVERAGE_HTML_DIR)'; \
	$(RM) -- '$(COVERAGE_TRACEFILE)' '$(COVERAGE_SUMMARY)' \
		'$(COVERAGE_SUMMARY).tmp'; \
	if test -z "$$(find "$$root" -mindepth 1 -maxdepth 1 \
		! -name '.ads-coverage-root' -print -quit)"; then \
		$(RM) -- "$$marker"; \
		rmdir -- "$$root" 2>/dev/null || :; \
	fi

clean-docs:
	$(RM) -r -- doxygen

clean: clean-build clean-tests clean-docs clean-coverage

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
		'MPIEXEC_FLAGS' '$(MPIEXEC_FLAGS)' \
		'MPI_NP_FLAG' '$(MPI_NP_FLAG)' \
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
		'SKIP_MPI_CASES' '$(SKIP_MPI_CASES)' \
		'PERFORMANCE_BUILD_ROOT' '$(PERFORMANCE_BUILD_ROOT_PATH)' \
		'PERFORMANCE_TIMEOUT' '$(PERFORMANCE_TIMEOUT)' \
		'PERFORMANCE_SUITE_TIMEOUT' '$(PERFORMANCE_SUITE_TIMEOUT)' \
		'PERFORMANCE_WARMUPS' '$(PERFORMANCE_WARMUPS)' \
		'PERFORMANCE_SAMPLES' '$(PERFORMANCE_SAMPLES)' \
		'PERFORMANCE_MIN_SPEEDUP' '$(PERFORMANCE_MIN_SPEEDUP)' \
		'PERFORMANCE_MAX_REGRESSION' '$(PERFORMANCE_MAX_REGRESSION)' \
		'PERFORMANCE_BASELINE' '$(PERFORMANCE_BASELINE_PATH)' \
		'COVERAGE_ROOT' '$(COVERAGE_ROOT_PATH)' \
		'COVERAGE_BUILD_ROOT' '$(COVERAGE_BUILD_ROOT_PATH)' \
		'COVERAGE_FLAGS' '$(COVERAGE_FLAGS)' \
		'COVERAGE_MIN_LINES' '$(COVERAGE_MIN_LINES)' \
		'COVERAGE_MIN_FUNCTIONS' '$(COVERAGE_MIN_FUNCTIONS)' \
		'COVERAGE_MIN_BRANCHES' '$(COVERAGE_MIN_BRANCHES)' \
		'GCOV' '$(GCOV)' \
		'LCOV' '$(LCOV)' \
		'GENHTML' '$(GENHTML)' \
		'OMP_PROC_BIND' '$(OMP_PROC_BIND)' \
		'OMP_PLACES' '$(OMP_PLACES)'

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
		'  make test | make check             complete run, including performance gate' \
		'  make test-build                    build all tests without running' \
		'  make test-layout                   verify one-to-one src/test layout' \
		'  make test-src | test-problems | test-driver' \
		'  make test-build-system             test the hierarchical Make interface' \
		'  make test-cli | test-smoke | test-integration' \
		'  make test-performance              release OMP1/OMP4 scaling gate' \
		'  make test-performance-self-test    test timing-gate logic without MPI' \
		'  make test-coverage                 GNU gcov/lcov report + coverage gates' \
		'    COVERAGE_MIN_LINES/FUNCTIONS/BRANCHES override the three gates' \
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
		'  make clean-tests | clean-performance | clean-coverage | clean-docs' \
		'  make show-config                   print effective configuration' \
		'  make list-problems | list-configs' \
		'  BUILD=debug|release and all paths/tools are configured in root m_options.' \
		'  Select an example with CONFIG=makeconfig/<name>.mk.'
