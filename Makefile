.DEFAULT_GOAL := all

# Problem modules intentionally reuse the same Fortran module and object names.
# Keep every recursive build and test suite strictly sequential.
.NOTPARALLEL:

CONFIG ?= m_options
CONFIG_PATH := $(abspath $(CONFIG))
include $(CONFIG_PATH)

MYMAKE_DIR := mymake
TESTS_DIR := tests
EXEC_DIR := $(abspath $(MYMAKE_DIR)/EXEC)
DOXYFILE ?= doxygen.conf

PROBLEMS := l2 heat eriksson pure_diffusion_igrm oil igrm_l2 igrm_heat \
	igrm_eirksson igrm_stokes igrm_pollution
PROBLEM ?= l2

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

l2_ARGS ?= 2 2 2 1 1 1 1
heat_ARGS ?= 2 1 1 0.01 1 1 1
eriksson_ARGS ?= 2 1 1 0.01 1 1 1
pure_diffusion_igrm_ARGS ?= 3 1 1 1 1 1 0.1 dg
oil_ARGS ?= 2 1 1 1 1 1 0.1 1 0.5 0.5 0.5 1 0.25 0.25 0.25
igrm_l2_ARGS ?= 2 2 2 3 3 3 1 1 1 1 1 1 dg
igrm_heat_ARGS ?= 2 2 2 3 3 3 2 2 2 1 1 1 1 0.001 dg
igrm_eirksson_ARGS ?= 4 4 4 3 3 3 2 2 2 1 1 1
igrm_stokes_ARGS ?= 4 4 4 2 2 2 2 2 2 1 1 1
igrm_pollution_ARGS ?= 4 0 1 0 2 1 1 1 1 1

TEST_OPTIONS = \
	CONFIG="$(CONFIG_PATH)" \
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
	clean clean-build clean-problems clean-library clean-tests clean-docs \
	distclean

all: library problems

build-all: all

build: build-$(PROBLEM)

library:
	+$(MAKE) --no-print-directory -j1 -C $(MYMAKE_DIR) \
		CONFIG="$(CONFIG_PATH)" library

problems: $(addprefix build-,$(PROBLEMS))

define DEFINE_BUILD_TARGET
.PHONY: build-$(1) $(1)
build-$(1):
	+$$(MAKE) --no-print-directory -j1 -C $$(MYMAKE_DIR) \
		CONFIG="$$(CONFIG_PATH)" $(1)

$(1): build-$(1)
endef
$(foreach problem,$(PROBLEMS),$(eval $(call DEFINE_BUILD_TARGET,$(problem))))

rebuild: clean-build all

list-problems:
	@printf '%s\n' $(PROBLEMS)

list-configs:
	@printf '%s\n' 'm_options (active local configuration)' \
		$(sort $(wildcard makeconfig/*.mk))

define DEFINE_RUN_TARGET
.PHONY: run-$(1)
run-$(1): build-$(1)
	@run_dir='$$(if $$(strip $$(RUN_DIR)),$$(RUN_DIR),output/$(1))'; \
	mkdir -p "$$$$run_dir"; \
	printf 'Running $(1) in %s with %s MPI rank(s), %s OpenMP thread(s)\n' \
		"$$$$run_dir" '$$(NP)' '$$(OMP_NUM_THREADS)'; \
	cd "$$$$run_dir" && \
	env OMP_DYNAMIC='$$(OMP_DYNAMIC)' \
		OMP_NUM_THREADS='$$(OMP_NUM_THREADS)' \
		OMP_PROC_BIND='$$(OMP_PROC_BIND)' \
		$$(if $$(and $$(filter oil,$(1)),$$(strip $$(OIL_SEED))),ADS_OIL_RANDOM_SEED='$$(OIL_SEED)') \
		$$(RUN_ENV) \
		"$$(MPIEXEC)" $$(MPIEXEC_FLAGS) "$$(MPI_NP_FLAG)" "$$(NP)" \
		"$$(EXEC_DIR)/$(1)" \
		$$(if $$(strip $$(ARGS)),$$(ARGS),$$($(1)_ARGS))
endef
$(foreach problem,$(PROBLEMS),$(eval $(call DEFINE_RUN_TARGET,$(problem))))

run: run-$(PROBLEM)

run-help:
	@printf '%s\n' \
		'Problem argument syntax (the values in parentheses are the defaults):' \
		'' \
		'l2:' \
		'  <isizex> <isizey> <isizez> <order> <procx> <procy> <procz>' \
		'  ($(l2_ARGS))' \
		'heat, eriksson:' \
		'  <size> <order> <steps> <dt> <procx> <procy> <procz>' \
		'  heat:     ($(heat_ARGS))' \
		'  eriksson: ($(eriksson_ARGS))' \
		'pure_diffusion_igrm:' \
		'  <size> <order> <procx> <procy> <procz> <steps> <dt> [dg|pr|be]' \
		'  ($(pure_diffusion_igrm_ARGS))' \
		'oil:' \
		'  <size> <order> <procx> <procy> <procz> <steps> <dt>' \
		'  <npumps> (<pump_x> <pump_y> <pump_z>)...' \
		'  <ndrains> (<drain_x> <drain_y> <drain_z>)...' \
		'  ($(oil_ARGS))' \
		'igrm_l2:' \
		'  <nelem_x> <nelem_y> <nelem_z> <ptest_x> <ptest_y> <ptest_z>' \
		'  <ptrial_x> <ptrial_y> <ptrial_z> <procx> <procy> <procz>' \
		'  [[tau] dg|pr|be]' \
		'  ($(igrm_l2_ARGS))' \
		'igrm_heat:' \
		'  <nelem_x> <nelem_y> <nelem_z> <ptest_x> <ptest_y> <ptest_z>' \
		'  <ptrial_x> <ptrial_y> <ptrial_z> <procx> <procy> <procz>' \
		'  <steps> <dt> [dg|pr|be]' \
		'  ($(igrm_heat_ARGS))' \
		'igrm_eirksson:' \
		'  <nelem_x> <nelem_y> <nelem_z> <ptest_x> <ptest_y> <ptest_z>' \
		'  <ptrial_x> <ptrial_y> <ptrial_z> <procx> <procy> <procz>' \
		'  ($(igrm_eirksson_ARGS))' \
		'igrm_stokes:' \
		'  <nelem_x> <nelem_y> <nelem_z> <ptest_x> <ptest_y> <ptest_z>' \
		'  <ptrial_x> <ptrial_y> <ptrial_z> <procx> <procy> <procz>' \
		'  ($(igrm_stokes_ARGS))' \
		'igrm_pollution:' \
		'  <N> <adapt:0|1> <p_trial> <C_trial> <p_test> <C_test>' \
		'  <steps> <procx> <procy> <procz>' \
		'  ($(igrm_pollution_ARGS))' \
		'' \
		'Use: make run[-<problem>] ARGS="..." NP=<procx*procy*procz>' \
		'Run make help for launcher, OpenMP, environment, and output controls.'

show-run:
	@case ' $(PROBLEMS) ' in \
		*' $(PROBLEM) '*) ;; \
		*) printf 'Unknown PROBLEM: %s\nAvailable: %s\n' '$(PROBLEM)' '$(PROBLEMS)'; exit 2 ;; \
	esac
	@printf '%-20s %s\n' \
		'PROBLEM' '$(PROBLEM)' \
		'EXECUTABLE' '$(EXEC_DIR)/$(PROBLEM)' \
		'ARGS' '$(if $(strip $(ARGS)),$(ARGS),$($(PROBLEM)_ARGS))' \
		'NP' '$(NP)' \
		'MPIEXEC' '$(MPIEXEC)' \
		'MPIEXEC_FLAGS' '$(MPIEXEC_FLAGS)' \
		'MPI_NP_FLAG' '$(MPI_NP_FLAG)' \
		'OMP_NUM_THREADS' '$(OMP_NUM_THREADS)' \
		'OMP_DYNAMIC' '$(OMP_DYNAMIC)' \
		'OMP_PROC_BIND' '$(OMP_PROC_BIND)' \
		'RUN_DIR' '$(if $(strip $(RUN_DIR)),$(RUN_DIR),output/$(PROBLEM))' \
		'RUN_ENV' '$(RUN_ENV)' \
		'OIL_SEED' '$(OIL_SEED)'

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
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR)/driver_cli \
		$(TEST_OPTIONS) run-cli

test-smoke:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR)/driver_cli \
		$(TEST_OPTIONS) run-smoke

test-integration:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) \
		$(TEST_OPTIONS) run-integration

test-list:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) $(TEST_OPTIONS) list

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

clean-build:
	+$(MAKE) --no-print-directory -j1 -C $(MYMAKE_DIR) \
		CONFIG="$(CONFIG_PATH)" clean-all

clean-problems:
	+$(MAKE) --no-print-directory -j1 -C $(MYMAKE_DIR) \
		CONFIG="$(CONFIG_PATH)" clean-problems

clean-library:
	+$(MAKE) --no-print-directory -j1 -C $(MYMAKE_DIR) \
		CONFIG="$(CONFIG_PATH)" clean-library

clean-tests:
	+$(MAKE) --no-print-directory -j1 -C $(TESTS_DIR) $(TEST_OPTIONS) clean

clean-docs:
	$(RM) -r -- doxygen

clean: clean-build clean-tests clean-docs

# Keep the selected configuration file. distclean is intentionally only an
# alias for generated-artifact cleanup.
distclean: clean

show-config config:
	@printf '%-28s %s\n' \
		'CONFIG' '$(CONFIG_PATH)' \
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
		'ADS MPI unified make interface' \
		'' \
		'Build:' \
		'  make | make all | make build-all  static library + all ten problems' \
		'  make library                      build mymake/LIB/libads.a' \
		'  make problems                     build all problem executables' \
		'  make build PROBLEM=heat           build one selected problem' \
		'  make build-heat                    named build target (available for all problems)' \
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
		'  make test-src                      library unit/MPI suites' \
		'  make test-problems                 problem callback suites' \
		'  make test-driver                   CLI + smoke + numerical integration' \
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
		'  make clean-build | clean-problems | clean-library | clean-tests | clean-docs' \
		'  make show-config                   print effective configuration' \
		'  make list-problems | list-configs' \
		'  BUILD=debug|release and all paths/tools are configured in root m_options.' \
		'  Select an example with CONFIG=makeconfig/<name>.mk.'
