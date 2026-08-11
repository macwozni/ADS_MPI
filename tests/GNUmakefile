.DEFAULT_GOAL := all

.NOTPARALLEL:

.PHONY: all tests run run-src run-problems run-driver run-integration check-layout list clean

PFUNIT_ROOT ?= /opt/lib/pfunit/PFUNIT-4.16
MPIEXEC     ?= /opt/lib/mpich-5.0.0/bin/mpiexec
MPIFC       ?= /opt/lib/mpich-5.0.0/bin/mpif90
MUMPS_DIR   ?= /opt/lib/MUMPS_5.8.2
SUITE_TIMEOUT ?= 600s

# Exactly one primary suite for every active source listed in mymake/m_files.
SRC_SUITES := \
	setup \
	interfaces \
	argument_parser \
	knot_vector \
	parallelism \
	communicators \
	gauss \
	basis \
	math \
	reorder_rhs \
	my_mpi \
	sparse \
	rhs_eq \
	igrm_space \
	operator_assembly \
	rhs_assembly \
	solution_reconstruction \
	projection_engine \
	ads_lifecycle \
	mumps_solver \
	ads_directional_solve \
	utils \
	plot \
	gnuplot \
	vtk \
	solution_output \
	ads \
	time_scheme

# Problem modules deliberately live in separate build directories because
# several problems use the same Fortran module names (input_data/RHS_fun).
PROBLEM_SUITES := \
	oil_input_data \
	oil_rhs_fun \
	heat_rhs_fun \
	igrm_heat_rhs_fun

DRIVER_SUITES := driver_cli
SUITES := $(SRC_SUITES) $(PROBLEM_SUITES) $(DRIVER_SUITES)

define run_suites
	@set -eu; \
	for suite in $(1); do \
		printf '\n==> tests/%s: %s\n' "$$suite" "$(2)"; \
		timeout $(SUITE_TIMEOUT) $(MAKE) -j1 -C "$$suite" \
			PFUNIT_ROOT="$(PFUNIT_ROOT)" MPIEXEC="$(MPIEXEC)" \
			MPIFC="$(MPIFC)" FC="$(MPIFC)" \
			MUMPS_DIR="$(MUMPS_DIR)" "$(2)"; \
	done
endef

define rebuild_and_run_suites
	@set -eu; \
	for suite in $(1); do \
		printf '\n==> tests/%s: clean + run\n' "$$suite"; \
		timeout $(SUITE_TIMEOUT) $(MAKE) -j1 -C "$$suite" \
			PFUNIT_ROOT="$(PFUNIT_ROOT)" MPIEXEC="$(MPIEXEC)" \
			MPIFC="$(MPIFC)" FC="$(MPIFC)" \
			MUMPS_DIR="$(MUMPS_DIR)" clean; \
		timeout $(SUITE_TIMEOUT) $(MAKE) -j1 -C "$$suite" \
			PFUNIT_ROOT="$(PFUNIT_ROOT)" MPIEXEC="$(MPIEXEC)" \
			MPIFC="$(MPIFC)" FC="$(MPIFC)" \
			MUMPS_DIR="$(MUMPS_DIR)" run; \
	done
endef

check-layout:
	python3 ./check_test_layout.py

all: check-layout
	$(call run_suites,$(SUITES),all)

tests: all

run-src: check-layout
	$(call rebuild_and_run_suites,$(SRC_SUITES))

run-problems: check-layout
	$(call rebuild_and_run_suites,$(PROBLEM_SUITES))

run-driver:
	$(call rebuild_and_run_suites,$(DRIVER_SUITES))

run-integration: check-layout
	$(call run_suites,$(DRIVER_SUITES),run-integration)

run: run-src run-problems run-driver

list:
	@printf '%s\n' $(SUITES)

clean:
	$(call run_suites,$(SUITES),clean)
