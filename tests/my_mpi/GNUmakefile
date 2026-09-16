.DEFAULT_GOAL := all

.PHONY: all tests run run-positive fault clean

SRC_DIR   ?= ../../src
MPIEXEC   ?= /opt/lib/mpich-5.0.0/bin/mpiexec
FC        := /opt/lib/mpich-5.0.0/bin/mpif90
OBJ_DIR   := _OBJ
TEST_NAME := test_my_mpi
FAULT_OBJ_DIR := _OBJ_fault
FAULT_TEST_NAME := distribute_failure_probe
FAULT_CONTRACT_TEST_NAME := distribute_error_contract_probe
FAULT_CONTRACT_OBJ := $(FAULT_OBJ_DIR)/distribute_error_contract_probe.o
FAULT_FC := /usr/bin/gfortran
FAULT_TIMEOUT ?= 5s

FFLAGS += -std=f2003 -fimplicit-none -pedantic -O0 -g
FFLAGS += -fcheck=all -fbounds-check -fopenmp
FFLAGS += -J$(OBJ_DIR) -I$(OBJ_DIR) -I$(SRC_DIR)

OBJS := \
	$(OBJ_DIR)/Setup.o \
	$(OBJ_DIR)/knot_vector.o \
	$(OBJ_DIR)/parallelism.o \
	$(OBJ_DIR)/communicators.o \
	$(OBJ_DIR)/gauss.o \
	$(OBJ_DIR)/basis.o \
	$(OBJ_DIR)/my_mpi.o \
	$(OBJ_DIR)/ads_lifecycle.o \
	$(OBJ_DIR)/test_my_mpi.o

FAULT_FFLAGS := -std=f2003 -fimplicit-none -pedantic -O0 -g
FAULT_FFLAGS += -fcheck=all -fbounds-check
FAULT_FFLAGS += -cpp
FAULT_FFLAGS += -J$(FAULT_OBJ_DIR) -I$(FAULT_OBJ_DIR) -I$(SRC_DIR)
FAULT_FFLAGS += $(COVERAGE_FLAGS)

FAULT_OBJS := \
	$(FAULT_OBJ_DIR)/distribute_mpi_stub.o \
	$(FAULT_OBJ_DIR)/Setup.o \
	$(FAULT_OBJ_DIR)/distribute_dependency_stubs.o \
	$(FAULT_OBJ_DIR)/my_mpi.o \
	$(FAULT_OBJ_DIR)/distribute_failure_probe.o

$(OBJ_DIR):
	mkdir -p $@

$(FAULT_OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/Setup.o: $(SRC_DIR)/Setup.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/knot_vector.o: $(SRC_DIR)/knot_vector.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/parallelism.o: $(SRC_DIR)/parallelism.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/communicators.o: $(SRC_DIR)/communicators.F90 $(OBJ_DIR)/parallelism.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/gauss.o: $(SRC_DIR)/gauss.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/basis.o: $(SRC_DIR)/basis.F90 $(OBJ_DIR)/gauss.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/my_mpi.o: $(SRC_DIR)/my_mpi.F90 \
		$(OBJ_DIR)/Setup.o $(OBJ_DIR)/parallelism.o $(OBJ_DIR)/communicators.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/ads_lifecycle.o: $(SRC_DIR)/ads_lifecycle.F90 \
		$(OBJ_DIR)/Setup.o $(OBJ_DIR)/knot_vector.o $(OBJ_DIR)/parallelism.o \
		$(OBJ_DIR)/basis.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/test_my_mpi.o: test_my_mpi.F90 \
		$(OBJ_DIR)/Setup.o $(OBJ_DIR)/parallelism.o $(OBJ_DIR)/communicators.o \
		$(OBJ_DIR)/my_mpi.o $(OBJ_DIR)/ads_lifecycle.o
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_NAME): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

$(FAULT_OBJ_DIR)/distribute_mpi_stub.o: distribute_mpi_stub.F90 | $(FAULT_OBJ_DIR)
	$(FAULT_FC) $(FAULT_FFLAGS) -c $< -o $@

$(FAULT_OBJ_DIR)/Setup.o: $(SRC_DIR)/Setup.F90 | $(FAULT_OBJ_DIR)
	$(FAULT_FC) $(FAULT_FFLAGS) -c $< -o $@

$(FAULT_OBJ_DIR)/distribute_dependency_stubs.o: distribute_dependency_stubs.F90 | $(FAULT_OBJ_DIR)
	$(FAULT_FC) $(FAULT_FFLAGS) -c $< -o $@

$(FAULT_OBJ_DIR)/my_mpi.o: $(SRC_DIR)/my_mpi.F90 \
		$(FAULT_OBJ_DIR)/distribute_mpi_stub.o \
		$(FAULT_OBJ_DIR)/Setup.o \
		$(FAULT_OBJ_DIR)/distribute_dependency_stubs.o
	$(FAULT_FC) $(FAULT_FFLAGS) -c $< -o $@

$(FAULT_OBJ_DIR)/distribute_failure_probe.o: distribute_failure_probe.F90 \
		$(FAULT_OBJ_DIR)/distribute_mpi_stub.o \
		$(FAULT_OBJ_DIR)/Setup.o \
		$(FAULT_OBJ_DIR)/my_mpi.o
	$(FAULT_FC) $(FAULT_FFLAGS) -c $< -o $@

$(FAULT_CONTRACT_OBJ): distribute_failure_probe.F90 \
		$(FAULT_OBJ_DIR)/distribute_mpi_stub.o \
		$(FAULT_OBJ_DIR)/Setup.o \
		$(FAULT_OBJ_DIR)/my_mpi.o
	$(FAULT_FC) $(FAULT_FFLAGS) -DDISTRIBUTE_ERROR_STATUS_API -c $< -o $@

$(FAULT_TEST_NAME): $(FAULT_OBJS)
	$(FAULT_FC) $(FAULT_FFLAGS) -o $@ $(FAULT_OBJS)

$(FAULT_CONTRACT_TEST_NAME): \
		$(filter-out $(FAULT_OBJ_DIR)/distribute_failure_probe.o,$(FAULT_OBJS)) \
		$(FAULT_CONTRACT_OBJ)
	$(FAULT_FC) $(FAULT_FFLAGS) -o $@ $^

all: $(TEST_NAME) $(FAULT_TEST_NAME)

tests: all

run-positive: $(TEST_NAME)
	$(MPIEXEC) -n 1 ./$(TEST_NAME) 1 1 1
	$(MPIEXEC) -n 8 ./$(TEST_NAME) 2 2 2
	$(MPIEXEC) -n 4 ./$(TEST_NAME) 4 1 1
	$(MPIEXEC) -n 4 ./$(TEST_NAME) 1 4 1
	$(MPIEXEC) -n 4 ./$(TEST_NAME) 1 1 4
	$(MPIEXEC) -n 8 ./$(TEST_NAME) 4 2 1
	$(MPIEXEC) -n 8 ./$(TEST_NAME) 8 1 1
	$(MPIEXEC) -n 8 ./$(TEST_NAME) 1 8 1
	$(MPIEXEC) -n 8 ./$(TEST_NAME) 1 1 8

run: run-positive fault

fault: $(FAULT_TEST_NAME)
	FAULT_TIMEOUT=$(FAULT_TIMEOUT) bash ./run_distribute_failures.sh \
		./$(FAULT_TEST_NAME) $(FAULT_CONTRACT_TEST_NAME)

clean:
	$(RM) $(OBJS) $(OBJ_DIR)/*.mod $(TEST_NAME)
	$(RM) $(FAULT_OBJS) $(FAULT_CONTRACT_OBJ) $(FAULT_OBJ_DIR)/*.mod \
		$(FAULT_TEST_NAME) $(FAULT_CONTRACT_TEST_NAME)
