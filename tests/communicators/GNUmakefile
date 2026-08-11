.DEFAULT_GOAL := all

.PHONY: all tests run run-positive run-stub clean

SRC_DIR   ?= ../../src
MPIEXEC   ?= /opt/lib/mpich-5.0.0/bin/mpiexec
FC        := /opt/lib/mpich-5.0.0/bin/mpif90
STUB_FC   := gfortran
OBJ_DIR   := _OBJ
STUB_OBJ_DIR := _STUB_OBJ
TEST_NAME := test_communicators
ERROR_PROBE := communicators_error_probe

FFLAGS += -std=f2003 -fimplicit-none -pedantic -O0 -g
FFLAGS += -fcheck=all -fbounds-check
FFLAGS += -J$(OBJ_DIR) -I$(OBJ_DIR) -I$(SRC_DIR)

STUB_FFLAGS := -std=f2003 -fimplicit-none -pedantic -O0 -g
STUB_FFLAGS += -fcheck=all -fbounds-check
STUB_FFLAGS += -J$(STUB_OBJ_DIR) -I$(STUB_OBJ_DIR) -I$(SRC_DIR)

OBJS := \
	$(OBJ_DIR)/parallelism.o \
	$(OBJ_DIR)/communicators.o \
	$(OBJ_DIR)/test_communicators.o

STUB_OBJS := \
	$(STUB_OBJ_DIR)/communicators_mpi_stub.o \
	$(STUB_OBJ_DIR)/communicators_parallelism_stub.o \
	$(STUB_OBJ_DIR)/communicators.o \
	$(STUB_OBJ_DIR)/communicators_error_probe.o

$(OBJ_DIR):
	mkdir -p $@

$(STUB_OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/parallelism.o: $(SRC_DIR)/parallelism.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/communicators.o: $(SRC_DIR)/communicators.F90 \
		$(OBJ_DIR)/parallelism.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/test_communicators.o: test_communicators.F90 \
		$(OBJ_DIR)/parallelism.o $(OBJ_DIR)/communicators.o
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_NAME): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

$(STUB_OBJ_DIR)/communicators_mpi_stub.o: communicators_mpi_stub.F90 | $(STUB_OBJ_DIR)
	$(STUB_FC) $(STUB_FFLAGS) -c $< -o $@

$(STUB_OBJ_DIR)/communicators_parallelism_stub.o: \
		communicators_parallelism_stub.F90 | $(STUB_OBJ_DIR)
	$(STUB_FC) $(STUB_FFLAGS) -c $< -o $@

$(STUB_OBJ_DIR)/communicators.o: $(SRC_DIR)/communicators.F90 \
		$(STUB_OBJ_DIR)/communicators_mpi_stub.o \
		$(STUB_OBJ_DIR)/communicators_parallelism_stub.o
	$(STUB_FC) $(STUB_FFLAGS) -c $< -o $@

$(STUB_OBJ_DIR)/communicators_error_probe.o: communicators_error_probe.F90 \
		$(STUB_OBJ_DIR)/communicators_mpi_stub.o \
		$(STUB_OBJ_DIR)/communicators.o
	$(STUB_FC) $(STUB_FFLAGS) -c $< -o $@

$(ERROR_PROBE): $(STUB_OBJS)
	$(STUB_FC) $(STUB_FFLAGS) -o $@ $^

all: $(TEST_NAME) $(ERROR_PROBE)

tests: all

run-positive: $(TEST_NAME)
	$(MPIEXEC) -n 1 ./$(TEST_NAME) 1 1 1
	$(MPIEXEC) -n 6 ./$(TEST_NAME) 3 2 1
	$(MPIEXEC) -n 6 ./$(TEST_NAME) 2 1 3
	$(MPIEXEC) -n 8 ./$(TEST_NAME) 2 2 2

run-stub: $(ERROR_PROBE)
	PROBE=./$(ERROR_PROBE) ./test_communicators_cli.sh

run: run-positive run-stub

clean:
	$(RM) $(OBJS) $(OBJ_DIR)/*.mod $(TEST_NAME)
	$(RM) $(STUB_OBJS) $(STUB_OBJ_DIR)/*.mod $(ERROR_PROBE)
