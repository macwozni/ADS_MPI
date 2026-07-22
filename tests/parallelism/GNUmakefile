.DEFAULT_GOAL := all

.PHONY: all tests run run-positive run-cli clean

SRC_DIR   ?= ../../src
MPIEXEC   ?= /opt/lib/mpich-5.0.0/bin/mpiexec
FC        := /opt/lib/mpich-5.0.0/bin/mpif90
OBJ_DIR   := _OBJ
TEST_NAME := test_parallelism
PROBE     := parallelism_probe

FFLAGS += -std=f2003 -fimplicit-none -pedantic -O0 -g
FFLAGS += -fcheck=all -fbounds-check
FFLAGS += -J$(OBJ_DIR) -I$(OBJ_DIR) -I$(SRC_DIR)

PROD_OBJ  := $(OBJ_DIR)/parallelism.o
TEST_OBJ  := $(OBJ_DIR)/test_parallelism.o
PROBE_OBJ := $(OBJ_DIR)/parallelism_probe.o

$(OBJ_DIR):
	mkdir -p $@

$(PROD_OBJ): $(SRC_DIR)/parallelism.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_OBJ): test_parallelism.F90 $(PROD_OBJ)
	$(FC) $(FFLAGS) -c $< -o $@

$(PROBE_OBJ): parallelism_probe.f90 $(PROD_OBJ)
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_NAME): $(PROD_OBJ) $(TEST_OBJ)
	$(FC) $(FFLAGS) -o $@ $^

$(PROBE): $(PROD_OBJ) $(PROBE_OBJ)
	$(FC) $(FFLAGS) -o $@ $^

all: $(TEST_NAME) $(PROBE)

tests: all

run-positive: $(TEST_NAME)
	$(MPIEXEC) -n 1 ./$(TEST_NAME) 1 1 1
	$(MPIEXEC) -n 6 ./$(TEST_NAME) 3 2 1
	$(MPIEXEC) -n 8 ./$(TEST_NAME) 2 2 2

run-cli: $(PROBE)
	MPIEXEC=$(MPIEXEC) PROBE=./$(PROBE) ./test_parallelism_cli.sh

run: run-positive run-cli

clean:
	$(RM) $(PROD_OBJ) $(TEST_OBJ) $(PROBE_OBJ)
	$(RM) $(OBJ_DIR)/*.mod $(TEST_NAME) $(PROBE)
