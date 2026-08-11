.DEFAULT_GOAL := all

.PHONY: all tests run clean

SRC_DIR   ?= ../../src
MUMPS_DIR ?= /opt/lib/MUMPS_5.8.2
MPIEXEC   ?= /opt/lib/mpich-5.0.0/bin/mpiexec
FC        := /opt/lib/mpich-5.0.0/bin/mpif90
OBJ_DIR   := _OBJ
TEST_NAME := test_mumps_solver
ABORT_PROBE := legacy_abort_probe
TIMEOUT   ?= timeout 20s

FFLAGS += -std=f2003 -fimplicit-none -pedantic -cpp -O0 -g
FFLAGS += -fcheck=all -fbounds-check -fopenmp -finit-real=snan
FFLAGS += -J$(OBJ_DIR) -I$(OBJ_DIR) -I$(SRC_DIR) -I$(MUMPS_DIR)/include

OBJS := \
	$(OBJ_DIR)/mumps_test_support.o \
	$(OBJ_DIR)/sparse.o \
	$(OBJ_DIR)/mumps_solver.o \
	$(OBJ_DIR)/dmumps_stub.o \
	$(OBJ_DIR)/test_mumps_solver.o

$(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/mumps_test_support.o: mumps_test_support.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/sparse.o: $(SRC_DIR)/sparse.F90 \
                      $(MUMPS_DIR)/include/dmumps_struc.h | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/mumps_solver.o: $(SRC_DIR)/mumps_solver.F90 \
                            $(OBJ_DIR)/sparse.o \
                            $(MUMPS_DIR)/include/dmumps_struc.h | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/dmumps_stub.o: dmumps_stub.F90 \
                           $(OBJ_DIR)/mumps_test_support.o \
                           $(MUMPS_DIR)/include/dmumps_struc.h | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/test_mumps_solver.o: test_mumps_solver.F90 \
                                $(OBJ_DIR)/mumps_test_support.o \
                                $(OBJ_DIR)/sparse.o \
                                $(OBJ_DIR)/mumps_solver.o | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/$(ABORT_PROBE).o: $(ABORT_PROBE).F90 \
                              $(OBJ_DIR)/mumps_test_support.o \
                              $(OBJ_DIR)/sparse.o \
                              $(OBJ_DIR)/mumps_solver.o | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_NAME): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(ABORT_PROBE): $(filter-out $(OBJ_DIR)/test_mumps_solver.o,$(OBJS)) \
                $(OBJ_DIR)/$(ABORT_PROBE).o
	$(FC) $(FFLAGS) -o $@ $^

all: $(TEST_NAME) $(ABORT_PROBE)

tests: $(TEST_NAME)

run: all
	$(TIMEOUT) ./$(TEST_NAME)
	MPIEXEC='$(MPIEXEC)' bash ./test_legacy_abort.sh ./$(ABORT_PROBE)

clean:
	$(RM) $(OBJS) $(OBJ_DIR)/$(ABORT_PROBE).o $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.smod \
		$(TEST_NAME) $(ABORT_PROBE)
