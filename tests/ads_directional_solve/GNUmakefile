.DEFAULT_GOAL := all

.PHONY: all tests run run-positive clean

SRC_DIR   ?= ../../src
MUMPS_DIR ?= /opt/lib/MUMPS_5.8.2
MPIEXEC   ?= /opt/lib/mpich-5.0.0/bin/mpiexec
TIMEOUT   ?= timeout 30s
FC        := /opt/lib/mpich-5.0.0/bin/mpif90
OBJ_DIR   := _OBJ
TEST_NAME := test_ads_directional_solve

FFLAGS += -std=f2003 -fimplicit-none -pedantic -cpp -O0 -g
FFLAGS += -fcheck=all -fbounds-check -fopenmp -finit-real=snan
FFLAGS += -J$(OBJ_DIR) -I$(OBJ_DIR) -I$(SRC_DIR) -I$(MUMPS_DIR)/include

OBJS := \
	$(OBJ_DIR)/Setup.o \
	$(OBJ_DIR)/sparse.o \
	$(OBJ_DIR)/parallelism.o \
	$(OBJ_DIR)/communicators.o \
	$(OBJ_DIR)/my_mpi.o \
	$(OBJ_DIR)/reorderRHS.o \
	$(OBJ_DIR)/directional_test_support.o \
	$(OBJ_DIR)/projection_engine_stub.o \
	$(OBJ_DIR)/mumps_solver_stub.o \
	$(OBJ_DIR)/ads_directional_solve.o \
	$(OBJ_DIR)/test_ads_directional_solve.o

$(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/Setup.o: $(SRC_DIR)/Setup.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/sparse.o: $(SRC_DIR)/sparse.F90 \
                      $(MUMPS_DIR)/include/dmumps_struc.h | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/parallelism.o: $(SRC_DIR)/parallelism.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/communicators.o: $(SRC_DIR)/communicators.F90 \
                             $(OBJ_DIR)/parallelism.o | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/my_mpi.o: $(SRC_DIR)/my_mpi.F90 \
                      $(OBJ_DIR)/parallelism.o \
                      $(OBJ_DIR)/communicators.o | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/reorderRHS.o: $(SRC_DIR)/reorderRHS.F90 \
                          $(OBJ_DIR)/Setup.o | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/directional_test_support.o: directional_test_support.F90 \
                                       $(OBJ_DIR)/Setup.o | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/projection_engine_stub.o: projection_engine_stub.F90 \
                                     $(OBJ_DIR)/sparse.o \
                                     $(OBJ_DIR)/directional_test_support.o | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/mumps_solver_stub.o: mumps_solver_stub.F90 \
                                $(OBJ_DIR)/sparse.o \
                                $(OBJ_DIR)/directional_test_support.o | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/ads_directional_solve.o: $(SRC_DIR)/ads_directional_solve.F90 \
                                    $(OBJ_DIR)/Setup.o \
                                    $(OBJ_DIR)/sparse.o \
                                    $(OBJ_DIR)/parallelism.o \
                                    $(OBJ_DIR)/communicators.o \
                                    $(OBJ_DIR)/my_mpi.o \
                                    $(OBJ_DIR)/reorderRHS.o \
                                    $(OBJ_DIR)/projection_engine_stub.o \
                                    $(OBJ_DIR)/mumps_solver_stub.o | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/test_ads_directional_solve.o: test_ads_directional_solve.F90 \
                                         $(OBJ_DIR)/Setup.o \
                                         $(OBJ_DIR)/parallelism.o \
                                         $(OBJ_DIR)/communicators.o \
                                         $(OBJ_DIR)/directional_test_support.o \
                                         $(OBJ_DIR)/ads_directional_solve.o | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_NAME): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

all: $(TEST_NAME)

tests: all

run-positive: $(TEST_NAME)
	OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 $(TIMEOUT) $(MPIEXEC) -n 1 ./$(TEST_NAME) 1 1 1
	OMP_DYNAMIC=FALSE OMP_NUM_THREADS=2 $(TIMEOUT) $(MPIEXEC) -n 6 ./$(TEST_NAME) 3 2 1
	OMP_DYNAMIC=FALSE OMP_NUM_THREADS=4 $(TIMEOUT) $(MPIEXEC) -n 8 ./$(TEST_NAME) 2 2 2
	OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 $(TIMEOUT) $(MPIEXEC) -n 8 ./$(TEST_NAME) 8 1 1
	OMP_DYNAMIC=FALSE OMP_NUM_THREADS=2 $(TIMEOUT) $(MPIEXEC) -n 8 ./$(TEST_NAME) 1 8 1
	OMP_DYNAMIC=FALSE OMP_NUM_THREADS=4 $(TIMEOUT) $(MPIEXEC) -n 8 ./$(TEST_NAME) 1 1 8

run: run-positive

clean:
	$(RM) $(OBJS) $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.smod $(TEST_NAME)
