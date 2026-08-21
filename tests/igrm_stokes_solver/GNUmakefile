.DEFAULT_GOAL := all

.PHONY: all tests run clean

CONFIG      ?= ../../m_options
CONFIG_PATH := $(abspath $(CONFIG))
include $(CONFIG_PATH)

SRC_DIR     ?= ../../src
PROBLEM_DIR ?= ../../problems/igrm_stokes
MUMPS_DIR   ?= /opt/lib/MUMPS_5.8.2
MPIEXEC     ?= /opt/lib/mpich-5.0.0/bin/mpiexec
MPIFC       ?= /opt/lib/mpich-5.0.0/bin/mpif90
TIMEOUT     ?= timeout 120s

FC        := $(MPIFC)
OBJ_DIR   := _OBJ
TEST_NAME := test_igrm_stokes_solver

FFLAGS += -std=f2003 -fimplicit-none -pedantic -cpp -O0 -g
FFLAGS += -fcheck=all -fbounds-check -fopenmp -finit-real=snan
FFLAGS += -J$(OBJ_DIR) -I$(OBJ_DIR) -I$(SRC_DIR) -I$(PROBLEM_DIR)
FFLAGS += -I$(MUMPS_DIR)/include

OBJS := \
	$(OBJ_DIR)/Setup.o \
	$(OBJ_DIR)/Interfaces.o \
	$(OBJ_DIR)/knot_vector.o \
	$(OBJ_DIR)/parallelism.o \
	$(OBJ_DIR)/gauss.o \
	$(OBJ_DIR)/basis.o \
	$(OBJ_DIR)/sparse.o \
	$(OBJ_DIR)/mumps_solver.o \
	$(OBJ_DIR)/igrm_stokes_solver.o \
	$(OBJ_DIR)/test_igrm_stokes_solver.o

$(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/Setup.o: $(SRC_DIR)/Setup.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/Interfaces.o: $(SRC_DIR)/Interfaces.F90 $(OBJ_DIR)/Setup.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/knot_vector.o: $(SRC_DIR)/knot_vector.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/parallelism.o: $(SRC_DIR)/parallelism.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/gauss.o: $(SRC_DIR)/gauss.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/basis.o: $(SRC_DIR)/basis.F90 $(OBJ_DIR)/gauss.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/sparse.o: $(SRC_DIR)/sparse.F90 \
		$(MUMPS_DIR)/include/dmumps_struc.h | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/mumps_solver.o: $(SRC_DIR)/mumps_solver.F90 \
		$(OBJ_DIR)/sparse.o $(MUMPS_DIR)/include/dmumps_struc.h
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/igrm_stokes_solver.o: $(PROBLEM_DIR)/igrm_stokes_solver.F90 \
		$(OBJ_DIR)/Setup.o $(OBJ_DIR)/Interfaces.o \
		$(OBJ_DIR)/knot_vector.o $(OBJ_DIR)/parallelism.o \
		$(OBJ_DIR)/basis.o $(OBJ_DIR)/sparse.o $(OBJ_DIR)/mumps_solver.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/test_igrm_stokes_solver.o: test_igrm_stokes_solver.F90 \
		$(OBJ_DIR)/parallelism.o $(OBJ_DIR)/sparse.o \
		$(OBJ_DIR)/igrm_stokes_solver.o
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_NAME): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(USER_LIB)

all: $(TEST_NAME)

tests: all

run: all
	OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 $(TIMEOUT) \
		$(MPIEXEC) -n 1 ./$(TEST_NAME) 1
	OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 $(TIMEOUT) \
		$(MPIEXEC) -n 2 ./$(TEST_NAME) 2

clean:
	$(RM) $(OBJS) $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.smod $(TEST_NAME)
	$(RM) test_igrm_stokes_result.vti
