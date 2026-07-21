.DEFAULT_GOAL := all

.PHONY: all tests run clean

SRC_DIR  ?= ../../src
MPIEXEC  ?= /opt/lib/mpich-5.0.0/bin/mpiexec
FC       := /opt/lib/mpich-5.0.0/bin/mpif90
OBJ_DIR  := _OBJ
TEST_NAME := test_norm_l2_mpi

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
   $(OBJ_DIR)/projection_engine_stub.o \
   $(OBJ_DIR)/ads_lifecycle.o \
   $(OBJ_DIR)/utils.o \
   $(OBJ_DIR)/test_norm_l2_mpi.o

$(OBJ_DIR):
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
                     $(OBJ_DIR)/parallelism.o $(OBJ_DIR)/communicators.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/projection_engine_stub.o: projection_engine_stub.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/ads_lifecycle.o: $(SRC_DIR)/ads_lifecycle.F90 \
                            $(OBJ_DIR)/Setup.o $(OBJ_DIR)/knot_vector.o \
                            $(OBJ_DIR)/parallelism.o $(OBJ_DIR)/basis.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/utils.o: $(SRC_DIR)/utils.F90 \
                    $(OBJ_DIR)/Setup.o $(OBJ_DIR)/basis.o \
                    $(OBJ_DIR)/my_mpi.o $(OBJ_DIR)/parallelism.o \
                    $(OBJ_DIR)/projection_engine_stub.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/test_norm_l2_mpi.o: test_norm_l2_mpi.F90 \
                              $(OBJ_DIR)/Setup.o $(OBJ_DIR)/basis.o \
                              $(OBJ_DIR)/parallelism.o $(OBJ_DIR)/communicators.o \
                              $(OBJ_DIR)/ads_lifecycle.o $(OBJ_DIR)/utils.o
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_NAME): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

all: $(TEST_NAME)

tests: all

run: $(TEST_NAME)
	$(MPIEXEC) -n 1 ./$(TEST_NAME) 1 1 1
	$(MPIEXEC) -n 8 ./$(TEST_NAME) 2 2 2

clean:
	$(RM) $(OBJS) $(OBJ_DIR)/*.mod $(TEST_NAME)
