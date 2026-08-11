.DEFAULT_GOAL := all

.PHONY: all tests run run-positive run-cli clean

SRC_DIR   ?= ../../src
MPIEXEC   ?= /opt/lib/mpich-5.0.0/bin/mpiexec
FC        := /opt/lib/mpich-5.0.0/bin/mpif90
OBJ_DIR   := _OBJ
TEST_NAME := test_ads_lifecycle
PROBE     := ads_lifecycle_probe

FFLAGS += -std=f2003 -fimplicit-none -pedantic -O0 -g
FFLAGS += -fcheck=all -fbounds-check -fopenmp -finit-real=snan
FFLAGS += -J$(OBJ_DIR) -I$(OBJ_DIR) -I$(SRC_DIR)

OBJS := \
   $(OBJ_DIR)/Setup.o \
   $(OBJ_DIR)/knot_vector.o \
   $(OBJ_DIR)/parallelism.o \
   $(OBJ_DIR)/gauss.o \
   $(OBJ_DIR)/basis.o \
   $(OBJ_DIR)/ads_lifecycle.o \
   $(OBJ_DIR)/test_ads_lifecycle.o

PROBE_OBJ := $(OBJ_DIR)/ads_lifecycle_probe.o

$(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/Setup.o: $(SRC_DIR)/Setup.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/knot_vector.o: $(SRC_DIR)/knot_vector.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/parallelism.o: $(SRC_DIR)/parallelism.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/gauss.o: $(SRC_DIR)/gauss.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/basis.o: $(SRC_DIR)/basis.F90 $(OBJ_DIR)/gauss.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/ads_lifecycle.o: $(SRC_DIR)/ads_lifecycle.F90 \
                            $(OBJ_DIR)/Setup.o $(OBJ_DIR)/knot_vector.o \
                            $(OBJ_DIR)/parallelism.o $(OBJ_DIR)/basis.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/test_ads_lifecycle.o: test_ads_lifecycle.F90 \
                                 $(OBJ_DIR)/Setup.o \
                                 $(OBJ_DIR)/parallelism.o \
                                 $(OBJ_DIR)/ads_lifecycle.o
	$(FC) $(FFLAGS) -c $< -o $@

$(PROBE_OBJ): ads_lifecycle_probe.F90 \
              $(OBJ_DIR)/Setup.o \
              $(OBJ_DIR)/parallelism.o \
              $(OBJ_DIR)/ads_lifecycle.o
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_NAME): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(PROBE): $(filter-out $(OBJ_DIR)/test_ads_lifecycle.o,$(OBJS)) $(PROBE_OBJ)
	$(FC) $(FFLAGS) -o $@ $^

all: $(TEST_NAME) $(PROBE)

tests: all

run-positive: $(TEST_NAME)
	OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 $(MPIEXEC) -n 1 ./$(TEST_NAME) 1 1 1
	OMP_DYNAMIC=FALSE OMP_NUM_THREADS=2 $(MPIEXEC) -n 6 ./$(TEST_NAME) 3 2 1
	OMP_DYNAMIC=FALSE OMP_NUM_THREADS=4 $(MPIEXEC) -n 8 ./$(TEST_NAME) 2 2 2

run-cli: $(PROBE)
	MPIEXEC=$(MPIEXEC) PROBE=./$(PROBE) ./test_ads_lifecycle_cli.sh

run: run-positive run-cli

clean:
	$(RM) $(OBJS) $(PROBE_OBJ) $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.smod \
		$(TEST_NAME) $(PROBE)
