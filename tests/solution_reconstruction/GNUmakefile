.DEFAULT_GOAL := all

.PHONY: all tests run clean

PFUNIT_ROOT ?= /opt/lib/pfunit/PFUNIT-4.16
SRC_DIR     ?= ../../src

FC := /opt/lib/mpich-5.0.0/bin/mpif90

include $(PFUNIT_ROOT)/include/PFUNIT.mk

OBJ_DIR       := _OBJ
DEBUG_OBJ_DIR := _DEBUG_OBJ
FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I. -I$(OBJ_DIR) -I$(SRC_DIR) -J$(OBJ_DIR)
FFLAGS += -O0 -g -fcheck=all -fbounds-check -fopenmp -finit-real=snan
DEBUG_FFLAGS = $(filter-out -J%,$(FFLAGS)) -I$(DEBUG_OBJ_DIR) -J$(DEBUG_OBJ_DIR)

STRICT_FFLAGS := -std=f2003 -pedantic -fimplicit-none
TEST_NAME     := test_solution_reconstruction
PROBE         := solution_reconstruction_probe
PROD_OBJS     := $(OBJ_DIR)/Setup.o $(OBJ_DIR)/Interfaces.o \
                 $(OBJ_DIR)/solution_reconstruction.o

$(OBJ_DIR):
	mkdir -p $@

$(DEBUG_OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/Setup.o: $(SRC_DIR)/Setup.F90 | $(OBJ_DIR)
	$(FC) -c $(FFLAGS) $(STRICT_FFLAGS) $< -o $@

$(OBJ_DIR)/Interfaces.o: $(SRC_DIR)/Interfaces.F90 $(OBJ_DIR)/Setup.o
	$(FC) -c $(FFLAGS) $(STRICT_FFLAGS) $< -o $@

$(OBJ_DIR)/solution_reconstruction.o: $(SRC_DIR)/solution_reconstruction.F90 \
                                      $(OBJ_DIR)/Setup.o $(OBJ_DIR)/Interfaces.o
	$(FC) -c $(FFLAGS) $(STRICT_FFLAGS) $< -o $@

$(DEBUG_OBJ_DIR)/parallelism.o: $(SRC_DIR)/parallelism.F90 | $(DEBUG_OBJ_DIR)
	$(FC) -c $(DEBUG_FFLAGS) $(STRICT_FFLAGS) $< -o $@

$(DEBUG_OBJ_DIR)/solution_reconstruction.o: \
      $(SRC_DIR)/solution_reconstruction.F90 $(OBJ_DIR)/Setup.o \
      $(OBJ_DIR)/Interfaces.o $(DEBUG_OBJ_DIR)/parallelism.o
	$(FC) -c $(DEBUG_FFLAGS) $(STRICT_FFLAGS) -DIDEBUG $< -o $@

$(OBJ_DIR)/$(PROBE).o: $(PROBE).f90 $(PROD_OBJS)
	$(FC) -c $(FFLAGS) $(STRICT_FFLAGS) $< -o $@

$(PROBE): $(OBJ_DIR)/$(PROBE).o $(PROD_OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(TEST_NAME)_driver.o: | $(OBJ_DIR)
$(TEST_NAME).o: $(PROD_OBJS) | $(OBJ_DIR)

%.o: %.F90
	$(FC) -c $(FFLAGS) $< -o $@

$(TEST_NAME)_TESTS := $(TEST_NAME).pf
$(TEST_NAME)_REGISTRY :=
$(TEST_NAME)_OTHER_SOURCES := $(PROD_OBJS)
$(TEST_NAME)_OTHER_LIBRARIES :=
$(TEST_NAME)_OTHER_INCS := -I. -I$(OBJ_DIR) -I$(SRC_DIR)

$(eval $(call make_pfunit_test,$(TEST_NAME)))

all: $(TEST_NAME) $(PROBE) $(DEBUG_OBJ_DIR)/solution_reconstruction.o

tests: all

run: all
	OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 ./$(TEST_NAME)
	bash ./test_solution_reconstruction_cli.sh ./$(PROBE)

clean:
	$(RM) *.o *.mod *.smod *.inc $(TEST_NAME).F90 $(TEST_NAME) $(PROBE)
	$(RM) $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.smod
	$(RM) $(DEBUG_OBJ_DIR)/*.o $(DEBUG_OBJ_DIR)/*.mod $(DEBUG_OBJ_DIR)/*.smod
