.DEFAULT_GOAL := all

.PHONY: all tests run clean

PFUNIT_ROOT ?= /opt/lib/pfunit/PFUNIT-4.16
SRC_DIR     ?= ../../src

FC := /opt/lib/mpich-5.0.0/bin/mpif90

include $(PFUNIT_ROOT)/include/PFUNIT.mk

OBJ_DIR := _OBJ
FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I. -I$(OBJ_DIR) -I$(SRC_DIR) -J$(OBJ_DIR)
FFLAGS += -O0 -g -fcheck=all -fbounds-check -fopenmp -finit-real=snan

TEST_NAME := test_gauss
PROD_OBJS := $(OBJ_DIR)/gauss.o

$(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/gauss.o: $(SRC_DIR)/gauss.F90 | $(OBJ_DIR)
	$(FC) -c $(FFLAGS) $< -o $@

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

all: $(TEST_NAME)

tests: all

run: all
	OMP_DYNAMIC=FALSE OMP_NUM_THREADS=8 ./$(TEST_NAME)

clean:
	$(RM) *.o *.mod *.smod *.inc $(TEST_NAME).F90 $(TEST_NAME)
	$(RM) $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.smod
