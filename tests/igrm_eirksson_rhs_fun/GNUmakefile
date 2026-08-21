.DEFAULT_GOAL := all

.PHONY: all tests run clean

CONFIG      ?= ../../m_options
CONFIG_PATH := $(abspath $(CONFIG))
include $(CONFIG_PATH)

PFUNIT_ROOT ?= /opt/lib/pfunit/PFUNIT-4.16
SRC_DIR     ?= ../../src
PROBLEM_DIR ?= ../../problems/igrm_eirksson

FC := $(MPIFC)

include $(PFUNIT_ROOT)/include/PFUNIT.mk

OBJ_DIR := _OBJ
FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I. -I$(OBJ_DIR) -I$(SRC_DIR) -I$(PROBLEM_DIR) -J$(OBJ_DIR)
FFLAGS += -O0 -g -fcheck=all -fbounds-check -finit-real=snan

STRICT_FFLAGS := -std=f2003 -pedantic -fimplicit-none
TEST_NAME     := test_igrm_eirksson_rhs_fun
PROD_OBJS     := $(OBJ_DIR)/input_data_stub.o $(OBJ_DIR)/RHS_fun.o

$(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/input_data_stub.o: input_data_stub.F90 | $(OBJ_DIR)
	$(FC) -c $(FFLAGS) $(STRICT_FFLAGS) $< -o $@

$(OBJ_DIR)/RHS_fun.o: $(PROBLEM_DIR)/RHS_fun.F90 \
                      $(OBJ_DIR)/input_data_stub.o
	$(FC) -c $(FFLAGS) $(STRICT_FFLAGS) $< -o $@

$(TEST_NAME)_driver.o: | $(OBJ_DIR)
$(TEST_NAME).o: $(PROD_OBJS) | $(OBJ_DIR)

%.o: %.F90
	$(FC) -c $(FFLAGS) $< -o $@

$(TEST_NAME)_TESTS := $(TEST_NAME).pf
$(TEST_NAME)_REGISTRY :=
$(TEST_NAME)_OTHER_SOURCES := $(PROD_OBJS)
$(TEST_NAME)_OTHER_LIBRARIES :=
$(TEST_NAME)_OTHER_INCS := -I. -I$(OBJ_DIR) -I$(SRC_DIR) -I$(PROBLEM_DIR)

$(eval $(call make_pfunit_test,$(TEST_NAME)))

all: $(TEST_NAME)

tests: all

run: all
	./$(TEST_NAME)

clean:
	$(RM) *.o *.mod *.smod *.inc $(TEST_NAME).F90 $(TEST_NAME)
	$(RM) $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.smod
