.DEFAULT_GOAL := all

.PHONY: all tests run clean

PFUNIT_ROOT ?= /opt/lib/pfunit/PFUNIT-4.16
SRC_DIR     ?= ../../src

FC := /opt/lib/mpich-5.0.0/bin/mpif90

include $(PFUNIT_ROOT)/include/PFUNIT.mk

OBJ_DIR := _OBJ
FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I. -I$(SRC_DIR) -I$(OBJ_DIR) -J$(OBJ_DIR)
FFLAGS += -O0 -g -fcheck=all -fbounds-check

STRICT_FFLAGS := -std=f2003 -pedantic -fimplicit-none
TEST_NAME     := test_rhs_eq
PROBE         := rhs_eq_probe
PROD_OBJS     := $(OBJ_DIR)/Setup.o $(OBJ_DIR)/Interfaces.o $(OBJ_DIR)/RHS_eq.o

$(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/Setup.o: $(SRC_DIR)/Setup.F90 | $(OBJ_DIR)
	$(FC) -c $(FFLAGS) $(STRICT_FFLAGS) $< -o $@

$(OBJ_DIR)/Interfaces.o: $(SRC_DIR)/Interfaces.F90 $(OBJ_DIR)/Setup.o
	$(FC) -c $(FFLAGS) $(STRICT_FFLAGS) $< -o $@

$(OBJ_DIR)/RHS_eq.o: $(SRC_DIR)/RHS_eq.F90 $(OBJ_DIR)/Setup.o $(OBJ_DIR)/Interfaces.o
	$(FC) -c $(FFLAGS) $(STRICT_FFLAGS) $< -o $@

$(OBJ_DIR)/$(PROBE).o: $(PROBE).f90 $(PROD_OBJS)
	$(FC) -c $(FFLAGS) $(STRICT_FFLAGS) $< -o $@

$(PROBE): $(OBJ_DIR)/$(PROBE).o $(PROD_OBJS)
	$(FC) -o $@ $^

$(TEST_NAME)_driver.o: | $(OBJ_DIR)
$(TEST_NAME).o: $(PROD_OBJS) | $(OBJ_DIR)

%.o: %.F90
	$(FC) -c $(FFLAGS) $< -o $@

$(TEST_NAME)_TESTS := $(TEST_NAME).pf
$(TEST_NAME)_REGISTRY :=
$(TEST_NAME)_OTHER_SOURCES := $(PROD_OBJS)
$(TEST_NAME)_OTHER_LIBRARIES :=
$(TEST_NAME)_OTHER_INCS := -I. -I$(SRC_DIR) -I$(OBJ_DIR)

$(eval $(call make_pfunit_test,$(TEST_NAME)))

all: $(TEST_NAME) $(PROBE)

tests: all

run: all
	./$(TEST_NAME)
	bash ./test_rhs_eq_cli.sh ./$(PROBE)

clean:
	$(RM) *.o *.mod *.smod *.inc $(TEST_NAME).F90 $(TEST_NAME) $(PROBE)
	$(RM) $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.smod
