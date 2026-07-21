.DEFAULT_GOAL := all

.PHONY: all tests run clean

PFUNIT_ROOT ?= /opt/lib/pfunit/PFUNIT-4.16
SRC_DIR     ?= ../../src

FC = /opt/lib/mpich-5.0.0/bin/mpif90

include $(PFUNIT_ROOT)/include/PFUNIT.mk

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I. -I$(SRC_DIR)
FFLAGS += -O0 -g -fcheck=all -fbounds-check -fopenmp -finit-real=snan

TEST_NAME := test_basis
PROD_OBJS := knot_vector.o math.o gauss.o basis.o

knot_vector.o: $(SRC_DIR)/knot_vector.F90
	$(FC) -c $(FFLAGS) $< -o $@

math.o: $(SRC_DIR)/math.F90
	$(FC) -c $(FFLAGS) $< -o $@

gauss.o: $(SRC_DIR)/gauss.F90
	$(FC) -c $(FFLAGS) $< -o $@

basis.o: $(SRC_DIR)/basis.F90 gauss.o
	$(FC) -c $(FFLAGS) $< -o $@

$(TEST_NAME).o: $(PROD_OBJS)

%.o: %.F90
	$(FC) -c $(FFLAGS) $< -o $@

$(TEST_NAME)_TESTS := $(TEST_NAME).pf
$(TEST_NAME)_REGISTRY :=
$(TEST_NAME)_OTHER_SOURCES := $(PROD_OBJS)
$(TEST_NAME)_OTHER_LIBRARIES :=
$(TEST_NAME)_OTHER_INCS := -I. -I$(SRC_DIR)

$(eval $(call make_pfunit_test,$(TEST_NAME)))

all: $(TEST_NAME)

tests: $(TEST_NAME)

run: $(TEST_NAME)
	./$(TEST_NAME)

clean:
	$(RM) *.o *.mod *.smod *.inc $(TEST_NAME).F90 $(TEST_NAME)
