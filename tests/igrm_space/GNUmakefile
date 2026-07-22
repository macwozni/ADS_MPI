.DEFAULT_GOAL := all

.PHONY: all tests run clean

PFUNIT_ROOT ?= /opt/lib/pfunit/PFUNIT-4.16
SRC_DIR     ?= ../../src

FC := /opt/lib/mpich-5.0.0/bin/mpif90

include $(PFUNIT_ROOT)/include/PFUNIT.mk

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I. -I$(SRC_DIR)
FFLAGS += -O0 -g -fcheck=all -fbounds-check

STRICT_FFLAGS := -std=f2003 -pedantic -fimplicit-none
TEST_NAME     := test_igrm_space
PROBE         := igrm_space_probe
PROD_OBJS     := igrm_space.o

igrm_space.o: $(SRC_DIR)/igrm_space.F90
	$(FC) -c $(FFLAGS) $(STRICT_FFLAGS) $< -o $@

$(PROBE): $(PROBE).o $(PROD_OBJS)
	$(FC) -o $@ $^

$(PROBE).o: $(PROBE).f90 igrm_space.o
	$(FC) -c $(FFLAGS) $(STRICT_FFLAGS) $< -o $@

$(TEST_NAME).o: $(PROD_OBJS)

%.o: %.F90
	$(FC) -c $(FFLAGS) $< -o $@

$(TEST_NAME)_TESTS := $(TEST_NAME).pf
$(TEST_NAME)_REGISTRY :=
$(TEST_NAME)_OTHER_SOURCES := $(PROD_OBJS)
$(TEST_NAME)_OTHER_LIBRARIES :=
$(TEST_NAME)_OTHER_INCS := -I. -I$(SRC_DIR)

$(eval $(call make_pfunit_test,$(TEST_NAME)))

all: $(TEST_NAME) $(PROBE)

tests: all

run: all
	./$(TEST_NAME)
	bash ./test_igrm_space_cli.sh ./$(PROBE)

clean:
	$(RM) *.o *.mod *.smod *.inc $(TEST_NAME).F90 $(TEST_NAME) $(PROBE)
