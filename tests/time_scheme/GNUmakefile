.DEFAULT_GOAL := all

.PHONY: all tests run clean

PFUNIT_ROOT ?= /opt/lib/pfunit/PFUNIT-4.16
SRC_DIR     ?= ../../src

FC = /opt/lib/mpich-5.0.0/bin/mpif90

include $(PFUNIT_ROOT)/include/PFUNIT.mk

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I. -I$(SRC_DIR)

TEST_NAME := test_time_scheme
PROBE := time_scheme_validation_probe
PROD_OBJS := setup.o interfaces.o gauss.o time_scheme_test_doubles.o time_scheme.o

setup.o: $(SRC_DIR)/Setup.F90
	$(FC) -c $(FFLAGS) $< -o $@

interfaces.o: $(SRC_DIR)/Interfaces.F90 setup.o
	$(FC) -c $(FFLAGS) $< -o $@

gauss.o: $(SRC_DIR)/gauss.F90
	$(FC) -c $(FFLAGS) $< -o $@

time_scheme_test_doubles.o: time_scheme_test_doubles.f90 setup.o interfaces.o
	$(FC) -c $(FFLAGS) $< -o $@

time_scheme.o: $(SRC_DIR)/time_scheme.F90 gauss.o time_scheme_test_doubles.o
	$(FC) -c $(FFLAGS) $< -o $@

$(TEST_NAME).o: $(PROD_OBJS)
$(PROBE).o: $(PROBE).f90 $(PROD_OBJS)

$(PROBE): $(PROBE).o $(PROD_OBJS)
	$(FC) -o $@ $^

%.o: %.F90
	$(FC) -c $(FFLAGS) $< -o $@

%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@

$(TEST_NAME)_TESTS := $(TEST_NAME).pf
$(TEST_NAME)_REGISTRY :=
$(TEST_NAME)_OTHER_SOURCES := $(PROD_OBJS)
$(TEST_NAME)_OTHER_LIBRARIES :=
$(TEST_NAME)_OTHER_INCS := -I. -I$(SRC_DIR)

$(eval $(call make_pfunit_test,$(TEST_NAME)))

all: $(TEST_NAME) $(PROBE)

tests: $(TEST_NAME)

run: all
	./$(TEST_NAME)
	bash ./test_time_scheme_validation_cli.sh ./$(PROBE)

clean:
	$(RM) *.o *.mod *.smod *.inc $(TEST_NAME).F90 $(TEST_NAME) $(PROBE)
