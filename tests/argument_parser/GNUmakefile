.DEFAULT_GOAL := all

.PHONY: all tests run clean

PFUNIT_ROOT ?= /opt/lib/pfunit/PFUNIT-4.16
SRC_DIR     ?= ../../src

FC = /opt/lib/mpich-5.0.0/bin/mpif90

include $(PFUNIT_ROOT)/include/PFUNIT.mk

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I. -I$(SRC_DIR)

TEST_NAME := test_argument_parser
PROBE     := argument_parser_probe
PROD_OBJS := argument_parser.o

argument_parser.o: $(SRC_DIR)/argument_parser.F90
	$(FC) -c $(FFLAGS) $< -o $@

$(PROBE): $(PROBE).o argument_parser.o
	$(FC) -o $@ $^

$(PROBE).o: argument_parser.o
$(TEST_NAME).o: $(PROD_OBJS)

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

tests: all

run: all
	./$(TEST_NAME)
	bash ./test_argument_parser_cli.sh ./$(PROBE)

clean:
	$(RM) *.o *.mod *.smod *.inc $(TEST_NAME).F90 $(TEST_NAME) $(PROBE)
