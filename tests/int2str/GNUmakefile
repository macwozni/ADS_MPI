.DEFAULT_GOAL := all

.PHONY: all tests run clean

PFUNIT_ROOT ?= /opt/lib/pfunit/PFUNIT-4.16
SRC_DIR     ?= ../../src

FC       := /opt/lib/mpich-5.0.0/bin/mpif90
OBJ_DIR  := _OBJ
TEST_NAME := test_int2str
PROBE     := int2str_probe

include $(PFUNIT_ROOT)/include/PFUNIT.mk

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -fimplicit-none
FFLAGS += -O0 -g -fcheck=all -fbounds-check -fopenmp
FFLAGS += -J$(OBJ_DIR) -I$(OBJ_DIR) -I$(SRC_DIR)

PROD_OBJS := \
   $(OBJ_DIR)/Setup.o \
   $(OBJ_DIR)/gauss.o \
   $(OBJ_DIR)/basis.o \
   $(OBJ_DIR)/parallelism.o \
   $(OBJ_DIR)/communicators.o \
   $(OBJ_DIR)/my_mpi.o \
   $(OBJ_DIR)/projection_engine_stub.o \
   $(OBJ_DIR)/utils.o

$(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/Setup.o: $(SRC_DIR)/Setup.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/gauss.o: $(SRC_DIR)/gauss.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/basis.o: $(SRC_DIR)/basis.F90 $(OBJ_DIR)/gauss.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/parallelism.o: $(SRC_DIR)/parallelism.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/communicators.o: $(SRC_DIR)/communicators.F90 \
                            $(OBJ_DIR)/parallelism.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/my_mpi.o: $(SRC_DIR)/my_mpi.F90 \
                     $(OBJ_DIR)/parallelism.o $(OBJ_DIR)/communicators.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/projection_engine_stub.o: projection_engine_stub.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/utils.o: $(SRC_DIR)/utils.F90 \
                    $(OBJ_DIR)/Setup.o $(OBJ_DIR)/basis.o \
                    $(OBJ_DIR)/my_mpi.o $(OBJ_DIR)/parallelism.o \
                    $(OBJ_DIR)/projection_engine_stub.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/$(PROBE).o: $(PROBE).f90 $(OBJ_DIR)/utils.o
	$(FC) $(FFLAGS) -c $< -o $@

$(PROBE): $(OBJ_DIR)/$(PROBE).o $(PROD_OBJS)
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.F90
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_NAME)_driver.o: | $(OBJ_DIR)
$(TEST_NAME).o: $(PROD_OBJS)

$(TEST_NAME)_TESTS := $(TEST_NAME).pf
$(TEST_NAME)_REGISTRY :=
$(TEST_NAME)_OTHER_SOURCES := $(PROD_OBJS)
$(TEST_NAME)_OTHER_LIBRARIES :=
$(TEST_NAME)_OTHER_INCS := -I$(OBJ_DIR) -I$(SRC_DIR)

$(eval $(call make_pfunit_test,$(TEST_NAME)))

all: $(TEST_NAME) $(PROBE)

tests: all

run: all
	./$(TEST_NAME)
	bash ./test_int2str_cli.sh ./$(PROBE)

clean:
	$(RM) *.o *.mod *.smod *.inc $(TEST_NAME).F90 $(TEST_NAME) $(PROBE)
	$(RM) $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.smod
