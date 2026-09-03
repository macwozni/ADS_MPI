.DEFAULT_GOAL := all

.PHONY: all tests run clean

SRC_DIR   ?= ../../src
FC        := gfortran
OBJ_DIR   := _OBJ
TEST_NAME := test_ads

FFLAGS += -std=f2008 -pedantic -O0 -g -cpp
FFLAGS += -fcheck=all -fbounds-check -finit-real=snan
FFLAGS += -ffree-line-length-none -J$(OBJ_DIR) -I$(OBJ_DIR) -I$(SRC_DIR)

OBJS := \
   $(OBJ_DIR)/Setup.o \
   $(OBJ_DIR)/Interfaces.o \
   $(OBJ_DIR)/workflow_test_doubles.o \
   $(OBJ_DIR)/ADS.o \
   $(OBJ_DIR)/test_ads.o

MODULES := $(addprefix $(OBJ_DIR)/, \
   setup.mod \
   interfaces.mod \
   workflow_test_support.mod \
   mpi.mod \
   sparse.mod \
   ads_lifecycle.mod \
   mumps_solver.mod \
   solution_output.mod \
   projection_engine.mod \
   ads_directional_solve.mod \
   reorderrhs.mod \
   my_mpi.mod \
   adss.mod)

# Object name used before the ADS error-propagation test was consolidated.
STALE_OBJS := $(OBJ_DIR)/test_ads_error_propagation.o

$(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/Setup.o: $(SRC_DIR)/Setup.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/Interfaces.o: $(SRC_DIR)/Interfaces.F90 $(OBJ_DIR)/Setup.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/workflow_test_doubles.o: workflow_test_doubles.F90 \
                                    $(OBJ_DIR)/Setup.o $(OBJ_DIR)/Interfaces.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/ADS.o: $(SRC_DIR)/ADS.F90 $(OBJ_DIR)/workflow_test_doubles.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/test_ads.o: test_ads.F90 $(OBJ_DIR)/ADS.o
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_NAME): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

all: $(TEST_NAME)

tests: all

run: $(TEST_NAME)
	./$(TEST_NAME)

clean:
	$(RM) -- $(OBJS) $(STALE_OBJS) $(MODULES) $(TEST_NAME)
