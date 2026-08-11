.DEFAULT_GOAL := all

.PHONY: all tests run clean

SRC_DIR   ?= ../../src
MPIEXEC   ?= /opt/lib/mpich-5.0.0/bin/mpiexec
FC        := /opt/lib/mpich-5.0.0/bin/mpif90
OBJ_DIR   := _OBJ
TEST_NAME := test_plot

FFLAGS += -std=f2003 -fimplicit-none -pedantic -O0 -g
FFLAGS += -fcheck=all -fbounds-check -fopenmp -finit-real=snan
FFLAGS += -cpp -J$(OBJ_DIR) -I$(OBJ_DIR) -I$(SRC_DIR)

OBJS := \
	$(OBJ_DIR)/math_stub.o \
	$(OBJ_DIR)/basis_stub.o \
	$(OBJ_DIR)/plot.o \
	$(OBJ_DIR)/plot_test_callbacks.o \
	$(OBJ_DIR)/test_plot.o

$(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/math_stub.o: math_stub.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/basis_stub.o: basis_stub.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/plot.o: $(SRC_DIR)/plot.F90 $(OBJ_DIR)/math_stub.o \
		$(OBJ_DIR)/basis_stub.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/plot_test_callbacks.o: plot_test_callbacks.F90 $(OBJ_DIR)/plot.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/test_plot.o: test_plot.F90 $(OBJ_DIR)/plot.o \
		$(OBJ_DIR)/plot_test_callbacks.o
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_NAME): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

all: $(TEST_NAME)

tests: all

run: $(TEST_NAME)
	OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 $(MPIEXEC) -n 1 ./$(TEST_NAME)
	OMP_DYNAMIC=FALSE OMP_NUM_THREADS=4 $(MPIEXEC) -n 4 ./$(TEST_NAME)

clean:
	$(RM) $(OBJS) $(OBJ_DIR)/*.mod $(TEST_NAME)
