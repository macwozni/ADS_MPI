.DEFAULT_GOAL := all

.PHONY: all tests run run-positive run-mpi-failure clean

SRC_DIR   ?= ../../src
FC        := gfortran
OBJ_DIR   := _OBJ
TEST_NAME := test_solution_output
ERROR_CONTRACT_NAME := solution_output_error_contract_probe
ERROR_CONTRACT_OBJ := $(OBJ_DIR)/solution_output_error_contract_probe.o
TIMEOUT   ?= timeout 10s

FFLAGS += -std=f2003 -fimplicit-none -pedantic -O0 -g
FFLAGS += -fcheck=all -fbounds-check -finit-real=snan
FFLAGS += -cpp -J$(OBJ_DIR) -I$(OBJ_DIR) -I$(SRC_DIR)

OBJS := \
	$(OBJ_DIR)/setup_stub.o \
	$(OBJ_DIR)/parallelism_stub.o \
	$(OBJ_DIR)/mpi_stub.o \
	$(OBJ_DIR)/plot_stub.o \
	$(OBJ_DIR)/vtk_stub.o \
	$(OBJ_DIR)/my_mpi_stub.o \
	$(OBJ_DIR)/solution_output.o \
	$(OBJ_DIR)/test_solution_output.o

$(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/setup_stub.o: setup_stub.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/parallelism_stub.o: parallelism_stub.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/mpi_stub.o: mpi_stub.F90 $(OBJ_DIR)/parallelism_stub.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/plot_stub.o: plot_stub.F90 $(OBJ_DIR)/parallelism_stub.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/vtk_stub.o: vtk_stub.F90 $(OBJ_DIR)/plot_stub.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/my_mpi_stub.o: my_mpi_stub.F90 $(OBJ_DIR)/setup_stub.o \
		$(OBJ_DIR)/parallelism_stub.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/solution_output.o: $(SRC_DIR)/solution_output.F90 \
		$(OBJ_DIR)/setup_stub.o $(OBJ_DIR)/parallelism_stub.o \
		$(OBJ_DIR)/mpi_stub.o $(OBJ_DIR)/plot_stub.o $(OBJ_DIR)/vtk_stub.o \
		$(OBJ_DIR)/my_mpi_stub.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/test_solution_output.o: test_solution_output.F90 \
		$(OBJ_DIR)/solution_output.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ERROR_CONTRACT_OBJ): solution_output_error_contract_probe.F90 \
		$(OBJ_DIR)/solution_output.o
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_NAME): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(ERROR_CONTRACT_NAME): \
		$(filter-out $(OBJ_DIR)/test_solution_output.o,$(OBJS)) \
		$(ERROR_CONTRACT_OBJ)
	$(FC) $(FFLAGS) -o $@ $^

all: $(TEST_NAME)

tests: all

run-positive: $(TEST_NAME)
	$(TIMEOUT) ./$(TEST_NAME) positive

run-mpi-failure: $(TEST_NAME)
	FAULT_TIMEOUT=10s bash ./run_solution_output_failures.sh \
		./$(TEST_NAME) $(ERROR_CONTRACT_NAME)

run: run-positive run-mpi-failure

clean:
	$(RM) $(OBJS) $(ERROR_CONTRACT_OBJ) $(OBJ_DIR)/*.mod \
		$(TEST_NAME) $(ERROR_CONTRACT_NAME)
