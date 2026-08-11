.DEFAULT_GOAL := all

.PHONY: all tests run clean

SRC_DIR   ?= ../../src
FC        := gfortran
OBJ_DIR   := _OBJ
TEST_NAME := test_vtk

FFLAGS += -std=f2003 -fimplicit-none -pedantic -O0 -g
FFLAGS += -fcheck=all -fbounds-check -finit-real=snan
FFLAGS += -cpp -J$(OBJ_DIR) -I$(OBJ_DIR) -I$(SRC_DIR)

OBJS := \
	$(OBJ_DIR)/plot_stub.o \
	$(OBJ_DIR)/vtk.o \
	$(OBJ_DIR)/test_vtk.o

$(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/plot_stub.o: plot_stub.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/vtk.o: $(SRC_DIR)/vtk.F90 $(OBJ_DIR)/plot_stub.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/test_vtk.o: test_vtk.F90 $(OBJ_DIR)/vtk.o
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_NAME): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

all: $(TEST_NAME)

tests: all

run: $(TEST_NAME)
	./$(TEST_NAME)

clean:
	$(RM) $(OBJS) $(OBJ_DIR)/*.mod $(TEST_NAME)
	$(RM) vtk_image.vti vtk_single.vti vtk_structured.vts
