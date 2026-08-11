.DEFAULT_GOAL := all

.PHONY: all tests run clean

SRC_DIR   ?= ../../src
FC        := gfortran
OBJ_DIR   := _OBJ
TEST_NAME := test_gnuplot

FFLAGS += -std=f2003 -fimplicit-none -pedantic -O0 -g
FFLAGS += -fcheck=all -fbounds-check -finit-real=snan
FFLAGS += -cpp -J$(OBJ_DIR) -I$(OBJ_DIR) -I$(SRC_DIR)

OBJS := \
	$(OBJ_DIR)/plot_stub.o \
	$(OBJ_DIR)/gnuplot.o \
	$(OBJ_DIR)/test_gnuplot.o

$(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/plot_stub.o: plot_stub.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/gnuplot.o: $(SRC_DIR)/gnuplot.F90 $(OBJ_DIR)/plot_stub.o
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/test_gnuplot.o: test_gnuplot.F90 $(OBJ_DIR)/gnuplot.o
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_NAME): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

all: $(TEST_NAME)

tests: all

run: $(TEST_NAME)
	./$(TEST_NAME)

clean:
	$(RM) $(OBJS) $(OBJ_DIR)/*.mod $(TEST_NAME)
	$(RM) gnuplot_single_12.plot gnuplot_stack_0.plot gnuplot_stack_1.plot
