.DEFAULT_GOAL := all

.PHONY: all tests run clean

SRC_DIR   ?= ../../src
FC        := gfortran
OBJ_DIR   := _OBJ
TEST_NAME := test_setup

FFLAGS += -std=f2008 -pedantic -fimplicit-none -O0 -g
FFLAGS += -fcheck=all -fbounds-check -finit-real=snan
FFLAGS += -J$(OBJ_DIR) -I$(OBJ_DIR)

OBJS := $(OBJ_DIR)/Setup.o $(OBJ_DIR)/test_setup.o

$(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/Setup.o: $(SRC_DIR)/Setup.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/test_setup.o: test_setup.F90 $(OBJ_DIR)/Setup.o | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_NAME): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

all: $(TEST_NAME)

tests: all

run: $(TEST_NAME)
	./$(TEST_NAME)

clean:
	$(RM) $(OBJS) $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.smod $(TEST_NAME)
