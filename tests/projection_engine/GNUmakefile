.DEFAULT_GOAL := all

.PHONY: all tests run clean

SRC_DIR   ?= ../../src
FC        := gfortran
OBJ_DIR   := _OBJ
TEST_NAME := test_projection_engine

FFLAGS += -std=f2008 -pedantic -fimplicit-none -O0 -g
FFLAGS += -fcheck=all -fbounds-check
FFLAGS += -J$(OBJ_DIR) -I$(OBJ_DIR)

OBJS := $(OBJ_DIR)/projection_dependencies_stub.o \
	$(OBJ_DIR)/projection_engine.o \
	$(OBJ_DIR)/test_projection_engine.o

$(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/projection_dependencies_stub.o: projection_dependencies_stub.F90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/projection_engine.o: $(SRC_DIR)/projection_engine.F90 \
		$(OBJ_DIR)/projection_dependencies_stub.o | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/test_projection_engine.o: test_projection_engine.F90 \
		$(OBJ_DIR)/projection_engine.o | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_NAME): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

all: $(TEST_NAME)

tests: all

run: $(TEST_NAME)
	./$(TEST_NAME)

clean:
	$(RM) $(OBJS) $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.smod $(TEST_NAME)
