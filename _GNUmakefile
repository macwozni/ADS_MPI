.PHONY: tests clean all
.DEFAULT_GOAL = tests


TOP_DIR := $(shell pwd)
SRC_DIR=$(TOP_DIR)/src
TEST_DIR=$(TOP_DIR)/tests

VPATH = . $(SRC_DIR) $(TEST_DIR)

include $(PFUNIT)/include/base.mk

EXE = tests$(EXE_EXT)
ifneq ($(UNAME),Windows)
	LIBS = -L$(PFUNIT)/lib -lpfunit 
else
	LIBS = $(PFUNIT)/lib/libpfunit$(LIB_EXT)
endif

all: $(EXE)


ifeq ($(USEMPI),YES)
	mpirun -np 1 ./$(EXE) -xml tests.xml
else
	./$(EXE) -xml tests.xml
endif

SUT:
	make -C $(SRC_DIR) SUT
	make -C $(TEST_DIR) tests

tests: all

$(EXE): testSuites.inc test_knot.pf knot_vector.F90 SUT
	$(FC) -o $@ -I$(PFUNIT)/mod -I$(PFUNIT)/include -Itests $(PFUNIT)/include/driver.F90 $(TEST_DIR)/*$(OBJ_EXT) $(SRC_DIR)/*$(OBJ_EXT) $(LIBS) $(FFLAGS) $(FPPFLAGS)

distclean: clean

clean: local-E0-clean

local-E0-clean:
	make -C $(SRC_DIR) clean
	make -C $(TEST_DIR) clean
	rm -f $(EXE) *$(OBJ_EXT) tests.xml

ifeq ($(UNAME),Windows)
	export PFUNIT
endif

export FC
export FPPFLAGS
export FFLAGS
export SRC_DIR
export TEST_DIR
export OBJ_EXT
export LIB_EXT
export EXE_EXT

