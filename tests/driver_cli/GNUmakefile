.DEFAULT_GOAL := all

.NOTPARALLEL:

.PHONY: all tests run run-cli run-smoke run-integration clean

MYMAKE_DIR ?= ../../mymake
SRC_DIR    ?= ../../src
MPIEXEC    ?= /opt/lib/mpich-5.0.0/bin/mpiexec
MPIFC      ?= /opt/lib/mpich-5.0.0/bin/mpif90
EXEC_DIR   := $(MYMAKE_DIR)/EXEC
MYMAKE_COMPILER := COMPILER='$(MPIFC)'

COMMON_SOURCES := $(wildcard $(SRC_DIR)/*.F90)
BUILD_INPUTS   := GNUmakefile $(COMMON_SOURCES) $(MYMAKE_DIR)/makefile \
                  $(MYMAKE_DIR)/m_files $(MYMAKE_DIR)/m_options

DRIVERS := \
   $(EXEC_DIR)/l2 \
   $(EXEC_DIR)/heat \
   $(EXEC_DIR)/eriksson \
   $(EXEC_DIR)/pure_diffusion_igrm \
   $(EXEC_DIR)/oil \
   $(EXEC_DIR)/igrm_l2 \
   $(EXEC_DIR)/igrm_heat

# A driver test is also a clean-build test. This deliberately avoids reusing
# an executable produced with a different compiler or build configuration.
.PHONY: $(DRIVERS)

$(EXEC_DIR)/l2: $(BUILD_INPUTS) $(wildcard ../../problems/L2Projection/*.F90)
	$(MAKE) -j1 -C $(MYMAKE_DIR) clean SOURCE_ALL='$$(SOURCES) $$(L2)' EXEC=l2
	$(MAKE) -j1 -C $(MYMAKE_DIR) $(MYMAKE_COMPILER) SOURCE_ALL='$$(SOURCES) $$(L2)' EXEC=l2

$(EXEC_DIR)/heat: $(BUILD_INPUTS) $(wildcard ../../problems/heat/*.F90)
	$(MAKE) -j1 -C $(MYMAKE_DIR) clean SOURCE_ALL='$$(SOURCES) $$(HEAT)' EXEC=heat
	$(MAKE) -j1 -C $(MYMAKE_DIR) $(MYMAKE_COMPILER) SOURCE_ALL='$$(SOURCES) $$(HEAT)' EXEC=heat

$(EXEC_DIR)/eriksson: $(BUILD_INPUTS) $(wildcard ../../problems/eriksson/*.F90)
	$(MAKE) -j1 -C $(MYMAKE_DIR) clean SOURCE_ALL='$$(SOURCES) $$(ERIKSSON)' EXEC=eriksson
	$(MAKE) -j1 -C $(MYMAKE_DIR) $(MYMAKE_COMPILER) SOURCE_ALL='$$(SOURCES) $$(ERIKSSON)' EXEC=eriksson

$(EXEC_DIR)/pure_diffusion_igrm: $(BUILD_INPUTS) $(wildcard ../../problems/pure_diffusion_igrm/*.F90)
	$(MAKE) -j1 -C $(MYMAKE_DIR) clean SOURCE_ALL='$$(SOURCES) $$(PURE_DIFFUSION_IGRM)' EXEC=pure_diffusion_igrm
	$(MAKE) -j1 -C $(MYMAKE_DIR) $(MYMAKE_COMPILER) SOURCE_ALL='$$(SOURCES) $$(PURE_DIFFUSION_IGRM)' EXEC=pure_diffusion_igrm

$(EXEC_DIR)/oil: $(BUILD_INPUTS) $(wildcard ../../problems/oil/*.F90)
	$(MAKE) -j1 -C $(MYMAKE_DIR) clean SOURCE_ALL='$$(SOURCES) $$(OIL)' EXEC=oil
	$(MAKE) -j1 -C $(MYMAKE_DIR) $(MYMAKE_COMPILER) SOURCE_ALL='$$(SOURCES) $$(OIL)' EXEC=oil

$(EXEC_DIR)/igrm_l2: $(BUILD_INPUTS) $(wildcard ../../problems/igrm_L2Projection/*.F90)
	$(MAKE) -j1 -C $(MYMAKE_DIR) clean SOURCE_ALL='$$(SOURCES) $$(IGRM_L2)' EXEC=igrm_l2
	$(MAKE) -j1 -C $(MYMAKE_DIR) $(MYMAKE_COMPILER) SOURCE_ALL='$$(SOURCES) $$(IGRM_L2)' EXEC=igrm_l2

$(EXEC_DIR)/igrm_heat: $(BUILD_INPUTS) $(wildcard ../../problems/igrm_heat/*.F90)
	$(MAKE) -j1 -C $(MYMAKE_DIR) clean SOURCE_ALL='$$(SOURCES) $$(IGRM_HEAT)' EXEC=igrm_heat
	$(MAKE) -j1 -C $(MYMAKE_DIR) $(MYMAKE_COMPILER) SOURCE_ALL='$$(SOURCES) $$(IGRM_HEAT)' EXEC=igrm_heat

all: $(DRIVERS)

tests: all

run-cli: all
	MPIEXEC='$(MPIEXEC)' EXEC_DIR='$(abspath $(EXEC_DIR))' bash ./test_driver_cli.sh

run-smoke: all
	MPIEXEC='$(MPIEXEC)' EXEC_DIR='$(abspath $(EXEC_DIR))' bash ./test_driver_smoke.sh

run-integration: all
	MPIEXEC='$(MPIEXEC)' EXEC_DIR='$(abspath $(EXEC_DIR))' python3 ./test_driver_integration.py

run: run-cli run-smoke run-integration

clean:
	@:
