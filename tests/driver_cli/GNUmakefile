.DEFAULT_GOAL := all

.NOTPARALLEL:

.PHONY: all tests run run-cli run-smoke run-integration clean

MYMAKE_DIR ?= ../../mymake
SRC_DIR    ?= ../../src
CONFIG     ?= ../../m_options
CONFIG_PATH := $(abspath $(CONFIG))
include $(CONFIG_PATH)

MPIEXEC    ?= /opt/lib/mpich-5.0.0/bin/mpiexec
MPIFC      ?= /opt/lib/mpich-5.0.0/bin/mpif90
PYTHON     ?= python3
DRIVER_CLI_TIMEOUT ?= 20s
DRIVER_SMOKE_TIMEOUT ?= 60s
DRIVER_INTEGRATION_TIMEOUT ?= 90
SKIP_MPI_CASES ?= 0
EXEC_DIR   := $(MYMAKE_DIR)/EXEC
MYMAKE_COMPILER := COMPILER='$(MPIFC)'

COMMON_SOURCES := $(wildcard $(SRC_DIR)/*.F90)
BUILD_INPUTS   := GNUmakefile $(COMMON_SOURCES) $(MYMAKE_DIR)/makefile \
                  $(MYMAKE_DIR)/m_files $(CONFIG_PATH)

DRIVERS := \
   $(EXEC_DIR)/l2 \
   $(EXEC_DIR)/heat \
   $(EXEC_DIR)/eriksson \
   $(EXEC_DIR)/pure_diffusion_igrm \
   $(EXEC_DIR)/oil \
   $(EXEC_DIR)/igrm_l2 \
   $(EXEC_DIR)/igrm_heat \
   $(EXEC_DIR)/igrm_eirksson \
   $(EXEC_DIR)/igrm_stokes \
   $(EXEC_DIR)/igrm_pollution

# A driver test is also a clean-build test. This deliberately avoids reusing
# an executable produced with a different compiler or build configuration.
.PHONY: $(DRIVERS)

$(EXEC_DIR)/l2: $(BUILD_INPUTS) $(wildcard ../../problems/L2Projection/*.F90)
	$(MAKE) -j1 -C $(MYMAKE_DIR) CONFIG='$(CONFIG_PATH)' $(MYMAKE_COMPILER) l2

$(EXEC_DIR)/heat: $(BUILD_INPUTS) $(wildcard ../../problems/heat/*.F90)
	$(MAKE) -j1 -C $(MYMAKE_DIR) CONFIG='$(CONFIG_PATH)' $(MYMAKE_COMPILER) heat

$(EXEC_DIR)/eriksson: $(BUILD_INPUTS) $(wildcard ../../problems/eriksson/*.F90)
	$(MAKE) -j1 -C $(MYMAKE_DIR) CONFIG='$(CONFIG_PATH)' $(MYMAKE_COMPILER) eriksson

$(EXEC_DIR)/pure_diffusion_igrm: $(BUILD_INPUTS) $(wildcard ../../problems/pure_diffusion_igrm/*.F90)
	$(MAKE) -j1 -C $(MYMAKE_DIR) CONFIG='$(CONFIG_PATH)' $(MYMAKE_COMPILER) pure_diffusion_igrm

$(EXEC_DIR)/oil: $(BUILD_INPUTS) $(wildcard ../../problems/oil/*.F90)
	$(MAKE) -j1 -C $(MYMAKE_DIR) CONFIG='$(CONFIG_PATH)' $(MYMAKE_COMPILER) oil

$(EXEC_DIR)/igrm_l2: $(BUILD_INPUTS) $(wildcard ../../problems/igrm_L2Projection/*.F90)
	$(MAKE) -j1 -C $(MYMAKE_DIR) CONFIG='$(CONFIG_PATH)' $(MYMAKE_COMPILER) igrm_l2

$(EXEC_DIR)/igrm_heat: $(BUILD_INPUTS) $(wildcard ../../problems/igrm_heat/*.F90)
	$(MAKE) -j1 -C $(MYMAKE_DIR) CONFIG='$(CONFIG_PATH)' $(MYMAKE_COMPILER) igrm_heat

$(EXEC_DIR)/igrm_eirksson: $(BUILD_INPUTS) $(wildcard ../../problems/igrm_eirksson/*.F90)
	$(MAKE) -j1 -C $(MYMAKE_DIR) CONFIG='$(CONFIG_PATH)' $(MYMAKE_COMPILER) igrm_eirksson

$(EXEC_DIR)/igrm_stokes: $(BUILD_INPUTS) $(wildcard ../../problems/igrm_stokes/*.F90)
	$(MAKE) -j1 -C $(MYMAKE_DIR) CONFIG='$(CONFIG_PATH)' $(MYMAKE_COMPILER) igrm_stokes

$(EXEC_DIR)/igrm_pollution: $(BUILD_INPUTS) $(wildcard ../../problems/igrm_pollution/*.F90)
	$(MAKE) -j1 -C $(MYMAKE_DIR) CONFIG='$(CONFIG_PATH)' $(MYMAKE_COMPILER) igrm_pollution

all: $(DRIVERS)

tests: all

run-cli: all
	MPIEXEC='$(MPIEXEC)' EXEC_DIR='$(abspath $(EXEC_DIR))' \
		DRIVER_CLI_TIMEOUT='$(DRIVER_CLI_TIMEOUT)' \
		SKIP_MPI_CASES='$(SKIP_MPI_CASES)' bash ./test_driver_cli.sh

run-smoke: all
	MPIEXEC='$(MPIEXEC)' EXEC_DIR='$(abspath $(EXEC_DIR))' \
		DRIVER_SMOKE_TIMEOUT='$(DRIVER_SMOKE_TIMEOUT)' \
		bash ./test_driver_smoke.sh

run-integration: all
	MPIEXEC='$(MPIEXEC)' EXEC_DIR='$(abspath $(EXEC_DIR))' \
		DRIVER_INTEGRATION_TIMEOUT='$(DRIVER_INTEGRATION_TIMEOUT)' \
		$(PYTHON) ./test_driver_integration.py

run: run-cli run-smoke run-integration

clean:
	@:
