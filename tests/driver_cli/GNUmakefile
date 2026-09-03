.DEFAULT_GOAL := all

.NOTPARALLEL:

.PHONY: all tests run run-cli run-smoke run-integration clean

CONFIG     ?= ../../m_options
CONFIG_PATH := $(abspath $(CONFIG))
include $(CONFIG_PATH)

MYMAKE_DIR ?= ../../mymake
SRC_DIR    ?= ../../src
PROBLEMS_DIR ?= ../../problems
MYMAKE_PATH := $(abspath $(MYMAKE_DIR))
SRC_PATH := $(abspath $(SRC_DIR))
PROBLEMS_PATH := $(abspath $(PROBLEMS_DIR))
BUILD_ROOT ?= $(MYMAKE_PATH)
BUILD_ROOT_PATH := $(abspath $(BUILD_ROOT))
include $(PROBLEMS_PATH)/problems.mk

MPIEXEC    ?= /opt/lib/mpich-5.0.0/bin/mpiexec
MPIFC      ?= /opt/lib/mpich-5.0.0/bin/mpif90
PYTHON     ?= python3
DRIVER_CLI_TIMEOUT ?= 20s
DRIVER_SMOKE_TIMEOUT ?= 60s
DRIVER_INTEGRATION_TIMEOUT ?= 90
SKIP_MPI_CASES ?= 0
EXEC_DIR   := $(BUILD_ROOT_PATH)/EXEC
PROBLEM_BUILD_OPTIONS := CONFIG='$(CONFIG_PATH)' \
	                         BUILD_ROOT='$(BUILD_ROOT_PATH)' \
	                         COMPILER='$(MPIFC)' MPIFC='$(MPIFC)'

COMMON_SOURCES := $(wildcard $(SRC_PATH)/*.F90)
BUILD_INPUTS   := GNUmakefile $(COMMON_SOURCES) $(SRC_PATH)/GNUmakefile \
	                  $(SRC_PATH)/sources.mk $(PROBLEMS_PATH)/GNUmakefile \
	                  $(PROBLEMS_PATH)/problem.mk $(PROBLEMS_PATH)/problems.mk \
	                  $(CONFIG_PATH)

DRIVERS := $(addprefix $(EXEC_DIR)/,$(PROBLEMS))

# A driver test is also a clean-build test.  Every executable is rebuilt via
# its owning problem makefile, with the test compiler/configuration forwarded.
.PHONY: $(DRIVERS)

define DEFINE_DRIVER_BUILD
$(EXEC_DIR)/$(1): $(BUILD_INPUTS) \
		$(wildcard $(PROBLEMS_PATH)/$($(1)_DIR)/*.F90) \
		$(PROBLEMS_PATH)/$($(1)_DIR)/GNUmakefile
	+$$(MAKE) --no-print-directory -j1 -C $$(PROBLEMS_PATH) \
		$$(PROBLEM_BUILD_OPTIONS) clean-$(1)
	+$$(MAKE) --no-print-directory -j1 -C $$(PROBLEMS_PATH) \
		$$(PROBLEM_BUILD_OPTIONS) build-$(1)
endef

$(foreach problem,$(PROBLEMS),$(eval $(call DEFINE_DRIVER_BUILD,$(problem))))

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
