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
ORACLE_OBJ_DIR := $(BUILD_ROOT_PATH)/driver_cli_test_OBJ
ORACLE_SOURCE := pure_diffusion_nonzero_oracle.F90
ORACLE_OBJECT := $(ORACLE_OBJ_DIR)/pure_diffusion_nonzero_oracle.o
ORACLE_MODULE := $(ORACLE_OBJ_DIR)/pure_diffusion_oracle_rhs.mod
ORACLE_EXECUTABLE := $(EXEC_DIR)/pure_diffusion_nonzero_oracle
PROBLEM_BUILD_OPTIONS := CONFIG='$(CONFIG_PATH)' \
	                         BUILD_ROOT='$(BUILD_ROOT_PATH)' \
	                         COMPILER='$(MPIFC)' MPIFC='$(MPIFC)'

MODULE_OUTPUT ?= -J.
ifneq ($(filter -J.%,$(MODULE_OUTPUT)),)
ORACLE_MODULE_FLAGS := -J$(ORACLE_OBJ_DIR)
else ifneq ($(filter -module,$(MODULE_OUTPUT)),)
ORACLE_MODULE_FLAGS := -module $(ORACLE_OBJ_DIR)
else
ORACLE_MODULE_FLAGS := $(MODULE_OUTPUT)
endif
ORACLE_FFLAGS := $(DIALECT) $(BUILD_OPTS) $(INCLUDE) \
	$(ORACLE_MODULE_FLAGS) -I$(ORACLE_OBJ_DIR) -I$(BUILD_ROOT_PATH)

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

$(ORACLE_OBJ_DIR) $(EXEC_DIR):
	mkdir -p $@

# Depending on one real driver ensures that the core library and its module
# files were rebuilt with the exact compiler/configuration used by this suite.
$(ORACLE_OBJECT): $(ORACLE_SOURCE) $(firstword $(DRIVERS)) GNUmakefile \
		$(CONFIG_PATH) | $(ORACLE_OBJ_DIR)
	$(MPIFC) $(ORACLE_FFLAGS) -c $< -o $@

$(ORACLE_EXECUTABLE): $(ORACLE_OBJECT) | $(EXEC_DIR)
	$(MPIFC) $(ORACLE_FFLAGS) -o $@ $< \
		$(BUILD_ROOT_PATH)/LIB/libads.a $(USER_LIB)

all: $(DRIVERS) $(ORACLE_EXECUTABLE)

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
	$(RM) -- $(ORACLE_OBJECT) $(ORACLE_MODULE) $(ORACLE_EXECUTABLE)
	@rmdir -- '$(ORACLE_OBJ_DIR)' 2>/dev/null || :
