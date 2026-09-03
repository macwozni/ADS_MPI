# Compatibility entry point for historical `make -C mymake ...` commands.
# Official build knowledge belongs to src/ and problems/.  Only an explicitly
# supplied SOURCE_ALL selects the isolated legacy adapter.
.DEFAULT_GOAL := all
.NOTPARALLEL:

MYMAKE_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
ROOT_DIR := $(abspath $(MYMAKE_DIR)/..)
SRC_DIR := $(ROOT_DIR)/src
PROBLEMS_DIR := $(ROOT_DIR)/problems

CONFIG ?= $(ROOT_DIR)/m_options
CONFIG_PATH := $(abspath $(CONFIG))
include $(CONFIG_PATH)

# A configuration may deliberately select another artifact tree, so resolve
# this only after loading CONFIG.
BUILD_ROOT ?= $(MYMAKE_DIR)
BUILD_ROOT_PATH := $(abspath $(BUILD_ROOT))

include $(PROBLEMS_DIR)/problems.mk
include $(MYMAKE_DIR)/legacy-source-build.mk

EXEC ?= l2
DELEGATE_OPTIONS = CONFIG="$(CONFIG_PATH)" BUILD_ROOT="$(BUILD_ROOT_PATH)"

.PHONY: all build library problems list-problems show-config clean clean-all \
	clean-problems clean-library clean-library-build help validate-exec \
	$(PROBLEMS)

validate-exec:
	@case ' $(PROBLEMS) ' in \
		*' $(EXEC) '*) ;; \
		*) printf 'Unknown EXEC: %s\nAvailable: %s\n' \
			'$(EXEC)' '$(PROBLEMS)' >&2; exit 2 ;; \
	esac

ifeq ($(origin SOURCE_ALL),undefined)
# Historical default: build the named EXEC (l2 unless overridden) through its
# owning problem makefile.
all build: validate-exec
	+$(MAKE) --no-print-directory -j1 -C $(PROBLEMS_DIR) \
		$(DELEGATE_OPTIONS) build-$(EXEC)

# Without SOURCE_ALL, clean retains the historical selected-EXEC/core-object
# scope, split between the two directories that now own those artifacts.
clean: validate-exec
	+$(MAKE) --no-print-directory -j1 -C $(PROBLEMS_DIR) \
		$(DELEGATE_OPTIONS) clean-$(EXEC)
	+$(MAKE) --no-print-directory -j1 -C $(SRC_DIR) \
		$(DELEGATE_OPTIONS) clean-objects
else
# Arbitrary legacy builds remain available, but they never populate or alter
# the official src/problem object directories.
all build: legacy-build

clean: legacy-clean
endif

# These public targets always delegate, even when SOURCE_ALL was inherited or
# supplied on the command line.
library:
	+$(MAKE) --no-print-directory -j1 -C $(SRC_DIR) \
		$(DELEGATE_OPTIONS) library

problems:
	+$(MAKE) --no-print-directory -j1 -C $(PROBLEMS_DIR) \
		$(DELEGATE_OPTIONS) problems

define DEFINE_PROBLEM_TARGET
$(1):
	+$$(MAKE) --no-print-directory -j1 -C $$(PROBLEMS_DIR) \
		$$(DELEGATE_OPTIONS) build-$(1)
endef
$(foreach problem,$(PROBLEMS),$(eval $(call DEFINE_PROBLEM_TARGET,$(problem))))

list-problems:
	+@$(MAKE) --no-print-directory -s -C $(PROBLEMS_DIR) \
		$(DELEGATE_OPTIONS) list-problems

show-config:
	+@$(MAKE) --no-print-directory -s -C $(SRC_DIR) \
		$(DELEGATE_OPTIONS) show-config

clean-problems:
	+$(MAKE) --no-print-directory -j1 -C $(PROBLEMS_DIR) \
		$(DELEGATE_OPTIONS) clean

clean-library:
	+$(MAKE) --no-print-directory -j1 -C $(SRC_DIR) \
		$(DELEGATE_OPTIONS) clean-library

# Historical alias: remove the complete core-library build, not configuration.
clean-library-build:
	+$(MAKE) --no-print-directory -j1 -C $(SRC_DIR) \
		$(DELEGATE_OPTIONS) clean

# Also remove every private legacy object tree and the executables registered
# by those trees.  The adapter verifies ownership markers before recursion.
clean-all: clean-problems clean-library-build legacy-clean-all

help:
	@printf '%s\n' \
		'mymake compatibility delegator' \
		'' \
		'  make                         build l2 (or a named EXEC override)' \
		'  make library                 delegate to ../src' \
		'  make problems                delegate all drivers to ../problems' \
		'  make <problem>               build one named driver' \
		'  make list-problems | show-config' \
		'  make clean                   clean the named EXEC and core objects' \
		'  make clean-problems|clean-library|clean-library-build|clean-all' \
		'' \
		'Legacy custom build:' \
		'  make SOURCE_ALL="../src/a.F90 ../app/main.F90" EXEC=my_app' \
		'  make clean SOURCE_ALL="..." EXEC=my_app' \
		'  Objects/modules: $(BUILD_ROOT_PATH)/_LEGACY_OBJ/<EXEC>' \
		'  Executable:      $(BUILD_ROOT_PATH)/EXEC/<EXEC>' \
		'' \
		'For new named builds, prefer the repository-root Makefile.' \
		'Configuration: $(CONFIG_PATH)' \
		'Build root:    $(BUILD_ROOT_PATH)'
