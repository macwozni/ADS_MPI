.DEFAULT_GOAL := all

# The core library owns its source list and build rules.  Generated files stay
# in the historical mymake locations so existing problem and external build
# commands keep working.
SRC_DIR := $(patsubst %/,%,$(dir $(abspath $(firstword $(MAKEFILE_LIST)))))
ROOT_DIR := $(abspath $(SRC_DIR)/..)

CONFIG ?= $(ROOT_DIR)/m_options
CONFIG_PATH := $(abspath $(CONFIG))
include $(CONFIG_PATH)
include $(SRC_DIR)/sources.mk

vpath %.F90 $(SRC_DIR)

.NOTPARALLEL:

BUILD_ROOT ?= $(ROOT_DIR)/mymake
BUILD_ROOT_PATH := $(abspath $(BUILD_ROOT))
CORE_OBJ_DIR ?= $(BUILD_ROOT_PATH)/_OBJ
CORE_OBJ_DIR_PATH := $(abspath $(CORE_OBJ_DIR))
CORE_LIB_DIR := $(BUILD_ROOT_PATH)/LIB
CORE_MODULE_DIR := $(BUILD_ROOT_PATH)

LIBRARY := $(CORE_LIB_DIR)/libads.a
BUILD_STAMP := $(CORE_OBJ_DIR_PATH)/.ads-core-build

OBJECTS := $(addprefix $(CORE_OBJ_DIR_PATH)/,$(SRC_FILES:.F90=.o))
MODULE_FILES := $(addprefix $(CORE_MODULE_DIR)/,$(addsuffix .mod,$(SRC_MODULES)))

# m_options describes the compiler-specific spelling using the old current
# directory ("-J." for GNU, "-module ." for Intel).  Redirect that output to
# the shared compatibility directory without hard-coding one compiler family.
MODULE_OUTPUT ?= -J.
ifneq ($(filter -J.%,$(MODULE_OUTPUT)),)
MODULE_FLAGS := -J$(CORE_MODULE_DIR)
else ifneq ($(filter -module,$(MODULE_OUTPUT)),)
MODULE_FLAGS := -module $(CORE_MODULE_DIR)
else
MODULE_FLAGS := $(MODULE_OUTPUT)
endif

override FFLAGS += $(MODULE_FLAGS) -I$(CORE_MODULE_DIR) -I$(CORE_OBJ_DIR_PATH)
ARFLAGS = rcs

.PHONY: all build library list-sources show-config clean clean-objects \
	clean-library distclean help validate-core-paths FORCE

all build library: validate-core-paths $(LIBRARY)

validate-core-paths:
	@case '$(BUILD_ROOT_PATH)' in \
		'/'|'$(ROOT_DIR)'|'$(SRC_DIR)') \
			printf 'Unsafe BUILD_ROOT: %s\n' '$(BUILD_ROOT_PATH)' >&2; exit 2 ;; \
	esac
	@case '$(CORE_OBJ_DIR_PATH)' in \
		'/'|'$(ROOT_DIR)'|'$(SRC_DIR)'|'$(BUILD_ROOT_PATH)') \
			printf 'Unsafe CORE_OBJ_DIR: %s\n' '$(CORE_OBJ_DIR_PATH)' >&2; exit 2 ;; \
	esac

$(LIBRARY): $(OBJECTS) | $(CORE_LIB_DIR)
	$(RM) -- $@
	$(AR) $(ARFLAGS) $@ $(OBJECTS)
	@printf 'Built %s\n' '$@'

$(CORE_OBJ_DIR_PATH)/%.o: %.F90 $(MAKEFILE_LIST) | $(CORE_OBJ_DIR_PATH) $(CORE_MODULE_DIR)
	$(FF) $(FFLAGS) -o $@ -c $<

# Command-line BUILD/FFLAGS and a newly selected CONFIG do not necessarily
# change a makefile timestamp.  Record the effective toolchain so switching
# configuration cannot silently reuse ABI-incompatible Fortran objects.
FORCE:

$(BUILD_STAMP): FORCE | $(CORE_OBJ_DIR_PATH)
	@build_stamp_tmp='$@.tmp'; \
	missing_module=0; \
	for module_file in $(MODULE_FILES); do \
		test -r "$$module_file" || missing_module=1; \
	done; \
	{ \
		printf '%s\n' 'FF=$(FF)'; \
		printf '%s\n' 'FFLAGS=$(FFLAGS)'; \
		printf '%s\n' 'AR=$(AR)'; \
		printf '%s\n' 'ARFLAGS=$(ARFLAGS)'; \
	} > "$$build_stamp_tmp"; \
	if test "$$missing_module" -eq 0 && test -r '$@' && \
		cmp -s "$$build_stamp_tmp" '$@'; then \
		$(RM) -- "$$build_stamp_tmp"; \
	else \
		mv -f -- "$$build_stamp_tmp" '$@'; \
	fi

$(CORE_OBJ_DIR_PATH) $(CORE_LIB_DIR) $(CORE_MODULE_DIR):
	mkdir -p $@

# Turn the documented source order into real prerequisites.  Besides making
# parallel invocations safe, this also recompiles all later module consumers
# after an earlier interface changes.
define add-ordered-dependencies
$(if $(word 2,$1),$(eval $(word 2,$1): $(word 1,$1))$(call add-ordered-dependencies,$(wordlist 2,$(words $1),$1)))
endef
$(call add-ordered-dependencies,$(OBJECTS))
$(firstword $(OBJECTS)): $(BUILD_STAMP)

list-sources:
	@printf '%s\n' $(SRC_FILES)

show-config:
	@printf '%-16s %s\n' \
		'CONFIG' '$(CONFIG_PATH)' \
		'BUILD_ROOT' '$(BUILD_ROOT_PATH)' \
		'BUILD' '$(BUILD)' \
		'COMPILER' '$(COMPILER)' \
		'OBJECTS' '$(CORE_OBJ_DIR_PATH)' \
		'MODULES' '$(CORE_MODULE_DIR)' \
		'LIBRARY' '$(LIBRARY)'

clean: clean-library clean-objects

clean-objects: validate-core-paths
	$(RM) -- $(OBJECTS) $(BUILD_STAMP) $(BUILD_STAMP).tmp $(MODULE_FILES)
	@rmdir -- '$(CORE_OBJ_DIR_PATH)' 2>/dev/null || :

clean-library: validate-core-paths
	$(RM) -- $(LIBRARY)
	@rmdir -- '$(CORE_LIB_DIR)' 2>/dev/null || :

# Configuration is user-owned and intentionally survives cleanup.
distclean: clean

help:
	@printf '%s\n' \
		'make                    build the ADS core library' \
		'make library            build mymake/LIB/libads.a' \
		'make list-sources       print the ordered core source list' \
		'make show-config        print the effective configuration and paths' \
		'make clean              remove core objects, modules, and library' \
		'make clean-library      remove only libads.a' \
		'make clean-objects      remove only core objects and modules' \
		'make CONFIG=...         select another root-compatible configuration'
