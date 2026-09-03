# Shared implementation for a single problem driver.
#
# Before including this file, a problem GNUmakefile must define EXEC_NAME,
# SOURCES, PROBLEM_MODULES, DEFAULT_ARGS, and USAGE. Source order is dependency
# order; PROBLEM_MODULES lists the module files owned by this driver.

ifndef EXEC_NAME
$(error EXEC_NAME must be defined before including problems/problem.mk)
endif
ifndef SOURCES
$(error SOURCES must be defined before including problems/problem.mk)
endif
ifndef PROBLEM_MODULES
$(error PROBLEM_MODULES must be defined before including problems/problem.mk)
endif
ifndef DEFAULT_ARGS
$(error DEFAULT_ARGS must be defined before including problems/problem.mk)
endif
ifndef USAGE
$(error USAGE must be defined before including problems/problem.mk)
endif

.DEFAULT_GOAL := all
.NOTPARALLEL:

PROBLEMS_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
ROOT_DIR := $(abspath $(PROBLEMS_DIR)/..)
PROBLEM_DIR := $(patsubst %/,%,$(dir $(abspath $(firstword $(MAKEFILE_LIST)))))

CONFIG ?= $(ROOT_DIR)/m_options
CONFIG_PATH := $(abspath $(CONFIG))
include $(CONFIG_PATH)

SRC_DIR := $(ROOT_DIR)/src
BUILD_ROOT ?= $(ROOT_DIR)/mymake
BUILD_ROOT_PATH := $(abspath $(BUILD_ROOT))
CORE_MODULE_DIR := $(BUILD_ROOT_PATH)
LIBRARY := $(CORE_MODULE_DIR)/LIB/libads.a
EXEC_DIR ?= $(BUILD_ROOT_PATH)/EXEC
EXEC_DIR_PATH := $(abspath $(EXEC_DIR))
PROBLEM_OBJ_DIR ?= _OBJ
PROBLEM_OBJ_DIR_PATH := $(abspath $(PROBLEM_OBJ_DIR))
EXECUTABLE := $(EXEC_DIR_PATH)/$(EXEC_NAME)
BUILD_STAMP := $(PROBLEM_OBJ_DIR_PATH)/.ads-problem-build

OBJECTS := $(addprefix $(PROBLEM_OBJ_DIR_PATH)/,$(SOURCES:.F90=.o))
MODULE_FILES := $(addprefix $(PROBLEM_OBJ_DIR_PATH)/,$(addsuffix .mod,$(PROBLEM_MODULES)))

# Redirect the compiler-specific module output spelling from the configuration
# to this problem's private object directory.  The explicit include paths make
# both core and earlier problem modules available to subsequent sources.
MODULE_OUTPUT ?= -J.
ifneq ($(filter -J.%,$(MODULE_OUTPUT)),)
MODULE_FLAGS := -J$(PROBLEM_OBJ_DIR_PATH)
else ifneq ($(filter -module,$(MODULE_OUTPUT)),)
MODULE_FLAGS := -module $(PROBLEM_OBJ_DIR_PATH)
else
MODULE_FLAGS := $(MODULE_OUTPUT)
endif

override FFLAGS += $(MODULE_FLAGS) -I$(PROBLEM_OBJ_DIR_PATH) -I$(CORE_MODULE_DIR)

NP ?= 1
MPIEXEC_FLAGS ?=
MPI_NP_FLAG ?= -n
OMP_NUM_THREADS ?= 1
OMP_DYNAMIC ?= FALSE
OMP_PROC_BIND ?= close
RUN_ENV ?=
RUN_DIR ?=
ARGS ?=
OIL_SEED ?=
PROBLEM_RUN_ENV ?=

EFFECTIVE_ARGS = $(if $(strip $(ARGS)),$(ARGS),$(if $(strip $($(EXEC_NAME)_ARGS)),$($(EXEC_NAME)_ARGS),$(DEFAULT_ARGS)))
EFFECTIVE_RUN_DIR = $(if $(strip $(RUN_DIR)),$(if $(filter /%,$(RUN_DIR)),$(RUN_DIR),$(ROOT_DIR)/$(RUN_DIR)),$(ROOT_DIR)/output/$(EXEC_NAME))

.PHONY: all build library run show-run show-usage run-help list-sources show-config \
	clean clean-objects clean-executable distclean help validate-problem-paths FORCE

all build: validate-problem-paths library $(EXECUTABLE)

validate-problem-paths:
	@case '$(BUILD_ROOT_PATH)' in \
		'/'|'$(ROOT_DIR)'|'$(PROBLEMS_DIR)'|'$(PROBLEM_DIR)') \
			printf 'Unsafe BUILD_ROOT: %s\n' '$(BUILD_ROOT_PATH)' >&2; exit 2 ;; \
	esac
	@case '$(PROBLEM_OBJ_DIR_PATH)' in \
		'/'|'$(ROOT_DIR)'|'$(SRC_DIR)'|'$(PROBLEMS_DIR)'|'$(PROBLEM_DIR)'|\
		'$(BUILD_ROOT_PATH)'|'$(EXEC_DIR_PATH)') \
			printf 'Unsafe PROBLEM_OBJ_DIR: %s\n' '$(PROBLEM_OBJ_DIR_PATH)' >&2; exit 2 ;; \
	esac
	@case '$(EXEC_DIR_PATH)' in \
		'/'|'$(ROOT_DIR)'|'$(SRC_DIR)'|'$(PROBLEMS_DIR)'|'$(PROBLEM_DIR)') \
			printf 'Unsafe EXEC_DIR: %s\n' '$(EXEC_DIR_PATH)' >&2; exit 2 ;; \
	esac

# Always ask src/ to update the core archive first.  Its own dependency graph
# makes an unchanged invocation cheap, while a changed module is rebuilt before
# any problem source is compiled.
library: validate-problem-paths
	+$(MAKE) --no-print-directory -j1 -C $(SRC_DIR) \
		CONFIG="$(CONFIG_PATH)" BUILD_ROOT="$(BUILD_ROOT_PATH)" library

# This rule also supports a direct request for the executable when the archive
# has not been built yet.  Normal builds reach it through the phony target above.
$(LIBRARY):
	+$(MAKE) --no-print-directory -j1 -C $(SRC_DIR) \
		CONFIG="$(CONFIG_PATH)" BUILD_ROOT="$(BUILD_ROOT_PATH)" library

$(EXECUTABLE): $(OBJECTS) $(LIBRARY) | $(EXEC_DIR_PATH)
	$(FF) $(FFLAGS) -o $@ $(OBJECTS) $(LIBRARY) $(USER_LIB)
	@printf 'Built %s\n' '$@'

$(PROBLEM_OBJ_DIR_PATH)/%.o: %.F90 $(MAKEFILE_LIST) | $(PROBLEM_OBJ_DIR_PATH)
	$(FF) $(FFLAGS) -o $@ -c $<

# A command-line BUILD/FFLAGS change does not alter a makefile timestamp.
# Persist the effective compilation and link settings so a profile/toolchain
# switch invalidates this problem's complete ordered object chain.
FORCE:

$(BUILD_STAMP): FORCE | $(PROBLEM_OBJ_DIR_PATH)
	@build_stamp_tmp='$@.tmp'; \
	missing_module=0; \
	for module_file in $(MODULE_FILES); do \
		test -r "$$module_file" || missing_module=1; \
	done; \
	{ \
		printf '%s\n' 'EXEC_NAME=$(EXEC_NAME)'; \
		printf '%s\n' 'FF=$(FF)'; \
		printf '%s\n' 'FFLAGS=$(FFLAGS)'; \
		printf '%s\n' 'USER_LIB=$(USER_LIB)'; \
	} > "$$build_stamp_tmp"; \
	if test "$$missing_module" -eq 0 && test -r '$@' && \
		cmp -s "$$build_stamp_tmp" '$@'; then \
		$(RM) -- "$$build_stamp_tmp"; \
	else \
		mv -f -- "$$build_stamp_tmp" '$@'; \
	fi

$(PROBLEM_OBJ_DIR_PATH) $(EXEC_DIR_PATH):
	mkdir -p $@

# Turn the local source order into real prerequisites.  This makes `make -j`
# safe and recompiles all later consumers after a module interface changes.
define add-ordered-dependencies
$(if $(word 2,$1),$(eval $(word 2,$1): $(word 1,$1))$(call add-ordered-dependencies,$(wordlist 2,$(words $1),$1)))
endef
$(call add-ordered-dependencies,$(OBJECTS))
$(firstword $(OBJECTS)): $(BUILD_STAMP) $(LIBRARY)

run: build
	@mkdir -p "$(EFFECTIVE_RUN_DIR)"
	@printf 'Running %s in %s with %s MPI rank(s), %s OpenMP thread(s)\n' \
		'$(EXEC_NAME)' '$(EFFECTIVE_RUN_DIR)' '$(NP)' '$(OMP_NUM_THREADS)'
	@cd "$(EFFECTIVE_RUN_DIR)" && \
	env OMP_DYNAMIC='$(OMP_DYNAMIC)' \
		OMP_NUM_THREADS='$(OMP_NUM_THREADS)' \
		OMP_PROC_BIND='$(OMP_PROC_BIND)' \
		$(PROBLEM_RUN_ENV) $(RUN_ENV) \
		"$(MPIEXEC)" $(MPIEXEC_FLAGS) "$(MPI_NP_FLAG)" "$(NP)" \
		"$(EXECUTABLE)" $(EFFECTIVE_ARGS)

show-run:
	@printf '%-20s %s\n' \
		'PROBLEM' '$(EXEC_NAME)' \
		'EXECUTABLE' '$(EXECUTABLE)' \
		'USAGE' '$(USAGE)' \
		'ARGS' '$(EFFECTIVE_ARGS)' \
		'NP' '$(NP)' \
		'MPIEXEC' '$(MPIEXEC)' \
		'MPIEXEC_FLAGS' '$(MPIEXEC_FLAGS)' \
		'MPI_NP_FLAG' '$(MPI_NP_FLAG)' \
		'OMP_NUM_THREADS' '$(OMP_NUM_THREADS)' \
		'OMP_DYNAMIC' '$(OMP_DYNAMIC)' \
		'OMP_PROC_BIND' '$(OMP_PROC_BIND)' \
		'RUN_DIR' '$(EFFECTIVE_RUN_DIR)' \
		'RUN_ENV' '$(RUN_ENV)' \
		'OIL_SEED' '$(OIL_SEED)' \
		'PROBLEM_RUN_ENV' '$(PROBLEM_RUN_ENV)'

show-usage run-help:
	@printf '%s:\n  %s\n  default: %s\n' \
		'$(EXEC_NAME)' '$(USAGE)' '$(EFFECTIVE_ARGS)'

list-sources:
	@printf '%s\n' $(SOURCES)

show-config:
	@printf '%-16s %s\n' \
		'CONFIG' '$(CONFIG_PATH)' \
		'BUILD_ROOT' '$(BUILD_ROOT_PATH)' \
		'BUILD' '$(BUILD)' \
		'COMPILER' '$(COMPILER)' \
		'OBJECTS' '$(PROBLEM_OBJ_DIR_PATH)' \
		'MODULES' '$(PROBLEM_OBJ_DIR_PATH)' \
		'LIBRARY' '$(LIBRARY)' \
		'EXECUTABLE' '$(EXECUTABLE)'

clean: clean-executable clean-objects

clean-objects: validate-problem-paths
	$(RM) -- $(OBJECTS) $(MODULE_FILES) $(BUILD_STAMP) $(BUILD_STAMP).tmp
	@rmdir -- '$(PROBLEM_OBJ_DIR_PATH)' 2>/dev/null || :

clean-executable: validate-problem-paths
	$(RM) -- $(EXECUTABLE)
	@rmdir -- '$(EXEC_DIR_PATH)' 2>/dev/null || :

# Configuration and the shared src library are not owned by a problem and are
# deliberately preserved by both cleanup targets.
distclean: clean

help:
	@printf '%s\n' \
		'Problem: $(EXEC_NAME)' \
		'  make | make build       update src library and build $(EXECUTABLE)' \
		'  make run                build and run with DEFAULT_ARGS' \
		'  make show-run           print effective runtime configuration' \
		'  make show-usage         print argument syntax and defaults' \
		'  make list-sources       print this problem source list' \
		'  make show-config        print build configuration and paths' \
		'  make clean              remove only this problem artifacts' \
		'Runtime overrides: ARGS NP MPIEXEC MPIEXEC_FLAGS MPI_NP_FLAG RUN_DIR' \
		'                   RUN_ENV OMP_NUM_THREADS OMP_DYNAMIC OMP_PROC_BIND' \
		'                   OIL_SEED (used only by oil)' \
		'Configuration: CONFIG=$(CONFIG_PATH)'
