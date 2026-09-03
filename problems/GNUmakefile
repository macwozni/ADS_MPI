.DEFAULT_GOAL := all
.NOTPARALLEL:

PROBLEMS_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
ROOT_DIR := $(abspath $(PROBLEMS_DIR)/..)

include $(PROBLEMS_DIR)/problems.mk

CONFIG ?= $(ROOT_DIR)/m_options
CONFIG_PATH := $(abspath $(CONFIG))
include $(CONFIG_PATH)
BUILD_ROOT ?= $(ROOT_DIR)/mymake
BUILD_ROOT_PATH := $(abspath $(BUILD_ROOT))
PROBLEM ?= l2

DELEGATE_OPTIONS = CONFIG="$(CONFIG_PATH)" BUILD_ROOT="$(BUILD_ROOT_PATH)"

.PHONY: all problems build run show-run show-usage run-help list \
	list-problems list-directories list-sources show-config clean \
	clean-problems distclean help validate-problem

all problems: $(addprefix build-,$(PROBLEMS))

validate-problem:
	@case ' $(PROBLEMS) ' in \
		*' $(PROBLEM) '*) ;; \
		*) printf 'Unknown PROBLEM: %s\nAvailable: %s\n' \
			'$(PROBLEM)' '$(PROBLEMS)' >&2; exit 2 ;; \
	esac

build: validate-problem build-$(PROBLEM)
run: validate-problem run-$(PROBLEM)
show-run: validate-problem show-run-$(PROBLEM)
show-usage: validate-problem usage-$(PROBLEM)
list-sources: validate-problem sources-$(PROBLEM)
show-config: validate-problem config-$(PROBLEM)

define DEFINE_PROBLEM_TARGETS
.PHONY: $(1) build-$(1) run-$(1) show-run-$(1) usage-$(1) \
	sources-$(1) config-$(1) clean-$(1)

$(1) build-$(1):
	+$$(MAKE) --no-print-directory -j1 -C "$$(PROBLEMS_DIR)/$$($(1)_DIR)" \
		$$(DELEGATE_OPTIONS) build

run-$(1):
	+$$(MAKE) --no-print-directory -j1 -C "$$(PROBLEMS_DIR)/$$($(1)_DIR)" \
		$$(DELEGATE_OPTIONS) run

show-run-$(1):
	+@$$(MAKE) --no-print-directory -j1 -C "$$(PROBLEMS_DIR)/$$($(1)_DIR)" \
		$$(DELEGATE_OPTIONS) show-run

usage-$(1):
	+@$$(MAKE) --no-print-directory -j1 -C "$$(PROBLEMS_DIR)/$$($(1)_DIR)" \
		$$(DELEGATE_OPTIONS) show-usage

sources-$(1):
	+@$$(MAKE) --no-print-directory -j1 -C "$$(PROBLEMS_DIR)/$$($(1)_DIR)" \
		$$(DELEGATE_OPTIONS) list-sources

config-$(1):
	+@$$(MAKE) --no-print-directory -j1 -C "$$(PROBLEMS_DIR)/$$($(1)_DIR)" \
		$$(DELEGATE_OPTIONS) show-config

clean-$(1):
	+$$(MAKE) --no-print-directory -j1 -C "$$(PROBLEMS_DIR)/$$($(1)_DIR)" \
		$$(DELEGATE_OPTIONS) clean
endef
$(foreach problem,$(PROBLEMS),$(eval $(call DEFINE_PROBLEM_TARGETS,$(problem))))

run-help: $(addprefix usage-,$(PROBLEMS))

list list-problems:
	@printf '%s\n' $(PROBLEMS)

list-directories:
	@$(foreach problem,$(PROBLEMS),printf '%-24s %s\n' \
		'$(problem)' '$($(problem)_DIR)';)

clean clean-problems: $(addprefix clean-,$(PROBLEMS))

# The core library and root configuration belong to src/ and are not removed
# by the problems subtree.
distclean: clean

help:
	@printf '%s\n' \
		'ADS MPI problem drivers' \
		'' \
		'  make | make problems                  build all problems' \
		'  make build PROBLEM=heat              build one problem' \
		'  make build-heat | make heat           named build target' \
		'  make run PROBLEM=heat [ARGS="..."]    build and run one problem' \
		'  make run-heat [ARGS="..."]            named run target' \
		'  make show-run PROBLEM=heat            print effective run settings' \
		'  make list-problems | list-directories' \
		'  make list-sources PROBLEM=heat' \
		'  make run-help                         list all argument syntaxes' \
		'  make clean | make clean-heat          clean all or one problem' \
		'' \
		'Runtime overrides: ARGS NP MPIEXEC MPIEXEC_FLAGS MPI_NP_FLAG RUN_DIR' \
		'                   RUN_ENV OMP_NUM_THREADS OMP_DYNAMIC OMP_PROC_BIND' \
		'                   OIL_SEED (used only by oil)' \
		'Configuration: CONFIG=$(CONFIG_PATH) BUILD_ROOT=$(BUILD_ROOT_PATH)'
