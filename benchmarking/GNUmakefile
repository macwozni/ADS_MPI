.DEFAULT_GOAL := help
.NOTPARALLEL:

BENCHMARKING_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
ROOT_DIR := $(abspath $(BENCHMARKING_DIR)/..)
CONFIG_DIR ?= $(BENCHMARKING_DIR)/configs
MPIEXEC_FLAGS ?=
BENCHMARK_LAUNCHER_TEMPLATE ?=
override BENCHMARK_LAUNCHER_TEMPLATE := $(value BENCHMARK_LAUNCHER_TEMPLATE)
export BENCHMARK_LAUNCHER_TEMPLATE

BENCHMARK_PROFILE ?= smoke
BENCHMARK_CONVERGENCE_PROFILE ?= temporal-full
BENCHMARK_H_CONVERGENCE_PROFILE ?= h-convergence-full
BENCHMARK_P_CONVERGENCE_PROFILE ?= p-convergence-full
BENCHMARK_VALIDATION_PROFILE ?= validation-full
BENCHMARK_STRONG_PROFILE ?= cluster-scaling
BENCHMARK_WEAK_PROFILE ?= cluster-weak-scaling
BENCHMARK_SHARD_PROFILE ?= $(BENCHMARK_WEAK_PROFILE)
BENCHMARK_SHARD_STRATEGY ?= problem-scheme-degree
BENCHMARK_SHARD_COUNT ?=
BENCHMARK_RESUME_PROFILE ?= $(BENCHMARK_CONVERGENCE_PROFILE)
BENCHMARK_PLAN_ARGS ?=
BENCHMARK_ANALYZE_ARGS ?=
BENCHMARK_RUN_ID ?=
RUN_ID ?=
PARENT_RUN_ID ?=
SHARD_INDEX ?=
SHARD_RUNS ?=
BENCHMARK_BUILD_ROOT ?= $(BENCHMARKING_DIR)/build
BENCHMARK_BUILD_ROOT_PATH := $(abspath $(BENCHMARK_BUILD_ROOT))
BENCHMARK_BUILD_MARKER := $(BENCHMARK_BUILD_ROOT_PATH)/.ads-benchmark-build

BENCHMARK_ADAPTERS := igrm_l2 igrm_heat pure_diffusion_igrm
BENCHMARK_CONFIG_GOALS := build build-adapters smoke convergence \
	h-convergence p-convergence validate strong weak run-shard resume \
	show-build-config \
	$(addprefix build-,$(BENCHMARK_ADAPTERS))
ifneq ($(filter plan,$(MAKECMDGOALS)),)
ifneq ($(findstring --show-commands,$(BENCHMARK_PLAN_ARGS)),)
ifneq ($(origin CONFIG),undefined)
BENCHMARK_CONFIG_GOALS += plan
endif
endif
endif
BENCHMARK_NEEDS_CONFIG := $(strip $(filter $(BENCHMARK_CONFIG_GOALS),$(MAKECMDGOALS)))

# Planning, self-tests, help, and cleanup must remain usable without a local
# Fortran toolchain configuration.  In particular, the cleanup self-test copies
# this makefile alone into a disposable repository.
ifneq ($(BENCHMARK_NEEDS_CONFIG),)
CONFIG ?= $(ROOT_DIR)/m_options
CONFIG_PATH := $(abspath $(CONFIG))

# The checked-in smoke profile selects the debug build.  Set this before the
# configuration is parsed because m_options chooses BUILD_OPTS in conditionals.
ifneq ($(filter smoke,$(MAKECMDGOALS)),)
override BUILD := debug
endif
ifneq ($(filter convergence h-convergence p-convergence validate strong weak run-shard resume,$(MAKECMDGOALS)),)
override BUILD := release
endif

include $(CONFIG_PATH)

MPI_NP_FLAG ?= -n
export MPIEXEC MPIEXEC_FLAGS MPI_NP_FLAG

ifeq ($(filter $(BUILD),debug release),)
$(error unsupported benchmark BUILD='$(BUILD)'; expected debug or release)
endif

FORTRAN_DIR := $(BENCHMARKING_DIR)/fortran
SRC_DIR := $(ROOT_DIR)/src
PROBLEMS_DIR := $(ROOT_DIR)/problems
BENCHMARK_PROFILE_BUILD_ROOT := $(BENCHMARK_BUILD_ROOT_PATH)/$(BUILD)
CORE_MODULE_DIR := $(BENCHMARK_PROFILE_BUILD_ROOT)
CORE_LIBRARY := $(BENCHMARK_PROFILE_BUILD_ROOT)/LIB/libads.a
BENCHMARK_OBJ_ROOT := $(BENCHMARK_PROFILE_BUILD_ROOT)/benchmark_OBJ
BENCHMARK_EXEC_DIR := $(BENCHMARK_PROFILE_BUILD_ROOT)/EXEC

BENCHMARK_COMMON_PREFIX_SOURCES := \
	$(FORTRAN_DIR)/benchmark_contract.F90 \
	$(FORTRAN_DIR)/manufactured_solution.F90 \
	$(FORTRAN_DIR)/benchmark_runtime.F90
BENCHMARK_COMMON_SUFFIX_SOURCES := \
	$(FORTRAN_DIR)/benchmark_cli.F90 \
	$(FORTRAN_DIR)/benchmark_harness.F90
IGRM_HEAT_SUPPORT_SOURCES := \
	$(PROBLEMS_DIR)/igrm_heat/input_data.F90 \
	$(PROBLEMS_DIR)/igrm_heat/RHS_fun.F90

benchmark_module_flags = $(if $(filter -J.%,$(MODULE_OUTPUT)),-J$(1),$(if $(filter -module,$(MODULE_OUTPUT)),-module $(1),$(MODULE_OUTPUT)))

.PHONY: library validate-build-profile FORCE FORCE_CORE_LIBRARY

FORCE:

FORCE_CORE_LIBRARY:

validate-build-profile: validate-build-ownership
	@set -eu; \
	root=$$(realpath -m -- '$(BENCHMARK_BUILD_ROOT_PATH)'); \
	profile=$$(realpath -m -- '$(BENCHMARK_PROFILE_BUILD_ROOT)'); \
	expected="$$root/$(BUILD)"; \
	if test "$$profile" != "$$expected"; then \
		printf 'Unsafe benchmark profile build root: %s\n' \
			'$(BENCHMARK_PROFILE_BUILD_ROOT)' >&2; exit 2; \
	fi; \
	for path in '$(BENCHMARK_PROFILE_BUILD_ROOT)' \
		'$(BENCHMARK_PROFILE_BUILD_ROOT)/_OBJ' \
		'$(BENCHMARK_PROFILE_BUILD_ROOT)/LIB' \
		'$(BENCHMARK_OBJ_ROOT)' '$(BENCHMARK_EXEC_DIR)' \
		$(foreach adapter,$(BENCHMARK_ADAPTERS),'$(BENCHMARK_OBJ_ROOT)/$(adapter)'); do \
		if test -L "$$path"; then \
			printf 'Unsafe benchmark build symlink: %s\n' "$$path" >&2; \
			exit 2; \
		fi; \
	done

library: $(CORE_LIBRARY)

# Refresh the isolated core exactly once per public build invocation.  The
# recursive owner keeps the archive timestamp unchanged when its inputs did
# not change, so normal downstream dependencies remain incremental.
$(CORE_LIBRARY): FORCE_CORE_LIBRARY validate-build-profile prepare-build-root
	+$(MAKE) --no-print-directory -j1 -C '$(SRC_DIR)' \
		CONFIG='$(CONFIG_PATH)' BUILD='$(BUILD)' \
		BUILD_ROOT='$(BENCHMARK_PROFILE_BUILD_ROOT)' library
	@test -f '$@' || { \
		printf 'Core benchmark library was not produced: %s\n' '$@' >&2; \
		exit 2; \
	}

$(BENCHMARK_EXEC_DIR):
	mkdir -p -- '$@'

define DEFINE_BENCHMARK_ADAPTER
$(1)_SOURCES := $$(BENCHMARK_COMMON_PREFIX_SOURCES) \
	$$(if $$(filter igrm_heat,$(1)),$$(IGRM_HEAT_SUPPORT_SOURCES)) \
	$$(FORTRAN_DIR)/$(1)_adapter.F90 \
	$$(BENCHMARK_COMMON_SUFFIX_SOURCES) \
	$$(FORTRAN_DIR)/main_$(1).F90
$(1)_OBJ_DIR := $$(BENCHMARK_OBJ_ROOT)/$(1)
$(1)_OBJECTS := $$(addprefix $$($(1)_OBJ_DIR)/,$$(notdir $$($(1)_SOURCES:.F90=.o)))
$(1)_BUILD_STAMP := $$($(1)_OBJ_DIR)/.ads-benchmark-adapter-build
$(1)_EXECUTABLE := $$(BENCHMARK_EXEC_DIR)/$(1)_manufactured
$(1)_COMPILE_FLAGS = $$(FFLAGS) \
	$$(call benchmark_module_flags,$$($(1)_OBJ_DIR)) \
	-I$$($(1)_OBJ_DIR) -I$$(CORE_MODULE_DIR)

.PHONY: build-$(1) validate-$(1)-sources

validate-$(1)-sources:
	@set -eu; \
	for source in $$($(1)_SOURCES); do \
		if ! test -f "$$$$source"; then \
			printf 'Missing benchmark source: %s\n' "$$$$source" >&2; \
			exit 2; \
		fi; \
	done

build-$(1): validate-build-profile validate-$(1)-sources $$($(1)_EXECUTABLE)

$$($(1)_OBJ_DIR):
	mkdir -p -- "$$@"

$$($(1)_BUILD_STAMP): FORCE | $$($(1)_OBJ_DIR)
	@{ \
		printf '%s\n' 'CONFIG=$$(CONFIG_PATH)'; \
		printf '%s\n' 'BUILD=$$(BUILD)'; \
		printf '%s\n' 'SOURCES=$$($(1)_SOURCES)'; \
		printf '%s\n' 'FF=$$(FF)'; \
		printf '%s\n' 'FFLAGS=$$($(1)_COMPILE_FLAGS)'; \
		printf '%s\n' 'USER_LIB=$$(USER_LIB)'; \
	} > "$$@.tmp"; \
	if test -r "$$@" && cmp -s "$$@.tmp" "$$@"; then \
		$$(RM) -- "$$@.tmp"; \
	else \
		mv -f -- "$$@.tmp" "$$@"; \
	fi

$$($(1)_OBJ_DIR)/%.o: $$(FORTRAN_DIR)/%.F90 $$(MAKEFILE_LIST) | $$($(1)_OBJ_DIR)
	$$(FF) $$($(1)_COMPILE_FLAGS) -o "$$@" -c "$$<"

# The archive is a normal prerequisite, not only a link input: when a public
# core module changes, the refreshed archive makes every adapter consumer
# recompile against the matching module interface before it is relinked.
$$($(1)_OBJECTS): $$($(1)_BUILD_STAMP) $$(CORE_LIBRARY)

$$($(1)_EXECUTABLE): $$($(1)_OBJECTS) $$(CORE_LIBRARY) | $$(BENCHMARK_EXEC_DIR)
	$$(FF) $$($(1)_COMPILE_FLAGS) -o "$$@" $$($(1)_OBJECTS) \
		"$$(CORE_LIBRARY)" $$(USER_LIB)
	@printf 'Built %s\n' "$$@"
endef

$(foreach adapter,$(BENCHMARK_ADAPTERS),$(eval $(call DEFINE_BENCHMARK_ADAPTER,$(adapter))))

$(igrm_heat_OBJ_DIR)/input_data.o: $(PROBLEMS_DIR)/igrm_heat/input_data.F90 \
		$(MAKEFILE_LIST) | $(igrm_heat_OBJ_DIR)
	$(FF) $(igrm_heat_COMPILE_FLAGS) -o "$@" -c "$<"

$(igrm_heat_OBJ_DIR)/RHS_fun.o: $(PROBLEMS_DIR)/igrm_heat/RHS_fun.F90 \
		$(MAKEFILE_LIST) | $(igrm_heat_OBJ_DIR)
	$(FF) $(igrm_heat_COMPILE_FLAGS) -o "$@" -c "$<"

define add-ordered-dependencies
$(if $(word 2,$1),$(eval $(word 2,$1): $(word 1,$1))$(call add-ordered-dependencies,$(wordlist 2,$(words $1),$1)))
endef
$(foreach adapter,$(BENCHMARK_ADAPTERS),$(call add-ordered-dependencies,$($(adapter)_OBJECTS)))

BENCHMARK_EXECUTABLES := $(foreach adapter,$(BENCHMARK_ADAPTERS),$($(adapter)_EXECUTABLE))

.PHONY: build build-adapters validate-benchmark-sources smoke convergence \
	h-convergence p-convergence validate strong weak run-shard resume \
	validate-smoke-run-id show-build-config

validate-benchmark-sources: $(addprefix validate-,$(addsuffix -sources,$(BENCHMARK_ADAPTERS)))

build build-adapters: validate-build-profile validate-benchmark-sources \
	$(BENCHMARK_EXECUTABLES)

validate-smoke-run-id:
	@set -eu; \
	run_id='$(BENCHMARK_RUN_ID)'; \
	case "$$run_id" in \
		''|[!A-Za-z0-9]*|*[!A-Za-z0-9._-]*) \
			printf '%s\n' \
				'BENCHMARK_RUN_ID must use 1-125 ASCII letters, digits, dots, underscores, or hyphens and start with a letter or digit' \
				>&2; \
			exit 2 ;; \
	esac; \
	if test "$${#run_id}" -gt 125; then \
		printf '%s\n' 'BENCHMARK_RUN_ID must contain at most 125 characters' >&2; \
		exit 2; \
	fi

smoke: validate-smoke-run-id build
	@MPIEXEC='$(MPIEXEC)' MPI_NP_FLAG='$(MPI_NP_FLAG)' \
		PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m ads_benchmark.cli run \
		--repository-root '$(ROOT_DIR)' --config-dir '$(CONFIG_DIR)' \
		--profile smoke --run-id '$(BENCHMARK_RUN_ID)-n4' $(BENCHMARK_PLAN_ARGS)
	@MPIEXEC='$(MPIEXEC)' MPI_NP_FLAG='$(MPI_NP_FLAG)' \
		PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m ads_benchmark.cli run \
		--repository-root '$(ROOT_DIR)' --config-dir '$(CONFIG_DIR)' \
		--profile smoke-refined --run-id '$(BENCHMARK_RUN_ID)-n8' \
		$(BENCHMARK_PLAN_ARGS)
	@PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m ads_benchmark.smoke_refinement \
		--repository-root '$(ROOT_DIR)' \
		--coarse-run-id '$(BENCHMARK_RUN_ID)-n4' \
		--refined-run-id '$(BENCHMARK_RUN_ID)-n8'

convergence: validate-run-id build
	@MPIEXEC='$(MPIEXEC)' MPI_NP_FLAG='$(MPI_NP_FLAG)' \
		PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m ads_benchmark.cli run \
		--repository-root '$(ROOT_DIR)' --config-dir '$(CONFIG_DIR)' \
		--profile '$(BENCHMARK_CONVERGENCE_PROFILE)' --run-id '$(RUN_ID)' \
		$(BENCHMARK_PLAN_ARGS)

h-convergence: validate-run-id build
	@MPIEXEC='$(MPIEXEC)' MPI_NP_FLAG='$(MPI_NP_FLAG)' \
		PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m ads_benchmark.cli run \
		--repository-root '$(ROOT_DIR)' --config-dir '$(CONFIG_DIR)' \
		--profile '$(BENCHMARK_H_CONVERGENCE_PROFILE)' --run-id '$(RUN_ID)' \
		$(BENCHMARK_PLAN_ARGS)

p-convergence: validate-run-id build
	@MPIEXEC='$(MPIEXEC)' MPI_NP_FLAG='$(MPI_NP_FLAG)' \
		PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m ads_benchmark.cli run \
		--repository-root '$(ROOT_DIR)' --config-dir '$(CONFIG_DIR)' \
		--profile '$(BENCHMARK_P_CONVERGENCE_PROFILE)' --run-id '$(RUN_ID)' \
		$(BENCHMARK_PLAN_ARGS)

validate: validate-run-id build
	@MPIEXEC='$(MPIEXEC)' MPI_NP_FLAG='$(MPI_NP_FLAG)' \
		PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m ads_benchmark.cli run \
		--repository-root '$(ROOT_DIR)' --config-dir '$(CONFIG_DIR)' \
		--profile '$(BENCHMARK_VALIDATION_PROFILE)' --run-id '$(RUN_ID)' \
		$(BENCHMARK_PLAN_ARGS)
	@PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m ads_benchmark.cli analyze \
		--repository-root '$(ROOT_DIR)' --run-id '$(RUN_ID)' \
		$(BENCHMARK_ANALYZE_ARGS)

strong: validate-run-id build
	@MPIEXEC='$(MPIEXEC)' MPI_NP_FLAG='$(MPI_NP_FLAG)' \
		MPIEXEC_FLAGS='$(MPIEXEC_FLAGS)' \
		PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m ads_benchmark.cli run \
		--repository-root '$(ROOT_DIR)' --config-dir '$(CONFIG_DIR)' \
		--profile '$(BENCHMARK_STRONG_PROFILE)' --run-id '$(RUN_ID)' \
		$(BENCHMARK_PLAN_ARGS)
	@PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m ads_benchmark.cli analyze \
		--repository-root '$(ROOT_DIR)' --run-id '$(RUN_ID)' --plot \
		$(BENCHMARK_ANALYZE_ARGS)

weak: validate-run-id build
	@MPIEXEC='$(MPIEXEC)' MPI_NP_FLAG='$(MPI_NP_FLAG)' \
		MPIEXEC_FLAGS='$(MPIEXEC_FLAGS)' \
		PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m ads_benchmark.cli run \
		--repository-root '$(ROOT_DIR)' --config-dir '$(CONFIG_DIR)' \
		--profile '$(BENCHMARK_WEAK_PROFILE)' --run-id '$(RUN_ID)' \
		$(BENCHMARK_PLAN_ARGS)
	@PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m ads_benchmark.cli analyze \
		--repository-root '$(ROOT_DIR)' --run-id '$(RUN_ID)' --plot \
		$(BENCHMARK_ANALYZE_ARGS)

run-shard: validate-run-id validate-parent-run-id validate-shard-index build
	@MPIEXEC='$(MPIEXEC)' MPI_NP_FLAG='$(MPI_NP_FLAG)' \
		MPIEXEC_FLAGS='$(MPIEXEC_FLAGS)' \
		PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m ads_benchmark.cli run-shard \
		--repository-root '$(ROOT_DIR)' \
		--parent-run-id '$(PARENT_RUN_ID)' \
		--shard-index '$(SHARD_INDEX)' --run-id '$(RUN_ID)' \
		$(BENCHMARK_PLAN_ARGS)

resume: validate-run-id build
	@MPIEXEC='$(MPIEXEC)' MPI_NP_FLAG='$(MPI_NP_FLAG)' \
		MPIEXEC_FLAGS='$(MPIEXEC_FLAGS)' \
		PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m ads_benchmark.cli resume \
		--repository-root '$(ROOT_DIR)' --config-dir '$(CONFIG_DIR)' \
		--profile '$(BENCHMARK_RESUME_PROFILE)' --run-id '$(RUN_ID)' \
		$(BENCHMARK_PLAN_ARGS)

show-build-config:
	@printf '%-24s %s\n' \
		'CONFIG' '$(CONFIG_PATH)' \
		'BUILD' '$(BUILD)' \
		'BENCHMARK_BUILD_ROOT' '$(BENCHMARK_BUILD_ROOT_PATH)' \
		'PROFILE_BUILD_ROOT' '$(BENCHMARK_PROFILE_BUILD_ROOT)' \
		'CORE_LIBRARY' '$(CORE_LIBRARY)' \
		'OBJECT_ROOT' '$(BENCHMARK_OBJ_ROOT)' \
		'EXECUTABLES' '$(BENCHMARK_EXECUTABLES)' \
		'COMPILER' '$(COMPILER)' \
		'MPIEXEC' '$(MPIEXEC)'
endif

PYTHON ?= python3

.PHONY: all plan shard-plan merge-shards analyze self-test clean clean-build \
	prepare-build-root validate-run-id validate-parent-run-id \
	validate-shard-index validate-shard-runs validate-build-root \
	validate-build-ownership help

all: self-test

plan:
	@PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m ads_benchmark.cli plan \
		--repository-root '$(ROOT_DIR)' --config-dir '$(CONFIG_DIR)' \
		--profile '$(BENCHMARK_PROFILE)' $(BENCHMARK_PLAN_ARGS)

shard-plan: validate-run-id
	@PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m ads_benchmark.cli shard \
		--repository-root '$(ROOT_DIR)' --config-dir '$(CONFIG_DIR)' \
		--profile '$(BENCHMARK_SHARD_PROFILE)' --run-id '$(RUN_ID)' \
		--strategy '$(BENCHMARK_SHARD_STRATEGY)' \
		$(if $(strip $(BENCHMARK_SHARD_COUNT)),--shard-count '$(BENCHMARK_SHARD_COUNT)') \
		$(BENCHMARK_PLAN_ARGS)

validate-run-id:
	@set -eu; \
	run_id='$(RUN_ID)'; \
	case "$$run_id" in \
		''|[!A-Za-z0-9]*|*[!A-Za-z0-9._-]*) \
			printf '%s\n' \
				'RUN_ID must use 1-128 ASCII letters, digits, dots, underscores, or hyphens and start with a letter or digit' \
				>&2; \
			exit 2 ;; \
	esac; \
	if test "$${#run_id}" -gt 128; then \
		printf '%s\n' 'RUN_ID must contain at most 128 characters' >&2; \
		exit 2; \
	fi

validate-parent-run-id:
	@set -eu; \
	parent_run_id='$(PARENT_RUN_ID)'; \
	case "$$parent_run_id" in \
		''|[!A-Za-z0-9]*|*[!A-Za-z0-9._-]*) \
			printf '%s\n' \
				'PARENT_RUN_ID must use 1-128 ASCII letters, digits, dots, underscores, or hyphens and start with a letter or digit' \
				>&2; \
			exit 2 ;; \
	esac; \
	if test "$${#parent_run_id}" -gt 128; then \
		printf '%s\n' 'PARENT_RUN_ID must contain at most 128 characters' >&2; \
		exit 2; \
	fi

validate-shard-index:
	@case '$(SHARD_INDEX)' in \
		''|*[!0-9]*) \
			printf '%s\n' 'SHARD_INDEX must be a nonnegative integer' >&2; \
			exit 2 ;; \
	esac

validate-shard-runs:
	@test -n '$(strip $(SHARD_RUNS))' || { \
		printf '%s\n' 'SHARD_RUNS must list every completed child run ID' >&2; \
		exit 2; \
	}

merge-shards: validate-parent-run-id validate-shard-runs
	@PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m ads_benchmark.cli merge-shards \
		--repository-root '$(ROOT_DIR)' \
		--parent-run-id '$(PARENT_RUN_ID)' \
		$(foreach run,$(SHARD_RUNS),--shard-run '$(run)') \
		--plot $(filter-out --plot,$(BENCHMARK_ANALYZE_ARGS))

analyze: validate-run-id
	@PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m ads_benchmark.cli analyze \
		--repository-root '$(ROOT_DIR)' --run-id '$(RUN_ID)' \
		$(BENCHMARK_ANALYZE_ARGS)

self-test:
	@PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m unittest discover \
		-s selftests -p 'test_*.py' -v

validate-build-root:
	@set -eu; \
	actual=$$(realpath -m -- '$(BENCHMARK_BUILD_ROOT_PATH)'); \
	expected=$$(realpath -m -- '$(BENCHMARKING_DIR)/build'); \
	if test "$$actual" != "$$expected"; then \
		printf 'Unsafe BENCHMARK_BUILD_ROOT: %s\n' \
			'$(BENCHMARK_BUILD_ROOT_PATH)' >&2; exit 2; \
	fi; \
	if test -L '$(BENCHMARK_BUILD_ROOT_PATH)'; then \
		printf 'Unsafe benchmark build symlink: %s\n' \
			'$(BENCHMARK_BUILD_ROOT_PATH)' >&2; exit 2; \
	fi

validate-build-ownership: validate-build-root
	@set -eu; \
	root='$(BENCHMARK_BUILD_ROOT_PATH)'; \
	marker='$(BENCHMARK_BUILD_MARKER)'; \
	if ! test -e "$$root"; then exit 0; fi; \
	if ! test -d "$$root"; then \
		printf 'Benchmark build root is not a directory: %s\n' "$$root" >&2; \
		exit 2; \
	fi; \
	if test -z "$$(find "$$root" -mindepth 1 -maxdepth 1 -print -quit)"; then \
		exit 0; \
	fi; \
	expected="ADS_MPI_BENCHMARK_BUILD=$$(realpath -m -- "$$root")"; \
	if test -L "$$marker" || ! test -f "$$marker" || \
		test "$$(cat -- "$$marker" 2>/dev/null || :)" != "$$expected"; then \
		printf 'Refusing unowned benchmark build root: %s\n' "$$root" >&2; \
		exit 2; \
	fi

prepare-build-root: validate-build-ownership
	@set -eu; \
	root='$(BENCHMARK_BUILD_ROOT_PATH)'; \
	marker='$(BENCHMARK_BUILD_MARKER)'; \
	expected="ADS_MPI_BENCHMARK_BUILD=$$(realpath -m -- "$$root")"; \
	mkdir -p -- "$$root"; \
	if ! test -f "$$marker"; then \
		marker_tmp="$$marker.tmp"; \
		printf '%s\n' "$$expected" > "$$marker_tmp"; \
		mv -f -- "$$marker_tmp" "$$marker"; \
	fi

clean clean-build: validate-build-ownership
	@find '$(BENCHMARKING_DIR)' -type f \( -name '*.pyc' -o -name '*.pyo' \) \
		-delete
	@find '$(BENCHMARKING_DIR)' -depth -type d -name __pycache__ -empty -delete
	@set -eu; \
	root='$(BENCHMARK_BUILD_ROOT_PATH)'; \
	marker='$(BENCHMARK_BUILD_MARKER)'; \
	if ! test -e "$$root"; then exit 0; fi; \
	if ! test -d "$$root"; then \
		printf 'Benchmark build root is not a directory: %s\n' "$$root" >&2; \
		exit 2; \
	fi; \
	expected="ADS_MPI_BENCHMARK_BUILD=$$(realpath -m -- "$$root")"; \
	if test -z "$$(find "$$root" -mindepth 1 -maxdepth 1 -print -quit)"; then \
		rmdir -- "$$root"; exit 0; \
	fi; \
	if test -L "$$marker" || ! test -f "$$marker" || \
		test "$$(cat -- "$$marker" 2>/dev/null || :)" != "$$expected"; then \
		printf 'Refusing unowned benchmark build root: %s\n' "$$root" >&2; \
		exit 2; \
	fi; \
	$(RM) -r -- "$$root"

help:
	@printf '%s\n' \
		'ADS MPI benchmark framework' \
		'' \
		'  make plan [BENCHMARK_PROFILE=smoke]' \
		'      validate and summarize a plan without writing results' \
		'  make plan BENCHMARK_PLAN_ARGS="--write-manifest --run-id NAME"' \
		'      create benchmarks/NAME/manifest.json exclusively' \
		'  make self-test' \
		'      run dependency-free planner/engine contract tests' \
		'  make build [BUILD=debug|release]' \
		'      build all three manufactured benchmark adapters' \
		'  make build-igrm_l2 | build-igrm_heat | build-pure_diffusion_igrm' \
		'      build one manufactured adapter in the selected profile' \
		'  make smoke BENCHMARK_RUN_ID=NAME' \
		'      run the real N=4/N=8 smoke pair and verify all nine refinements' \
		'  make convergence RUN_ID=NAME [BENCHMARK_CONVERGENCE_PROFILE=temporal-full]' \
		'      build release adapters and execute a new frozen temporal manifest' \
		'  make h-convergence RUN_ID=NAME [BENCHMARK_H_CONVERGENCE_PROFILE=h-convergence-full]' \
		'      execute a new frozen h-convergence manifest (use h-convergence-smoke for smoke)' \
		'  make p-convergence RUN_ID=NAME [BENCHMARK_P_CONVERGENCE_PROFILE=p-convergence-full]' \
		'      execute isotropic p convergence; p-anisotropic-{smoke,full} are also selectable' \
		'  make validate RUN_ID=NAME [BENCHMARK_VALIDATION_PROFILE=validation-full]' \
		'      run full-field validation; BENCHMARK_PLAN_ARGS must declare MPI/CPU slots' \
		'  make strong RUN_ID=NAME [BENCHMARK_STRONG_PROFILE=cluster-scaling]' \
		'      run release strong scaling with field gates; declare MPI/CPU slots' \
		'  make weak RUN_ID=NAME [BENCHMARK_WEAK_PROFILE=cluster-weak-scaling]' \
		'      run release weak scaling with field gates; declare MPI/CPU slots' \
		'  make shard-plan RUN_ID=PARENT [BENCHMARK_SHARD_PROFILE=cluster-weak-scaling]' \
		'      freeze a plan and generate deterministic group or index shards' \
		'  make run-shard PARENT_RUN_ID=PARENT SHARD_INDEX=N RUN_ID=CHILD' \
		'      build release adapters and execute one shard; pass capacity flags in BENCHMARK_PLAN_ARGS' \
		'  make merge-shards PARENT_RUN_ID=PARENT SHARD_RUNS="CHILD ..."' \
		'      verify the complete child-run set and analyze the parent plan' \
		'  make resume RUN_ID=NAME [BENCHMARK_RESUME_PROFILE=temporal-full]' \
		'      verify the manifest and retry only incomplete/failed cases' \
		'  make analyze RUN_ID=NAME [BENCHMARK_ANALYZE_ARGS=--plot]' \
		'      write the registered convergence JSON/CSV report and optional plot' \
		'  make show-build-config [BUILD=debug|release]' \
		'      print the selected isolated build paths and toolchain' \
		'  make clean-build' \
		'      remove only benchmark framework caches and marker-owned builds' \
		'' \
		'Profiles: smoke smoke-refined temporal-validation temporal-full h-convergence-smoke h-convergence-full' \
		'          p-convergence-smoke p-convergence-full p-anisotropic-smoke p-anisotropic-full' \
		'          validation-smoke validation-full strong-scaling-smoke local-scaling cluster-scaling' \
		'          weak-scaling-smoke local-weak-scaling cluster-weak-scaling'
