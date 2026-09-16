.DEFAULT_GOAL := help
.NOTPARALLEL:

BENCHMARKING_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
ROOT_DIR := $(abspath $(BENCHMARKING_DIR)/..)
CONFIG_DIR ?= $(BENCHMARKING_DIR)/configs
PYTHON ?= python3

BENCHMARK_PROFILE ?= smoke
BENCHMARK_PLAN_ARGS ?=
BENCHMARK_BUILD_ROOT ?= $(BENCHMARKING_DIR)/build
BENCHMARK_BUILD_ROOT_PATH := $(abspath $(BENCHMARK_BUILD_ROOT))
BENCHMARK_BUILD_MARKER := $(BENCHMARK_BUILD_ROOT_PATH)/.ads-benchmark-build

.PHONY: all plan self-test clean clean-build validate-build-root \
	validate-build-ownership help

all: self-test

plan:
	@PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m ads_benchmark.cli plan \
		--repository-root '$(ROOT_DIR)' --config-dir '$(CONFIG_DIR)' \
		--profile '$(BENCHMARK_PROFILE)' $(BENCHMARK_PLAN_ARGS)

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
		'ADS MPI benchmark planning framework' \
		'' \
		'  make plan [BENCHMARK_PROFILE=smoke]' \
		'      validate and summarize a plan without writing results' \
		'  make plan BENCHMARK_PLAN_ARGS="--write-manifest --run-id NAME"' \
		'      create benchmarks/NAME/manifest.json exclusively' \
		'  make self-test' \
		'      run dependency-free planner/engine contract tests' \
		'  make clean-build' \
		'      remove only benchmark framework caches and marker-owned builds' \
		'' \
		'Profiles: smoke temporal-full local-scaling cluster-scaling'
