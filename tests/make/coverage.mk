# Loaded through GNU Make's MAKEFILES mechanism only by `make test-coverage`.
# Keeping instrumentation out of the FFLAGS environment is important: several
# error-contract tests launch a nested make, and environment-derived FFLAGS
# would otherwise receive every suite-local flag a second time.
ifndef ADS_TEST_COVERAGE_FLAGS_INCLUDED
ADS_TEST_COVERAGE_FLAGS_INCLUDED := 1
FFLAGS += $(COVERAGE_FLAGS)
endif

# Production problem makefiles normally keep private objects beside their
# sources. During coverage, route them into the owned build tree as well, so a
# failed or successful instrumented run cannot contaminate the next normal
# problem build. The value is recursive because BUILD_ROOT is defined by the
# including makefile after GNU Make has loaded this injected file.
PROBLEM_OBJ_DIR ?= $(abspath $(BUILD_ROOT))/problem-objects/$(notdir $(CURDIR))
