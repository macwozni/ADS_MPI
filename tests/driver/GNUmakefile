TEST_GROUP_NAME := driver
SUITES := driver_cli
GROUP_RUN_TIMEOUT = $(PERFORMANCE_SUITE_TIMEOUT)

include ../make/group.mk

.PHONY: run-cli run-smoke run-integration run-performance \
	run-performance-self-test clean-performance

run-cli:
	$(call rebuild_and_run_suites,$(SUITES),run-cli)

run-smoke:
	$(call rebuild_and_run_suites,$(SUITES),run-smoke)

# Keep the historical integration behavior: build what is stale and run the
# integration checks without first forcing a clean rebuild.
run-integration:
	$(call run_suites,$(SUITES),run-integration)

run-performance:
	$(call run_suites_with_timeout,$(SUITES),run-performance,$(PERFORMANCE_SUITE_TIMEOUT))

run-performance-self-test:
	$(call run_suites,$(SUITES),run-performance-self-test)

clean-performance:
	$(call run_suites,$(SUITES),clean-performance)
