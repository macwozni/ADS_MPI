TEST_GROUP_NAME := driver
SUITES := driver_cli

include ../make/group.mk

.PHONY: run-cli run-smoke run-integration

run-cli:
	$(call rebuild_and_run_suites,$(SUITES),run-cli)

run-smoke:
	$(call rebuild_and_run_suites,$(SUITES),run-smoke)

# Keep the historical integration behavior: build what is stale and run the
# integration checks without first forcing a clean rebuild.
run-integration:
	$(call run_suites,$(SUITES),run-integration)
