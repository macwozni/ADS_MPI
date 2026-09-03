.DEFAULT_GOAL := all

PYTHON ?= python3
TEST_SCRIPT := test_make_hierarchy.py

.PHONY: all tests run clean

all:
	@$(PYTHON) -m py_compile $(TEST_SCRIPT)

tests: all

run: all
	@$(PYTHON) ./$(TEST_SCRIPT)

clean:
	$(RM) -r -- __pycache__
