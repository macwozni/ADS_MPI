.DEFAULT_GOAL := all

PYTHON ?= python3
SCRIPTS := check_coverage.py test_check_coverage.py

.PHONY: all tests run clean

all:
	@$(PYTHON) -m py_compile $(SCRIPTS)

tests: all

run: all
	@$(PYTHON) ./test_check_coverage.py

clean:
	$(RM) -r -- __pycache__
