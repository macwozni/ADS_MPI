.DEFAULT_GOAL := all

PYTHON ?= python3
REPOSITORY_ROOT := $(abspath ../..)
BENCHMARKING_ROOT := $(REPOSITORY_ROOT)/benchmarking
TEST_SCRIPTS := $(wildcard *.py)
PYTHON_ENV := env PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$(BENCHMARKING_ROOT)"

.PHONY: all tests run clean

all:
	@$(PYTHON_ENV) $(PYTHON) -m py_compile $(TEST_SCRIPTS)

tests run: all
	@$(PYTHON_ENV) $(PYTHON) -m unittest discover -s . -p 'test_*.py' -v

clean:
	$(RM) -r -- __pycache__
