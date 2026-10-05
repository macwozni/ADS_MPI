"""Repository paths shared by the benchmark framework tests."""

from pathlib import Path


TEST_DIRECTORY = Path(__file__).resolve().parent
REPOSITORY_ROOT = TEST_DIRECTORY.parent.parent
BENCHMARKING_ROOT = REPOSITORY_ROOT / "benchmarking"
CONFIG_DIRECTORY = BENCHMARKING_ROOT / "configs"
