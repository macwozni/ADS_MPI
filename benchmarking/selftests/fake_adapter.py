"""Executable extension used only by the framework contract self-tests."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import sys
from typing import Mapping, Sequence

from ads_benchmark.framework.model import CaseSpec, ExecutionContext


PROGRAM = Path(__file__).with_name("fake_program.py")
RESULT_PREFIX = "ADS_BENCHMARK_RESULT "


@dataclass(frozen=True)
class FakeAdapter:
    mode: str = "success"
    name: str = "fake"
    execution_ready: bool = True

    def validate_case(self, case: CaseSpec) -> None:
        if self.mode == "reject":
            raise ValueError("intentional fake validation rejection")
        if case.problem != self.name:
            raise ValueError("fake adapter received the wrong problem")

    def build_payload_command(
        self, case: CaseSpec, context: ExecutionContext
    ) -> Sequence[str]:
        program_mode = "success" if self.mode in {"nan", "nul-argv"} else self.mode
        if self.mode == "nul-argv":
            return (sys.executable, "bad\0argument")
        return (
            sys.executable,
            str(PROGRAM),
            "--mode",
            program_mode,
            "--steps",
            str(case.time.steps),
            "--expected-cwd",
            str(context.case_directory),
        )

    def parse_result(self, stdout: str, stderr: str) -> Mapping[str, object]:
        if self.mode == "nan":
            return {"invalid": float("nan")}
        matches = [line for line in stdout.splitlines() if line.startswith(RESULT_PREFIX)]
        if len(matches) != 1:
            raise ValueError("expected exactly one tagged fake result")
        document = json.loads(matches[0][len(RESULT_PREFIX) :])
        if not isinstance(document, dict):
            raise ValueError("fake result must be an object")
        return document
