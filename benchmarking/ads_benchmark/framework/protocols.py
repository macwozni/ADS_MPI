"""Structural extension contracts used by the generic engine."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, Mapping, Protocol, Sequence

from .model import CaseSpec, ExecutionContext


class ProblemAdapter(Protocol):
    """Problem-specific behavior; process management stays in the engine."""

    name: str
    execution_ready: bool

    def validate_case(self, case: CaseSpec) -> None:
        """Raise on a problem-specific incompatibility."""

    def build_payload_command(
        self, case: CaseSpec, context: ExecutionContext
    ) -> Sequence[str]:
        """Return an argv payload without shell syntax or launcher wrapping."""

    def parse_result(self, stdout: str, stderr: str) -> Mapping[str, object]:
        """Parse the problem-domain record from a successful process."""

    def validate_result(
        self, case: CaseSpec, result: Mapping[str, object]
    ) -> None:
        """Raise when a parsed result does not describe its planned case."""


class Launcher(Protocol):
    """Wrap a payload command and define runtime environment additions."""

    name: str

    def validate_case(self, case: CaseSpec) -> None:
        """Raise when the launcher cannot realize the recorded topology."""

    def command(self, payload: Sequence[str], case: CaseSpec) -> Sequence[str]:
        """Return the complete argv executed by the engine."""

    def validate_resources(self, case: CaseSpec) -> None:
        """Raise when a case exceeds the declared runtime allocation."""

    def environment(self, case: CaseSpec) -> Mapping[str, str]:
        """Return launcher-specific environment additions."""


class LauncherFailureClassifier(Protocol):
    """Optional launcher extension for explicit control-plane failures."""

    def classify_process_failure(
        self,
        *,
        return_code: int,
        stdout: str,
        stderr: str,
        case_directory: Path,
    ) -> Literal["mpi"] | None:
        """Return ``"mpi"`` only for a positively identified launcher failure."""
