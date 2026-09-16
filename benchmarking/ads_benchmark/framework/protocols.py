"""Structural extension contracts used by the generic engine."""

from __future__ import annotations

from typing import Mapping, Protocol, Sequence

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


class Launcher(Protocol):
    """Wrap a payload command and define runtime environment additions."""

    name: str

    def validate_case(self, case: CaseSpec) -> None:
        """Raise when the launcher cannot realize the recorded topology."""

    def command(self, payload: Sequence[str], case: CaseSpec) -> Sequence[str]:
        """Return the complete argv executed by the engine."""

    def environment(self, case: CaseSpec) -> Mapping[str, str]:
        """Return launcher-specific environment additions."""
