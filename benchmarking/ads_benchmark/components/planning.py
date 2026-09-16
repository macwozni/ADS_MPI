"""Stage-one planning adapters and reusable local/MPI launchers."""

from __future__ import annotations

from dataclasses import dataclass
import os
from typing import Mapping, Sequence

from ..framework.errors import ExecutionError
from ..framework.model import CaseSpec, ExecutionContext


@dataclass(frozen=True)
class PlanningAdapter:
    """Declare a real ADS problem while execution harnesses are still absent."""

    name: str
    execution_ready: bool = False

    def validate_case(self, case: CaseSpec) -> None:
        if case.problem != self.name:
            raise ValueError(f"adapter {self.name} received problem {case.problem}")

    def build_payload_command(
        self, case: CaseSpec, context: ExecutionContext
    ) -> Sequence[str]:
        raise ExecutionError(
            f"adapter {self.name} has no solver harness until benchmark stage 2"
        )

    def parse_result(self, stdout: str, stderr: str) -> Mapping[str, object]:
        raise ExecutionError(
            f"adapter {self.name} has no result parser until benchmark stage 2"
        )


@dataclass(frozen=True)
class DirectLauncher:
    name: str = "direct"

    def validate_case(self, case: CaseSpec) -> None:
        if case.mpi.ranks != 1:
            raise ValueError("direct launcher supports exactly one MPI rank")

    def command(self, payload: Sequence[str], case: CaseSpec) -> Sequence[str]:
        return tuple(payload)

    def environment(self, case: CaseSpec) -> Mapping[str, str]:
        return {}


@dataclass(frozen=True)
class MpiLauncher:
    """Data-driven MPI wrapper used once executable adapters are installed."""

    executable: str
    rank_flag: str
    name: str = "mpi"

    def validate_case(self, case: CaseSpec) -> None:
        if not self.executable or not self.rank_flag:
            raise ValueError("MPI launcher executable and rank flag must be nonempty")

    def command(self, payload: Sequence[str], case: CaseSpec) -> Sequence[str]:
        return (
            self.executable,
            self.rank_flag,
            str(case.mpi.ranks),
            *payload,
        )

    def environment(self, case: CaseSpec) -> Mapping[str, str]:
        return {}


def default_mpi_launcher() -> MpiLauncher:
    return MpiLauncher(
        executable=os.environ.get("MPIEXEC", "mpiexec"),
        rank_flag=os.environ.get("MPI_NP_FLAG", "-n"),
    )
