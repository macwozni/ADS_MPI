"""Planning-only adapters and reusable local/MPI launchers."""

from __future__ import annotations

from dataclasses import dataclass
import os
import shlex
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
            f"planning-only adapter {self.name} has no solver harness"
        )

    def parse_result(self, stdout: str, stderr: str) -> Mapping[str, object]:
        raise ExecutionError(
            f"planning-only adapter {self.name} has no result parser"
        )

    def validate_result(
        self, case: CaseSpec, result: Mapping[str, object]
    ) -> None:
        raise ExecutionError(
            f"planning-only adapter {self.name} has no result validator"
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
    available_mpi_slots: int | None = None
    available_cpu_slots: int | None = None
    extra_args: tuple[str, ...] = ()

    def validate_case(self, case: CaseSpec) -> None:
        if not self.executable or not self.rank_flag:
            raise ValueError("MPI launcher executable and rank flag must be nonempty")
        if any(
            not isinstance(argument, str) or not argument or "\0" in argument
            for argument in self.extra_args
        ):
            raise ValueError("MPI launcher extra arguments must be nonempty strings")
        if (
            self.available_mpi_slots is not None
            and (
                type(self.available_mpi_slots) is not int
                or self.available_mpi_slots <= 0
            )
        ):
            raise ValueError("available MPI slots must be a positive integer")
        if (
            self.available_cpu_slots is not None
            and (
                type(self.available_cpu_slots) is not int
                or self.available_cpu_slots <= 0
            )
        ):
            raise ValueError("available CPU slots must be a positive integer")

    def validate_resources(self, case: CaseSpec) -> None:
        """Reject a selected case that exceeds the declared allocation."""

        if (
            self.available_mpi_slots is not None
            and case.mpi.ranks > self.available_mpi_slots
        ):
            raise ValueError(
                f"case requests {case.mpi.ranks} MPI ranks but only "
                f"{self.available_mpi_slots} MPI slots are available"
            )
        requested_cpu_slots = case.mpi.ranks * case.openmp_threads
        if (
            self.available_cpu_slots is not None
            and requested_cpu_slots > self.available_cpu_slots
        ):
            raise ValueError(
                f"case requests {case.mpi.ranks} MPI ranks * "
                f"{case.openmp_threads} OpenMP threads = "
                f"{requested_cpu_slots} CPU slots but only "
                f"{self.available_cpu_slots} CPU slots are available"
            )

    def command(self, payload: Sequence[str], case: CaseSpec) -> Sequence[str]:
        return (
            self.executable,
            *self.extra_args,
            self.rank_flag,
            str(case.mpi.ranks),
            *payload,
        )

    def environment(self, case: CaseSpec) -> Mapping[str, str]:
        return {}


def default_mpi_launcher(
    *,
    available_mpi_slots: int | None = None,
    available_cpu_slots: int | None = None,
) -> MpiLauncher:
    try:
        executable_parts = tuple(shlex.split(os.environ.get("MPIEXEC", "mpiexec")))
        flag_parts = tuple(shlex.split(os.environ.get("MPIEXEC_FLAGS", "")))
    except ValueError as error:
        raise ValueError(f"invalid MPI launcher configuration: {error}") from error
    if not executable_parts:
        raise ValueError("MPIEXEC must name a launcher executable")
    return MpiLauncher(
        executable=executable_parts[0],
        rank_flag=os.environ.get("MPI_NP_FLAG", "-n"),
        available_mpi_slots=available_mpi_slots,
        available_cpu_slots=available_cpu_slots,
        extra_args=(*executable_parts[1:], *flag_parts),
    )
