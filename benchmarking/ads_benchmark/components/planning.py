"""Planning-only adapters and reusable local/MPI launchers."""

from __future__ import annotations

from dataclasses import dataclass
import os
import shlex
import string
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

    def validate_resources(self, case: CaseSpec) -> None:
        self.validate_case(case)

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


_TEMPLATE_FIELDS = frozenset(
    {"ranks", "threads", "procx", "procy", "procz", "payload"}
)


def _template_fields(arguments: Sequence[str]) -> tuple[str, ...]:
    formatter = string.Formatter()
    fields: list[str] = []
    for argument in arguments:
        for _, field_name, format_spec, conversion in formatter.parse(argument):
            if field_name is None:
                continue
            if field_name not in _TEMPLATE_FIELDS:
                raise ValueError(
                    f"unknown launcher template placeholder {{{field_name}}}"
                )
            if format_spec or conversion:
                raise ValueError(
                    "launcher template placeholders do not accept formatting"
                )
            fields.append(field_name)
    return tuple(fields)


@dataclass(frozen=True)
class TemplateLauncher:
    """Shell-free command template for schedulers such as SLURM ``srun``."""

    arguments: tuple[str, ...]
    name: str = "mpi"
    available_mpi_slots: int | None = None
    available_cpu_slots: int | None = None
    allocated_threads_per_rank: int | None = None

    def _fields(self) -> tuple[str, ...]:
        if not self.arguments:
            raise ValueError("launcher template must be nonempty")
        if any(not argument or "\0" in argument for argument in self.arguments):
            raise ValueError(
                "launcher template arguments must be nonempty NUL-free strings"
            )
        fields = _template_fields(self.arguments)
        if (
            fields.count("payload") != 1
            or not self.arguments
            or self.arguments[-1] != "{payload}"
        ):
            raise ValueError(
                "launcher template requires exactly one final standalone {payload}"
            )
        if fields.count("ranks") < 1:
            raise ValueError("launcher template must expose {ranks}")
        for capacity, label in (
            (self.available_mpi_slots, "available MPI slots"),
            (self.available_cpu_slots, "available CPU slots"),
            (self.allocated_threads_per_rank, "allocated threads per rank"),
        ):
            if capacity is not None and (
                type(capacity) is not int or capacity <= 0
            ):
                raise ValueError(f"{label} must be a positive integer")
        return fields

    def validate_case(self, case: CaseSpec) -> None:
        fields = self._fields()
        if case.openmp_threads > 1 and "threads" not in fields:
            raise ValueError(
                "hybrid execution requires {threads} in the launcher template"
            )
        if case.openmp_threads > 1 and (
            case.openmp_dynamic is not False
            or case.openmp_proc_bind is None
            or case.openmp_places is None
        ):
            raise ValueError(
                "hybrid execution requires OMP_DYNAMIC=FALSE and explicit binding"
            )

    def validate_resources(self, case: CaseSpec) -> None:
        self.validate_case(case)
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
        if (
            self.allocated_threads_per_rank is not None
            and case.openmp_threads > self.allocated_threads_per_rank
        ):
            raise ValueError(
                f"case requests {case.openmp_threads} OpenMP threads per rank "
                f"but the allocation provides only "
                f"{self.allocated_threads_per_rank}"
            )

    def command(self, payload: Sequence[str], case: CaseSpec) -> Sequence[str]:
        self.validate_case(case)
        procx, procy, procz = case.mpi.process_grid
        values = {
            "ranks": str(case.mpi.ranks),
            "threads": str(case.openmp_threads),
            "procx": str(procx),
            "procy": str(procy),
            "procz": str(procz),
        }
        command: list[str] = []
        for argument in self.arguments:
            if argument == "{payload}":
                command.extend(payload)
            else:
                command.append(argument.format_map(values))
        return tuple(command)

    def environment(self, case: CaseSpec) -> Mapping[str, str]:
        return {}


def _positive_environment_integer(name: str) -> int | None:
    text = os.environ.get(name)
    if text is None:
        return None
    try:
        value = int(text)
    except ValueError as error:
        raise ValueError(f"{name} must be a positive integer") from error
    if value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _constrained_capacity(
    declared: int | None,
    detected: int | None,
) -> int | None:
    """Use the stricter of CLI policy and the scheduler's real allocation."""

    if declared is None:
        return detected
    if detected is None:
        return declared
    return min(declared, detected)


def default_mpi_launcher(
    *,
    available_mpi_slots: int | None = None,
    available_cpu_slots: int | None = None,
    command_template: str | None = None,
) -> MpiLauncher | TemplateLauncher:
    template = (
        command_template
        if command_template is not None
        else os.environ.get("BENCHMARK_LAUNCHER_TEMPLATE")
    )
    if template is not None and template.strip():
        try:
            arguments = tuple(shlex.split(template))
        except ValueError as error:
            raise ValueError(f"invalid launcher template: {error}") from error
        slurm_tasks = _positive_environment_integer("SLURM_NTASKS")
        slurm_threads = _positive_environment_integer("SLURM_CPUS_PER_TASK")
        detected_mpi_slots = _constrained_capacity(
            available_mpi_slots, slurm_tasks
        )
        slurm_cpu_slots = (
            slurm_tasks * slurm_threads
            if slurm_tasks is not None and slurm_threads is not None
            else None
        )
        detected_cpu_slots = _constrained_capacity(
            available_cpu_slots, slurm_cpu_slots
        )
        launcher = TemplateLauncher(
            arguments=arguments,
            available_mpi_slots=detected_mpi_slots,
            available_cpu_slots=detected_cpu_slots,
            allocated_threads_per_rank=slurm_threads,
        )
        # Fail on malformed templates at composition time, before planning.
        launcher._fields()
        return launcher

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
