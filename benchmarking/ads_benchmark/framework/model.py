"""Immutable public data model shared by every benchmark family."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping


Vector3 = tuple[int, int, int]

# The executable accepts up to 257 points when no field artifact is retained.
# A written CSV row has a fixed 150-byte representation; 76^3 rows plus the
# header fit below ResultStore's 64 MiB generated-artifact ceiling, while 77^3
# rows do not.
MAX_SAMPLE_POINTS_PER_AXIS = 257
MAX_WRITTEN_SAMPLE_POINTS_PER_AXIS = 76


@dataclass(frozen=True)
class TimeSpec:
    """Exact rational time data serialized as decimal or fraction strings."""

    final_time: str
    time_step: str
    steps: int

    def to_dict(self) -> dict[str, object]:
        return {
            "final_time": self.final_time,
            "time_step": self.time_step,
            "steps": self.steps,
        }


@dataclass(frozen=True)
class MpiSpec:
    ranks: int
    process_grid: Vector3

    def to_dict(self) -> dict[str, object]:
        return {"ranks": self.ranks, "process_grid": list(self.process_grid)}


@dataclass(frozen=True)
class MeasurementSpec:
    warmups: int
    samples: int
    timeout_seconds: str
    minimum_sample_seconds: str | None = None

    def to_dict(self) -> dict[str, object]:
        document: dict[str, object] = {
            "warmups": self.warmups,
            "samples": self.samples,
            "timeout_seconds": self.timeout_seconds,
        }
        if self.minimum_sample_seconds is not None:
            document["minimum_sample_seconds"] = self.minimum_sample_seconds
        return document


@dataclass(frozen=True)
class SamplingSpec:
    """Regular-grid field sampling requested from the solver harness."""

    points_per_axis: int
    write_samples: bool

    def to_dict(self) -> dict[str, object]:
        return {
            "points_per_axis": self.points_per_axis,
            "write_samples": self.write_samples,
        }


@dataclass(frozen=True)
class WeakScalingSpec:
    """Weak-scaling semantics attached to an expanded case.

    ``local_elements`` is the constant spatial workload declared by the
    profile.  ``role`` keeps serial correctness helpers out of performance
    series without giving them a separate execution path.
    """

    local_elements: Vector3
    workload_basis: str
    role: str

    def to_dict(self) -> dict[str, object]:
        return {
            "local_elements": list(self.local_elements),
            "workload_basis": self.workload_basis,
            "role": self.role,
        }


@dataclass(frozen=True)
class CaseSpec:
    """One fully expanded, normalized benchmark configuration."""

    family: str
    problem: str
    scheme: str
    exact_case: str
    time: TimeSpec
    mesh: Vector3
    test_degree: Vector3
    trial_degree: Vector3
    mpi: MpiSpec
    openmp_threads: int
    sampling: SamplingSpec
    measurement: MeasurementSpec
    build_profile: str
    launcher: str
    openmp_dynamic: bool | None = None
    openmp_proc_bind: str | None = None
    openmp_places: str | None = None
    weak_scaling: WeakScalingSpec | None = None

    def to_dict(self) -> dict[str, object]:
        openmp: dict[str, object] = {"threads": self.openmp_threads}
        if self.openmp_dynamic is not None:
            openmp["dynamic"] = self.openmp_dynamic
        if self.openmp_proc_bind is not None:
            openmp["proc_bind"] = self.openmp_proc_bind
        if self.openmp_places is not None:
            openmp["places"] = self.openmp_places
        document: dict[str, object] = {
            "family": self.family,
            "problem": self.problem,
            "scheme": self.scheme,
            "exact_case": self.exact_case,
            "time": self.time.to_dict(),
            "mesh": {"elements": list(self.mesh)},
            "spaces": {
                "test_degree": list(self.test_degree),
                "trial_degree": list(self.trial_degree),
            },
            "mpi": self.mpi.to_dict(),
            "openmp": openmp,
            "sampling": self.sampling.to_dict(),
            "measurement": self.measurement.to_dict(),
            "build": {"profile": self.build_profile},
            "launcher": self.launcher,
        }
        if self.weak_scaling is not None:
            document["weak_scaling"] = self.weak_scaling.to_dict()
        return document


@dataclass(frozen=True)
class PlannedCase:
    case_id: str
    spec: CaseSpec

    def to_dict(self) -> dict[str, object]:
        return {"case_id": self.case_id, "configuration": self.spec.to_dict()}


@dataclass(frozen=True)
class RepositoryState:
    commit: str
    dirty: bool
    worktree_fingerprint: str | None = None

    def to_dict(self) -> dict[str, object]:
        document: dict[str, object] = {
            "commit": self.commit,
            "dirty": self.dirty,
        }
        if self.worktree_fingerprint is not None:
            document["worktree_fingerprint"] = self.worktree_fingerprint
        return document


@dataclass(frozen=True)
class ExecutionContext:
    repository_root: Path
    case_directory: Path


@dataclass(frozen=True)
class FamilyDefinition:
    name: str
    description: str
    analyzer: str | None = None


@dataclass(frozen=True)
class ExactCaseDefinition:
    name: str
    description: str


@dataclass(frozen=True)
class BuildProfileDefinition:
    name: str
    description: str


JsonObject = Mapping[str, Any]
