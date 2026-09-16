"""Immutable public data model shared by every benchmark family."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping


Vector3 = tuple[int, int, int]


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

    def to_dict(self) -> dict[str, object]:
        return {
            "warmups": self.warmups,
            "samples": self.samples,
            "timeout_seconds": self.timeout_seconds,
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
    measurement: MeasurementSpec
    build_profile: str
    launcher: str

    def to_dict(self) -> dict[str, object]:
        return {
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
            "openmp": {"threads": self.openmp_threads},
            "measurement": self.measurement.to_dict(),
            "build": {"profile": self.build_profile},
            "launcher": self.launcher,
        }


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

    def to_dict(self) -> dict[str, object]:
        return {"commit": self.commit, "dirty": self.dirty}


@dataclass(frozen=True)
class ExecutionContext:
    repository_root: Path
    case_directory: Path


@dataclass(frozen=True)
class FamilyDefinition:
    name: str
    description: str


@dataclass(frozen=True)
class ExactCaseDefinition:
    name: str
    description: str


@dataclass(frozen=True)
class BuildProfileDefinition:
    name: str
    description: str


JsonObject = Mapping[str, Any]
