"""Typed, conjunctive filters for expanded benchmark plans."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from .errors import ValidationError
from .model import PlannedCase, Vector3


DegreePair = tuple[Vector3, Vector3]


def parse_vector3(text: str, field: str) -> Vector3:
    separator = "x" if "x" in text.lower() else ","
    fields = text.lower().split(separator)
    if len(fields) != 3:
        raise ValidationError(f"{field} must use AxBxC or A,B,C syntax: {text}")
    try:
        values = tuple(int(value) for value in fields)
    except ValueError as error:
        raise ValidationError(f"{field} must contain integers: {text}") from error
    if any(value <= 0 for value in values):
        raise ValidationError(f"{field} values must be positive: {text}")
    return values  # type: ignore[return-value]


def parse_degree_pair(text: str) -> DegreePair:
    if ":" not in text:
        raise ValidationError(
            f"degree pair must use TEST:TRIAL syntax: {text}"
        )
    test_text, trial_text = text.split(":", 1)

    def degree_vector(value: str, field: str) -> Vector3:
        if "x" in value.lower() or "," in value:
            return parse_vector3(value, field)
        try:
            degree = int(value)
        except ValueError as error:
            raise ValidationError(f"{field} must be an integer or vector") from error
        if degree <= 0:
            raise ValidationError(f"{field} must be positive")
        return (degree, degree, degree)

    return degree_vector(test_text, "test degree"), degree_vector(
        trial_text, "trial degree"
    )


@dataclass(frozen=True)
class CaseFilters:
    problems: frozenset[str] = frozenset()
    schemes: frozenset[str] = frozenset()
    degree_pairs: frozenset[DegreePair] = frozenset()
    meshes: frozenset[Vector3] = frozenset()
    mpi_grids: frozenset[Vector3] = frozenset()
    mpi_ranks: frozenset[int] = frozenset()
    openmp_threads: frozenset[int] = frozenset()

    def matches(self, case: PlannedCase) -> bool:
        spec = case.spec
        checks = (
            not self.problems or spec.problem in self.problems,
            not self.schemes or spec.scheme in self.schemes,
            not self.degree_pairs
            or (spec.test_degree, spec.trial_degree) in self.degree_pairs,
            not self.meshes or spec.mesh in self.meshes,
            not self.mpi_grids or spec.mpi.process_grid in self.mpi_grids,
            not self.mpi_ranks or spec.mpi.ranks in self.mpi_ranks,
            not self.openmp_threads or spec.openmp_threads in self.openmp_threads,
        )
        return all(checks)

    def to_dict(self) -> dict[str, object]:
        return {
            "problems": sorted(self.problems),
            "schemes": sorted(self.schemes),
            "degree_pairs": [
                {"test": list(test), "trial": list(trial)}
                for test, trial in sorted(self.degree_pairs)
            ],
            "meshes": [list(value) for value in sorted(self.meshes)],
            "mpi_grids": [list(value) for value in sorted(self.mpi_grids)],
            "mpi_ranks": sorted(self.mpi_ranks),
            "openmp_threads": sorted(self.openmp_threads),
        }


def positive_integer_set(values: Iterable[int], field: str) -> frozenset[int]:
    result = frozenset(values)
    if any(type(value) is not int or value <= 0 for value in result):
        raise ValidationError(f"{field} filters must be positive integers")
    return result
