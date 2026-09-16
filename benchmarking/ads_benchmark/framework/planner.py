"""Deterministic matrix expansion, case identity, filtering, and manifests."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
from itertools import product
import json
from typing import Iterable

from .config import ProfileDefinition
from .errors import DuplicateCaseError, ValidationError
from .filtering import CaseFilters
from .model import CaseSpec, PlannedCase, RepositoryState
from .registry import Catalog, Registry
from .validation import validate_case


PLAN_SCHEMA_VERSION = 1


def canonical_json(value: object) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )


def case_identity(spec: CaseSpec) -> tuple[str, str]:
    canonical = canonical_json(spec.to_dict())
    digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:20]
    return f"{spec.problem}-{spec.scheme}-{digest}", canonical


@dataclass(frozen=True)
class Plan:
    profile: ProfileDefinition
    cases: tuple[PlannedCase, ...]
    filters: CaseFilters
    config_hash: str

    def manifest(
        self,
        *,
        run_id: str,
        repository: RepositoryState,
        created_at: str | None = None,
    ) -> dict[str, object]:
        timestamp = created_at or datetime.now(timezone.utc).isoformat()
        return {
            "schema_version": PLAN_SCHEMA_VERSION,
            "kind": "ads-benchmark-plan",
            "run_id": run_id,
            "created_at": timestamp,
            "profile": {
                "name": self.profile.name,
                "description": self.profile.description,
            },
            "filters": self.filters.to_dict(),
            "config_hash": f"sha256:{self.config_hash}",
            "repository": repository.to_dict(),
            "case_count": len(self.cases),
            "cases": [case.to_dict() for case in self.cases],
            "result_layout": {
                "case_directory": "cases/<case-id>",
                "status": "cases/<case-id>/status.json",
                "result": "cases/<case-id>/result.json",
                "stdout": "cases/<case-id>/stdout.log",
                "stderr": "cases/<case-id>/stderr.log",
            },
        }


class Planner:
    """Expand registered profile data without problem-name branches."""

    def __init__(
        self, profiles: Registry[ProfileDefinition], catalog: Catalog
    ) -> None:
        self.profiles = profiles
        self.catalog = catalog

    def _expanded_specs(self, profile: ProfileDefinition) -> Iterable[CaseSpec]:
        axes = product(
            profile.exact_cases,
            profile.problems,
            profile.schemes,
            profile.time_discretizations,
            profile.meshes,
            profile.degree_pairs,
            profile.process_layouts,
            profile.thread_counts,
            profile.build_profiles,
        )
        for (
            exact_case,
            problem_name,
            scheme,
            time_spec,
            mesh,
            degrees,
            mpi,
            threads,
            build_profile,
        ) in axes:
            test_degree, trial_degree = degrees
            yield CaseSpec(
                family=profile.family,
                problem=problem_name,
                scheme=scheme,
                exact_case=exact_case,
                time=time_spec,
                mesh=mesh,
                test_degree=test_degree,
                trial_degree=trial_degree,
                mpi=mpi,
                openmp_threads=threads,
                measurement=profile.measurement,
                build_profile=build_profile,
                launcher=profile.launcher,
            )

    def plan(
        self, profile_name: str, filters: CaseFilters | None = None
    ) -> Plan:
        profile = self.profiles.get(profile_name)
        active_filters = filters or CaseFilters()
        by_canonical: dict[str, str] = {}
        by_id: dict[str, str] = {}
        expanded: list[PlannedCase] = []

        for spec in self._expanded_specs(profile):
            validate_case(spec, self.catalog)
            case_id, canonical = case_identity(spec)
            if canonical in by_canonical:
                raise DuplicateCaseError(
                    f"profile {profile.name} expands duplicate case {case_id}"
                )
            if case_id in by_id and by_id[case_id] != canonical:
                raise DuplicateCaseError(f"case_id hash collision: {case_id}")
            by_canonical[canonical] = case_id
            by_id[case_id] = canonical
            expanded.append(PlannedCase(case_id=case_id, spec=spec))

        expanded.sort(key=lambda item: item.case_id)
        selected = tuple(case for case in expanded if active_filters.matches(case))
        if not selected:
            raise ValidationError(
                f"filters selected no cases from profile {profile.name}"
            )

        hash_input = [case.spec.to_dict() for case in selected]
        config_hash = hashlib.sha256(
            canonical_json(hash_input).encode("utf-8")
        ).hexdigest()
        return Plan(
            profile=profile,
            cases=selected,
            filters=active_filters,
            config_hash=config_hash,
        )
