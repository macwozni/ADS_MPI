"""Deterministic matrix expansion, case identity, filtering, and manifests."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
from itertools import product
import json
import re
from typing import Iterable, Mapping

from .config import ProfileDefinition
from .errors import DuplicateCaseError, ValidationError
from .filtering import CaseFilters
from .model import (
    CaseSpec,
    MeasurementSpec,
    MpiSpec,
    PlannedCase,
    RepositoryState,
    SamplingSpec,
    TimeSpec,
)
from .registry import Catalog, Registry
from .validation import validate_case


PLAN_SCHEMA_VERSION = 1
_SHA256 = re.compile(r"^sha256:([0-9a-f]{64})$")
_COMMIT_SHA = re.compile(r"^[0-9a-f]{40}(?:[0-9a-f]{24})?$")
_MANIFEST_KEYS = {
    "schema_version",
    "kind",
    "run_id",
    "created_at",
    "profile",
    "filters",
    "config_hash",
    "repository",
    "case_count",
    "cases",
    "result_layout",
}
_RESULT_LAYOUT = {
    "case_directory": "cases/<case-id>",
    "status": "cases/<case-id>/status.json",
    "result": "cases/<case-id>/result.json",
    "stdout": "cases/<case-id>/stdout.log",
    "stderr": "cases/<case-id>/stderr.log",
}
_FILTER_KEYS = {
    "problems",
    "schemes",
    "degree_pairs",
    "meshes",
    "mpi_grids",
    "mpi_ranks",
    "openmp_threads",
}


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
                **_RESULT_LAYOUT,
            },
        }


@dataclass(frozen=True)
class FrozenPlan:
    """Validated executable content read back from an immutable manifest."""

    run_id: str
    profile_name: str
    profile_description: str
    filters: Mapping[str, object]
    config_hash: str
    repository: RepositoryState
    cases: tuple[PlannedCase, ...]


def _object(value: object, keys: set[str], field: str) -> Mapping[str, object]:
    if not isinstance(value, dict):
        raise ValidationError(f"frozen manifest {field} must be an object")
    missing = sorted(keys - set(value))
    extra = sorted(set(value) - keys)
    if missing or extra:
        details = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if extra:
            details.append("unknown " + ", ".join(extra))
        raise ValidationError(
            f"frozen manifest {field} has invalid keys: {'; '.join(details)}"
        )
    return value


def _filters_object(value: object) -> Mapping[str, object]:
    if not isinstance(value, dict):
        raise ValidationError("frozen manifest filters must be an object")
    keys = set(value)
    missing = sorted(_FILTER_KEYS - keys)
    extra = sorted(keys - (_FILTER_KEYS | {"steps"}))
    if missing or extra:
        details = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if extra:
            details.append("unknown " + ", ".join(extra))
        raise ValidationError(
            "frozen manifest filters has invalid keys: " + "; ".join(details)
        )
    return value


def _repository_object(value: object) -> Mapping[str, object]:
    if not isinstance(value, dict):
        raise ValidationError("frozen manifest repository must be an object")
    required = {"commit", "dirty"}
    allowed = required | {"worktree_fingerprint"}
    missing = sorted(required - set(value))
    extra = sorted(set(value) - allowed)
    if missing or extra:
        details = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if extra:
            details.append("unknown " + ", ".join(extra))
        raise ValidationError(
            "frozen manifest repository has invalid keys: " + "; ".join(details)
        )
    return value


def _string(value: object, field: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValidationError(
            f"frozen manifest {field} must be a nonempty string"
        )
    return value


def _integer(value: object, field: str, *, minimum: int = 1) -> int:
    if type(value) is not int or value < minimum:
        raise ValidationError(
            f"frozen manifest {field} must be an integer >= {minimum}"
        )
    return value


def _boolean(value: object, field: str) -> bool:
    if type(value) is not bool:
        raise ValidationError(f"frozen manifest {field} must be a boolean")
    return value


def _vector3(value: object, field: str) -> tuple[int, int, int]:
    if not isinstance(value, list) or len(value) != 3:
        raise ValidationError(
            f"frozen manifest {field} must contain exactly three integers"
        )
    parsed = tuple(
        _integer(entry, f"{field}[{index}]") for index, entry in enumerate(value)
    )
    return parsed  # type: ignore[return-value]


def _case_spec(document: object, field: str) -> CaseSpec:
    case = _object(
        document,
        {
            "family",
            "problem",
            "scheme",
            "exact_case",
            "time",
            "mesh",
            "spaces",
            "mpi",
            "openmp",
            "sampling",
            "measurement",
            "build",
            "launcher",
        },
        field,
    )
    time = _object(case["time"], {"final_time", "time_step", "steps"}, f"{field}.time")
    mesh = _object(case["mesh"], {"elements"}, f"{field}.mesh")
    spaces = _object(
        case["spaces"], {"test_degree", "trial_degree"}, f"{field}.spaces"
    )
    mpi = _object(case["mpi"], {"ranks", "process_grid"}, f"{field}.mpi")
    openmp = _object(case["openmp"], {"threads"}, f"{field}.openmp")
    sampling = _object(
        case["sampling"],
        {"points_per_axis", "write_samples"},
        f"{field}.sampling",
    )
    measurement = _object(
        case["measurement"],
        {"warmups", "samples", "timeout_seconds"},
        f"{field}.measurement",
    )
    build = _object(case["build"], {"profile"}, f"{field}.build")
    return CaseSpec(
        family=_string(case["family"], f"{field}.family"),
        problem=_string(case["problem"], f"{field}.problem"),
        scheme=_string(case["scheme"], f"{field}.scheme"),
        exact_case=_string(case["exact_case"], f"{field}.exact_case"),
        time=TimeSpec(
            final_time=_string(time["final_time"], f"{field}.time.final_time"),
            time_step=_string(time["time_step"], f"{field}.time.time_step"),
            steps=_integer(time["steps"], f"{field}.time.steps"),
        ),
        mesh=_vector3(mesh["elements"], f"{field}.mesh.elements"),
        test_degree=_vector3(
            spaces["test_degree"], f"{field}.spaces.test_degree"
        ),
        trial_degree=_vector3(
            spaces["trial_degree"], f"{field}.spaces.trial_degree"
        ),
        mpi=MpiSpec(
            ranks=_integer(mpi["ranks"], f"{field}.mpi.ranks"),
            process_grid=_vector3(
                mpi["process_grid"], f"{field}.mpi.process_grid"
            ),
        ),
        openmp_threads=_integer(openmp["threads"], f"{field}.openmp.threads"),
        sampling=SamplingSpec(
            points_per_axis=_integer(
                sampling["points_per_axis"],
                f"{field}.sampling.points_per_axis",
            ),
            write_samples=_boolean(
                sampling["write_samples"], f"{field}.sampling.write_samples"
            ),
        ),
        measurement=MeasurementSpec(
            warmups=_integer(
                measurement["warmups"], f"{field}.measurement.warmups", minimum=0
            ),
            samples=_integer(
                measurement["samples"], f"{field}.measurement.samples"
            ),
            timeout_seconds=_string(
                measurement["timeout_seconds"],
                f"{field}.measurement.timeout_seconds",
            ),
        ),
        build_profile=_string(build["profile"], f"{field}.build.profile"),
        launcher=_string(case["launcher"], f"{field}.launcher"),
    )


def frozen_plan_from_manifest(
    document: Mapping[str, object], catalog: Catalog
) -> FrozenPlan:
    """Strictly validate and deserialize the only plan a runner may execute."""

    manifest = _object(document, _MANIFEST_KEYS, "root")
    if (
        type(manifest["schema_version"]) is not int
        or manifest["schema_version"] != PLAN_SCHEMA_VERSION
    ):
        raise ValidationError(
            f"frozen manifest schema_version must be {PLAN_SCHEMA_VERSION}"
        )
    if manifest["kind"] != "ads-benchmark-plan":
        raise ValidationError("frozen manifest kind is not ads-benchmark-plan")
    run_id = _string(manifest["run_id"], "run_id")
    _string(manifest["created_at"], "created_at")
    profile = _object(manifest["profile"], {"name", "description"}, "profile")
    filters = _filters_object(manifest["filters"])
    repository = _repository_object(manifest["repository"])
    commit = _string(repository["commit"], "repository.commit")
    if not _COMMIT_SHA.fullmatch(commit):
        raise ValidationError("frozen manifest repository.commit is not a Git SHA")
    fingerprint: str | None = None
    if "worktree_fingerprint" in repository:
        fingerprint = _string(
            repository["worktree_fingerprint"],
            "repository.worktree_fingerprint",
        )
        if _SHA256.fullmatch(fingerprint) is None:
            raise ValidationError(
                "frozen manifest repository.worktree_fingerprint is not a SHA-256 digest"
            )
    repository_state = RepositoryState(
        commit=commit,
        dirty=_boolean(repository["dirty"], "repository.dirty"),
        worktree_fingerprint=fingerprint,
    )

    config_hash = _string(manifest["config_hash"], "config_hash")
    hash_match = _SHA256.fullmatch(config_hash)
    if hash_match is None:
        raise ValidationError("frozen manifest config_hash is not a SHA-256 digest")

    raw_cases = manifest["cases"]
    if not isinstance(raw_cases, list) or not raw_cases:
        raise ValidationError("frozen manifest cases must be a nonempty array")
    cases: list[PlannedCase] = []
    seen: set[str] = set()
    for index, raw_case in enumerate(raw_cases):
        entry = _object(
            raw_case, {"case_id", "configuration"}, f"cases[{index}]"
        )
        case_id = _string(entry["case_id"], f"cases[{index}].case_id")
        spec = _case_spec(entry["configuration"], f"cases[{index}].configuration")
        validate_case(spec, catalog)
        expected_id, _ = case_identity(spec)
        if case_id != expected_id:
            raise ValidationError(
                f"frozen manifest case_id does not match configuration: {case_id}"
            )
        if case_id in seen:
            raise DuplicateCaseError(f"frozen manifest duplicates case {case_id}")
        seen.add(case_id)
        cases.append(PlannedCase(case_id=case_id, spec=spec))

    if (
        type(manifest["case_count"]) is not int
        or manifest["case_count"] != len(cases)
    ):
        raise ValidationError("frozen manifest case_count does not match cases")
    if [case.case_id for case in cases] != sorted(case.case_id for case in cases):
        raise ValidationError("frozen manifest cases are not in deterministic order")
    calculated_hash = hashlib.sha256(
        canonical_json([case.spec.to_dict() for case in cases]).encode("utf-8")
    ).hexdigest()
    if hash_match.group(1) != calculated_hash:
        raise ValidationError("frozen manifest config_hash does not match its cases")
    if manifest["result_layout"] != _RESULT_LAYOUT:
        raise ValidationError("frozen manifest result_layout is incompatible")

    return FrozenPlan(
        run_id=run_id,
        profile_name=_string(profile["name"], "profile.name"),
        profile_description=_string(profile["description"], "profile.description"),
        filters=dict(filters),
        config_hash=calculated_hash,
        repository=repository_state,
        cases=tuple(cases),
    )


def validate_resume_request(
    frozen: FrozenPlan,
    requested: Plan,
    repository: RepositoryState,
) -> None:
    """Refuse to append when either source revision or requested plan changed."""

    frozen_filters = dict(frozen.filters)
    requested_filters = requested.filters.to_dict()
    # ``steps`` was added as a selector after the first manifest schema.  An
    # absent value and an explicit empty selection carry the same semantics.
    if "steps" in frozen_filters or "steps" in requested_filters:
        frozen_filters.setdefault("steps", [])
        requested_filters.setdefault("steps", [])

    if frozen.repository.commit != repository.commit:
        raise ValidationError(
            "resume repository SHA mismatch: "
            f"manifest={frozen.repository.commit} current={repository.commit}"
        )
    if frozen.repository.dirty != repository.dirty:
        raise ValidationError(
            "resume repository dirty-state mismatch with frozen manifest"
        )
    if frozen.repository.worktree_fingerprint is None:
        if frozen.repository.dirty:
            raise ValidationError(
                "resume cannot verify legacy dirty manifest without worktree fingerprint"
            )
    elif (
        repository.worktree_fingerprint is None
        or frozen.repository.worktree_fingerprint
        != repository.worktree_fingerprint
    ):
        raise ValidationError(
            "resume worktree fingerprint mismatch with frozen manifest"
        )
    if (
        frozen.profile_name != requested.profile.name
        or frozen.profile_description != requested.profile.description
        or frozen_filters != requested_filters
        or frozen.config_hash != requested.config_hash
        or [case.to_dict() for case in frozen.cases]
        != [case.to_dict() for case in requested.cases]
    ):
        raise ValidationError(
            "resume configuration mismatch with frozen manifest"
        )


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
                sampling=profile.sampling,
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
