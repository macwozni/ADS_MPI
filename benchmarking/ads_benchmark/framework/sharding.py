"""Deterministic, self-validating partitions of frozen plan manifests.

The functions in this module are deliberately free of filesystem and runner
concerns.  A shard is a JSON-compatible envelope around a subset of the
``cases`` array from an immutable ``ads-benchmark-plan`` manifest.  Every
envelope repeats the complete parent metadata and cryptographic identities
needed to prove that a collection of shards reconstructs exactly that parent.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import hashlib
import re
from typing import Iterable, Mapping

from .errors import DuplicateCaseError, ValidationError
from .planner import PLAN_SCHEMA_VERSION, canonical_json


SHARD_SCHEMA_VERSION = 1
SHARD_KIND = "ads-benchmark-plan-shard"
SHARD_BY_INDEX = "index"
SHARD_BY_GROUP = "problem-scheme-degree"

_SHA256 = re.compile(r"^sha256:[0-9a-f]{64}$")
_PLAN_KEYS = {
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
_PARENT_KEYS = _PLAN_KEYS - {"cases"}
_SHARD_KEYS = {
    "schema_version",
    "kind",
    "parent_plan_hash",
    "parent_plan",
    "expected_case_count",
    "expected_cases_digest",
    "strategy",
    "shard_index",
    "shard_count",
    "group",
    "case_count",
    "subset_digest",
    "cases",
}
_CASE_KEYS = {"case_id", "configuration"}
_GROUP_KEYS = {"problem", "scheme", "test_degree", "trial_degree"}


def _digest(value: object) -> str:
    encoded = canonical_json(value).encode("utf-8")
    return "sha256:" + hashlib.sha256(encoded).hexdigest()


def _object(
    value: object,
    *,
    keys: set[str],
    field: str,
) -> Mapping[str, object]:
    if not isinstance(value, dict):
        raise ValidationError(f"{field} must be an object")
    missing = sorted(keys - set(value))
    extra = sorted(set(value) - keys)
    if missing or extra:
        details: list[str] = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if extra:
            details.append("unknown " + ", ".join(extra))
        raise ValidationError(f"{field} has invalid keys: {'; '.join(details)}")
    return value


def _positive_integer(value: object, field: str) -> int:
    if type(value) is not int or value < 1:
        raise ValidationError(f"{field} must be a positive integer")
    return value


def _nonnegative_integer(value: object, field: str) -> int:
    if type(value) is not int or value < 0:
        raise ValidationError(f"{field} must be a nonnegative integer")
    return value


def _digest_string(value: object, field: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        raise ValidationError(f"{field} must be a SHA-256 digest")
    return value


def _case_id(entry: object, field: str) -> str:
    case = _object(entry, keys=_CASE_KEYS, field=field)
    case_id = case["case_id"]
    if not isinstance(case_id, str) or not case_id:
        raise ValidationError(f"{field}.case_id must be a nonempty string")
    if not isinstance(case["configuration"], dict):
        raise ValidationError(f"{field}.configuration must be an object")
    return case_id


def _validated_cases(
    value: object,
    *,
    field: str,
    require_sorted: bool,
) -> list[Mapping[str, object]]:
    if not isinstance(value, list) or not value:
        raise ValidationError(f"{field} must be a nonempty array")
    cases: list[Mapping[str, object]] = []
    identifiers: list[str] = []
    seen: set[str] = set()
    for index, raw_case in enumerate(value):
        case_id = _case_id(raw_case, f"{field}[{index}]")
        if case_id in seen:
            raise DuplicateCaseError(f"{field} duplicates case {case_id}")
        seen.add(case_id)
        identifiers.append(case_id)
        assert isinstance(raw_case, dict)
        cases.append(raw_case)
    if require_sorted and identifiers != sorted(identifiers):
        raise ValidationError(f"{field} is not in deterministic case_id order")
    return cases


def _validate_parent_metadata(plan: Mapping[str, object], field: str) -> int:
    if (
        type(plan["schema_version"]) is not int
        or plan["schema_version"] != PLAN_SCHEMA_VERSION
    ):
        raise ValidationError(
            f"{field}.schema_version must be {PLAN_SCHEMA_VERSION}"
        )
    if plan["kind"] != "ads-benchmark-plan":
        raise ValidationError(f"{field}.kind must be ads-benchmark-plan")
    if not isinstance(plan["run_id"], str) or not plan["run_id"]:
        raise ValidationError(f"{field}.run_id must be a nonempty string")
    _digest_string(plan["config_hash"], f"{field}.config_hash")
    for metadata_field in ("profile", "filters", "repository", "result_layout"):
        if not isinstance(plan[metadata_field], dict):
            raise ValidationError(
                f"{field}.{metadata_field} must be an object"
            )
    if not isinstance(plan["created_at"], str) or not plan["created_at"]:
        raise ValidationError(f"{field}.created_at must be a nonempty string")
    return _positive_integer(plan["case_count"], f"{field}.case_count")


def _validate_parent_manifest(
    manifest: Mapping[str, object],
) -> list[Mapping[str, object]]:
    plan = _object(manifest, keys=_PLAN_KEYS, field="parent plan manifest")
    expected_count = _validate_parent_metadata(plan, "parent plan")

    cases = _validated_cases(
        plan["cases"], field="parent plan cases", require_sorted=True
    )
    if expected_count != len(cases):
        raise ValidationError("parent plan case_count does not match cases")
    configuration_digest = _digest([case["configuration"] for case in cases])
    if configuration_digest != plan["config_hash"]:
        raise ValidationError("parent plan config_hash does not match its cases")
    return cases


def _degree_vector(value: object, field: str) -> tuple[int, int, int]:
    if (
        not isinstance(value, list)
        or len(value) != 3
        or any(type(component) is not int or component < 1 for component in value)
    ):
        raise ValidationError(f"{field} must contain three positive integers")
    return value[0], value[1], value[2]


def _case_group(case: Mapping[str, object], field: str) -> tuple[object, ...]:
    configuration = case["configuration"]
    assert isinstance(configuration, dict)
    problem = configuration.get("problem")
    scheme = configuration.get("scheme")
    spaces = configuration.get("spaces")
    if not isinstance(problem, str) or not problem:
        raise ValidationError(f"{field}.configuration.problem must be a string")
    if not isinstance(scheme, str) or not scheme:
        raise ValidationError(f"{field}.configuration.scheme must be a string")
    if not isinstance(spaces, dict):
        raise ValidationError(f"{field}.configuration.spaces must be an object")
    test = _degree_vector(
        spaces.get("test_degree"), f"{field}.configuration.spaces.test_degree"
    )
    trial = _degree_vector(
        spaces.get("trial_degree"), f"{field}.configuration.spaces.trial_degree"
    )
    return problem, scheme, test, trial


def _group_document(group: tuple[object, ...]) -> dict[str, object]:
    problem, scheme, test, trial = group
    return {
        "problem": problem,
        "scheme": scheme,
        "test_degree": list(test),  # type: ignore[arg-type]
        "trial_degree": list(trial),  # type: ignore[arg-type]
    }


def _parse_group(value: object, field: str) -> tuple[object, ...]:
    group = _object(value, keys=_GROUP_KEYS, field=field)
    problem = group["problem"]
    scheme = group["scheme"]
    if not isinstance(problem, str) or not problem:
        raise ValidationError(f"{field}.problem must be a nonempty string")
    if not isinstance(scheme, str) or not scheme:
        raise ValidationError(f"{field}.scheme must be a nonempty string")
    return (
        problem,
        scheme,
        _degree_vector(group["test_degree"], f"{field}.test_degree"),
        _degree_vector(group["trial_degree"], f"{field}.trial_degree"),
    )


def _normalized_strategy(strategy: str) -> str:
    if strategy == "group":
        return SHARD_BY_GROUP
    if strategy not in {SHARD_BY_INDEX, SHARD_BY_GROUP}:
        raise ValidationError(
            "shard strategy must be 'index' or 'problem-scheme-degree'"
        )
    return strategy


def shard_manifest(
    manifest: Mapping[str, object],
    *,
    strategy: str,
    shard_count: int | None = None,
) -> tuple[dict[str, object], ...]:
    """Partition one complete plan manifest into deterministic envelopes.

    ``index`` uses round-robin assignment over the already sorted ``case_id``
    sequence and requires an explicit ``shard_count``.  The grouping strategy
    creates exactly one shard for every unique
    ``(problem, scheme, test_degree, trial_degree)`` tuple.
    """

    cases = _validate_parent_manifest(manifest)
    normalized = _normalized_strategy(strategy)
    partitions: list[list[Mapping[str, object]]]
    groups: list[tuple[object, ...] | None]

    if normalized == SHARD_BY_INDEX:
        count = _positive_integer(shard_count, "shard_count")
        if count > len(cases):
            raise ValidationError(
                f"shard_count {count} exceeds parent case count {len(cases)}"
            )
        partitions = [[] for _ in range(count)]
        for index, case in enumerate(cases):
            partitions[index % count].append(case)
        groups = [None] * count
    else:
        if shard_count is not None:
            raise ValidationError(
                "problem-scheme-degree sharding derives shard_count; do not set it"
            )
        by_group: dict[tuple[object, ...], list[Mapping[str, object]]] = {}
        for index, case in enumerate(cases):
            group = _case_group(case, f"parent plan cases[{index}]")
            by_group.setdefault(group, []).append(case)
        ordered_groups = sorted(by_group)
        partitions = [by_group[group] for group in ordered_groups]
        groups = list(ordered_groups)
        count = len(partitions)

    parent = {
        key: deepcopy(value) for key, value in manifest.items() if key != "cases"
    }
    parent_hash = _digest(manifest)
    all_cases_digest = _digest(cases)
    shards: list[dict[str, object]] = []
    for index, (partition, group) in enumerate(zip(partitions, groups, strict=True)):
        copied_cases = deepcopy(partition)
        shards.append(
            {
                "schema_version": SHARD_SCHEMA_VERSION,
                "kind": SHARD_KIND,
                "parent_plan_hash": parent_hash,
                "parent_plan": deepcopy(parent),
                "expected_case_count": len(cases),
                "expected_cases_digest": all_cases_digest,
                "strategy": normalized,
                "shard_index": index,
                "shard_count": count,
                "group": _group_document(group) if group is not None else None,
                "case_count": len(copied_cases),
                "subset_digest": _digest(copied_cases),
                "cases": copied_cases,
            }
        )
    return tuple(shards)


@dataclass(frozen=True)
class ValidatedShardSet:
    """A complete compatible shard set, ready for execution/result merging."""

    parent_manifest: dict[str, object]
    strategy: str
    shard_count: int
    case_ids_by_shard: tuple[tuple[str, ...], ...]

    @property
    def case_ids(self) -> tuple[str, ...]:
        return tuple(
            sorted(case_id for shard in self.case_ids_by_shard for case_id in shard)
        )


@dataclass(frozen=True)
class _ParsedShard:
    document: Mapping[str, object]
    parent: Mapping[str, object]
    strategy: str
    index: int
    count: int
    group: tuple[object, ...] | None
    cases: tuple[Mapping[str, object], ...]


def _parse_shard(document: Mapping[str, object], position: int) -> _ParsedShard:
    field = f"shards[{position}]"
    shard = _object(document, keys=_SHARD_KEYS, field=field)
    if shard["schema_version"] != SHARD_SCHEMA_VERSION:
        raise ValidationError(
            f"{field}.schema_version must be {SHARD_SCHEMA_VERSION}"
        )
    if shard["kind"] != SHARD_KIND:
        raise ValidationError(f"{field}.kind must be {SHARD_KIND}")
    _digest_string(shard["parent_plan_hash"], f"{field}.parent_plan_hash")
    _digest_string(shard["expected_cases_digest"], f"{field}.expected_cases_digest")
    _digest_string(shard["subset_digest"], f"{field}.subset_digest")
    parent = _object(
        shard["parent_plan"], keys=_PARENT_KEYS, field=f"{field}.parent_plan"
    )
    parent_case_count = _validate_parent_metadata(parent, f"{field}.parent_plan")

    strategy_value = shard["strategy"]
    if not isinstance(strategy_value, str):
        raise ValidationError(f"{field}.strategy must be a string")
    strategy = _normalized_strategy(strategy_value)
    if strategy_value != strategy:
        raise ValidationError(f"{field}.strategy must use its canonical name")
    index = _nonnegative_integer(shard["shard_index"], f"{field}.shard_index")
    count = _positive_integer(shard["shard_count"], f"{field}.shard_count")
    if index >= count:
        raise ValidationError(f"{field}.shard_index must be below shard_count")

    if strategy == SHARD_BY_INDEX:
        if shard["group"] is not None:
            raise ValidationError(f"{field}.group must be null for index sharding")
        group = None
    else:
        group = _parse_group(shard["group"], f"{field}.group")

    cases = _validated_cases(
        shard["cases"], field=f"{field}.cases", require_sorted=True
    )
    case_count = _positive_integer(shard["case_count"], f"{field}.case_count")
    if case_count != len(cases):
        raise ValidationError(f"{field}.case_count does not match cases")
    if _digest(cases) != shard["subset_digest"]:
        raise ValidationError(f"{field}.subset_digest does not match cases")
    expected_count = _positive_integer(
        shard["expected_case_count"], f"{field}.expected_case_count"
    )
    if parent_case_count != expected_count:
        raise ValidationError(
            f"{field}.expected_case_count conflicts with parent plan"
        )
    if count > expected_count:
        raise ValidationError(f"{field}.shard_count exceeds expected case count")
    if group is not None:
        for case_index, case in enumerate(cases):
            if _case_group(case, f"{field}.cases[{case_index}]") != group:
                raise ValidationError(
                    f"{field} contains a case outside its declared group"
                )
    return _ParsedShard(
        document=shard,
        parent=parent,
        strategy=strategy,
        index=index,
        count=count,
        group=group,
        cases=tuple(cases),
    )


def validate_shard_manifests(
    documents: Iterable[Mapping[str, object]],
) -> ValidatedShardSet:
    """Validate completeness, compatibility, uniqueness, and assignment.

    The helper is intentionally useful before either manifest merging or
    result merging: callers can compare contributed result case IDs with
    ``case_ids_by_shard`` without duplicating shard-set validation.
    """

    materialized = list(documents)
    if not materialized:
        raise ValidationError("cannot merge an empty shard collection")
    parsed = [
        _parse_shard(document, position)
        for position, document in enumerate(materialized)
    ]
    first = parsed[0]
    first_document = first.document
    parent_canonical = canonical_json(first.parent)
    seen_indices: set[int] = set()
    cases_by_id: dict[str, tuple[str, Mapping[str, object], int]] = {}

    for shard in parsed:
        if shard.count != first.count or shard.strategy != first.strategy:
            raise ValidationError("incompatible shard strategy or shard_count")
        if canonical_json(shard.parent) != parent_canonical:
            if shard.parent.get("config_hash") != first.parent.get("config_hash"):
                raise ValidationError("incompatible config hashes across shards")
            if shard.parent.get("repository") != first.parent.get("repository"):
                raise ValidationError(
                    "incompatible repository/source identities across shards"
                )
            raise ValidationError("incompatible parent plan metadata across shards")
        if (
            shard.document["expected_case_count"]
            != first_document["expected_case_count"]
        ):
            raise ValidationError("incompatible expected case totals across shards")
        if shard.document["parent_plan_hash"] != first_document["parent_plan_hash"]:
            raise ValidationError("incompatible parent plan hashes across shards")
        if (
            shard.document["expected_cases_digest"]
            != first_document["expected_cases_digest"]
        ):
            raise ValidationError("incompatible expected all-case digests across shards")
        if shard.index in seen_indices:
            raise ValidationError(f"duplicate shard index {shard.index}")
        seen_indices.add(shard.index)
        for case in shard.cases:
            case_id = str(case["case_id"])
            canonical = canonical_json(case)
            previous = cases_by_id.get(case_id)
            if previous is not None:
                if previous[0] == canonical:
                    raise DuplicateCaseError(
                        f"duplicate case {case_id} in shards "
                        f"{previous[2]} and {shard.index}"
                    )
                raise ValidationError(
                    f"conflicting definitions for case {case_id} in shards "
                    f"{previous[2]} and {shard.index}"
                )
            cases_by_id[case_id] = (canonical, case, shard.index)

    missing_indices = sorted(set(range(first.count)) - seen_indices)
    if missing_indices:
        raise ValidationError(
            "missing shard indices: "
            + ", ".join(str(index) for index in missing_indices)
        )
    expected_count = int(first_document["expected_case_count"])
    if len(cases_by_id) != expected_count:
        qualifier = "missing" if len(cases_by_id) < expected_count else "unexpected"
        raise ValidationError(
            f"{qualifier} cases: expected {expected_count}, found {len(cases_by_id)}"
        )

    ordered_cases = [cases_by_id[case_id][1] for case_id in sorted(cases_by_id)]
    if _digest(ordered_cases) != first_document["expected_cases_digest"]:
        raise ValidationError(
            "merged cases do not match the expected all-case digest"
        )
    reconstructed = deepcopy(dict(first.parent))
    reconstructed["cases"] = deepcopy(ordered_cases)
    if _digest(reconstructed) != first_document["parent_plan_hash"]:
        raise ValidationError("merged shards do not match the frozen parent plan hash")
    # Also re-check the plan-level configuration digest independently.  This
    # produces a precise error when a producer tampers with both shard digests.
    configuration_digest = _digest(
        [case["configuration"] for case in ordered_cases]
    )
    if configuration_digest != first.parent["config_hash"]:
        raise ValidationError("merged cases do not match the parent config hash")

    parsed_by_index = {shard.index: shard for shard in parsed}
    if first.strategy == SHARD_BY_INDEX:
        for global_index, case in enumerate(ordered_cases):
            actual = cases_by_id[str(case["case_id"])][2]
            expected = global_index % first.count
            if actual != expected:
                raise ValidationError(
                    f"case {case['case_id']} belongs to shard {expected}, not {actual}"
                )
    else:
        expected_groups = sorted(
            {
                _case_group(case, f"merged cases[{index}]")
                for index, case in enumerate(ordered_cases)
            }
        )
        if len(expected_groups) != first.count:
            raise ValidationError(
                "group shard_count does not match unique problem/scheme/degree groups"
            )
        for index, expected_group in enumerate(expected_groups):
            shard = parsed_by_index[index]
            if shard.group != expected_group:
                raise ValidationError(
                    f"group shard {index} has the wrong deterministic group assignment"
                )

    case_ids_by_shard = tuple(
        tuple(str(case["case_id"]) for case in parsed_by_index[index].cases)
        for index in range(first.count)
    )
    return ValidatedShardSet(
        parent_manifest=reconstructed,
        strategy=first.strategy,
        shard_count=first.count,
        case_ids_by_shard=case_ids_by_shard,
    )


def merge_shard_manifests(
    documents: Iterable[Mapping[str, object]],
) -> dict[str, object]:
    """Reconstruct the exact parent plan after strict shard-set validation."""

    return validate_shard_manifests(documents).parent_manifest


def shard_subset_manifest(
    document: Mapping[str, object],
    *,
    run_id: str,
) -> dict[str, object]:
    """Return a conventional executable plan for one validated shard.

    The returned manifest intentionally has a subset-specific ``config_hash``
    and run ID so the existing frozen-plan parser can execute it unchanged.
    The original shard envelope must remain next to the run as its proof of
    membership in the complete parent plan.
    """

    shard = _parse_shard(document, 0)
    if not isinstance(run_id, str) or not run_id:
        raise ValidationError("shard run_id must be a nonempty string")
    manifest = deepcopy(dict(shard.parent))
    manifest["run_id"] = run_id
    manifest["case_count"] = len(shard.cases)
    manifest["cases"] = deepcopy(list(shard.cases))
    manifest["config_hash"] = _digest(
        [case["configuration"] for case in shard.cases]
    )
    return manifest


def expected_shard_case_ids(
    documents: Iterable[Mapping[str, object]],
) -> tuple[tuple[str, ...], ...]:
    """Expose validated per-shard coverage for result-merging code."""

    return validate_shard_manifests(documents).case_ids_by_shard
