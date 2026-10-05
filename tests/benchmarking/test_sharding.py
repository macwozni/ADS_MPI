from __future__ import annotations

from copy import deepcopy
import hashlib
from pathlib import Path
import unittest

from ads_benchmark.catalog import build_catalog
from ads_benchmark.framework.config import load_profiles
from ads_benchmark.framework.errors import DuplicateCaseError, ValidationError
from ads_benchmark.framework.model import RepositoryState
from ads_benchmark.framework.planner import (
    Planner,
    canonical_json,
    frozen_plan_from_manifest,
)
from ads_benchmark.framework.sharding import (
    SHARD_BY_GROUP,
    SHARD_BY_INDEX,
    expected_shard_case_ids,
    merge_shard_manifests,
    shard_manifest,
    shard_subset_manifest,
    validate_shard_manifests,
)
from benchmark_paths import CONFIG_DIRECTORY


def digest(value: object) -> str:
    return "sha256:" + hashlib.sha256(
        canonical_json(value).encode("utf-8")
    ).hexdigest()


def reseal_subset(shard: dict[str, object]) -> None:
    cases = shard["cases"]
    assert isinstance(cases, list)
    cases.sort(key=lambda case: case["case_id"])
    shard["case_count"] = len(cases)
    shard["subset_digest"] = digest(cases)


class ManifestShardingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.catalog = build_catalog()
        profiles = load_profiles(CONFIG_DIRECTORY)
        cls.planner = Planner(profiles, cls.catalog)
        cls.repository = RepositoryState(commit="a" * 40, dirty=False)
        cls.small_manifest = cls.planner.plan("strong-scaling-smoke").manifest(
            run_id="shard-parent",
            repository=cls.repository,
            created_at="2026-10-05T12:00:00+00:00",
        )

    def test_index_sharding_is_deterministic_round_robin_and_reversible(self) -> None:
        first = shard_manifest(
            self.small_manifest, strategy=SHARD_BY_INDEX, shard_count=3
        )
        second = shard_manifest(
            deepcopy(self.small_manifest), strategy="index", shard_count=3
        )
        self.assertEqual(first, second)
        self.assertEqual(merge_shard_manifests(reversed(first)), self.small_manifest)

        source_ids = [case["case_id"] for case in self.small_manifest["cases"]]
        expected = tuple(tuple(source_ids[index::3]) for index in range(3))
        self.assertEqual(expected_shard_case_ids(first), expected)
        for index, shard in enumerate(first):
            self.assertEqual(shard["shard_index"], index)
            self.assertEqual(shard["shard_count"], 3)
            self.assertEqual(shard["expected_case_count"], len(source_ids))
            self.assertEqual(
                shard["parent_plan"]["config_hash"],
                self.small_manifest["config_hash"],
            )
            self.assertRegex(shard["parent_plan_hash"], r"^sha256:[0-9a-f]{64}$")
            self.assertRegex(shard["expected_cases_digest"], r"^sha256:[0-9a-f]{64}$")
            self.assertRegex(shard["subset_digest"], r"^sha256:[0-9a-f]{64}$")

    def test_group_sharding_scales_to_full_thousands_case_plan(self) -> None:
        manifest = self.planner.plan("cluster-scaling").manifest(
            run_id="cluster-shards",
            repository=self.repository,
            created_at="2026-10-05T12:00:00+00:00",
        )
        self.assertEqual(manifest["case_count"], 12_960)
        shards = shard_manifest(manifest, strategy="group")
        self.assertEqual(len(shards), 3 * 3 * 15)
        self.assertTrue(all(shard["strategy"] == SHARD_BY_GROUP for shard in shards))
        self.assertTrue(all(shard["case_count"] == 96 for shard in shards))
        for shard in shards:
            group = shard["group"]
            self.assertIsInstance(group, dict)
            cases = shard["cases"]
            self.assertTrue(
                all(
                    case["configuration"]["problem"] == group["problem"]
                    and case["configuration"]["scheme"] == group["scheme"]
                    and case["configuration"]["spaces"]["test_degree"]
                    == group["test_degree"]
                    and case["configuration"]["spaces"]["trial_degree"]
                    == group["trial_degree"]
                    for case in cases
                )
            )
        validated = validate_shard_manifests(reversed(shards))
        self.assertEqual(validated.parent_manifest, manifest)
        self.assertEqual(len(validated.case_ids), 12_960)

    def test_subset_manifest_round_trips_through_existing_frozen_parser(self) -> None:
        shard = shard_manifest(
            self.small_manifest, strategy="index", shard_count=2
        )[0]
        executable = shard_subset_manifest(shard, run_id="shard-parent-0000")
        frozen = frozen_plan_from_manifest(executable, self.catalog)
        self.assertEqual(frozen.run_id, "shard-parent-0000")
        self.assertEqual(
            [case.case_id for case in frozen.cases],
            [case["case_id"] for case in shard["cases"]],
        )
        self.assertNotEqual(executable["config_hash"], self.small_manifest["config_hash"])

    def test_invalid_shard_counts_and_group_options_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValidationError, "positive integer"):
            shard_manifest(self.small_manifest, strategy="index", shard_count=0)
        with self.assertRaisesRegex(ValidationError, "exceeds parent case count"):
            shard_manifest(self.small_manifest, strategy="index", shard_count=999)
        with self.assertRaisesRegex(ValidationError, "derives shard_count"):
            shard_manifest(self.small_manifest, strategy="group", shard_count=2)
        with self.assertRaisesRegex(ValidationError, "shard strategy"):
            shard_manifest(self.small_manifest, strategy="random", shard_count=2)

    def test_missing_and_duplicate_shard_indices_are_rejected(self) -> None:
        shards = list(
            shard_manifest(self.small_manifest, strategy="index", shard_count=3)
        )
        with self.assertRaisesRegex(ValidationError, "missing shard indices: 1"):
            validate_shard_manifests((shards[0], shards[2]))

        duplicate = deepcopy(shards)
        duplicate[1]["shard_index"] = 0
        with self.assertRaisesRegex(ValidationError, "duplicate shard index 0"):
            validate_shard_manifests(duplicate)

    def test_missing_duplicate_and_conflicting_cases_are_rejected(self) -> None:
        original = list(
            shard_manifest(self.small_manifest, strategy="index", shard_count=2)
        )

        missing = deepcopy(original)
        missing[0]["cases"].pop()
        reseal_subset(missing[0])
        with self.assertRaisesRegex(ValidationError, "missing cases"):
            validate_shard_manifests(missing)

        duplicated = deepcopy(original)
        duplicated[0]["cases"].append(deepcopy(duplicated[1]["cases"][0]))
        reseal_subset(duplicated[0])
        with self.assertRaisesRegex(DuplicateCaseError, "duplicate case"):
            validate_shard_manifests(duplicated)

        conflicting = deepcopy(original)
        altered = deepcopy(conflicting[1]["cases"][0])
        altered["configuration"]["mesh"]["elements"][0] += 1
        conflicting[0]["cases"].append(altered)
        reseal_subset(conflicting[0])
        with self.assertRaisesRegex(ValidationError, "conflicting definitions"):
            validate_shard_manifests(conflicting)

    def test_incompatible_config_source_and_parent_identities_are_rejected(self) -> None:
        original = list(
            shard_manifest(self.small_manifest, strategy="index", shard_count=2)
        )

        config = deepcopy(original)
        config[1]["parent_plan"]["config_hash"] = "sha256:" + "0" * 64
        with self.assertRaisesRegex(ValidationError, "incompatible config hashes"):
            validate_shard_manifests(config)

        source = deepcopy(original)
        source[1]["parent_plan"]["repository"]["commit"] = "b" * 40
        with self.assertRaisesRegex(ValidationError, "repository/source identities"):
            validate_shard_manifests(source)

        identity = deepcopy(original)
        identity[1]["parent_plan_hash"] = "sha256:" + "f" * 64
        with self.assertRaisesRegex(ValidationError, "parent plan hashes"):
            validate_shard_manifests(identity)

    def test_digest_tampering_and_wrong_index_assignment_are_rejected(self) -> None:
        original = list(
            shard_manifest(self.small_manifest, strategy="index", shard_count=2)
        )
        tampered = deepcopy(original)
        tampered[0]["cases"][0]["configuration"]["mesh"]["elements"][0] += 1
        with self.assertRaisesRegex(ValidationError, "subset_digest"):
            validate_shard_manifests(tampered)

        misplaced = deepcopy(original)
        first = misplaced[0]["cases"].pop(0)
        second = misplaced[1]["cases"].pop(0)
        misplaced[0]["cases"].append(second)
        misplaced[1]["cases"].append(first)
        reseal_subset(misplaced[0])
        reseal_subset(misplaced[1])
        with self.assertRaisesRegex(ValidationError, "belongs to shard"):
            validate_shard_manifests(misplaced)

    def test_group_assignment_and_complete_parent_metadata_are_enforced(self) -> None:
        shards = list(shard_manifest(self.small_manifest, strategy="group"))
        self.assertGreater(len(shards), 1)
        swapped = deepcopy(shards)
        swapped[0]["shard_index"], swapped[1]["shard_index"] = (
            swapped[1]["shard_index"],
            swapped[0]["shard_index"],
        )
        with self.assertRaisesRegex(ValidationError, "wrong deterministic group"):
            validate_shard_manifests(swapped)

        incomplete = deepcopy(shards)
        del incomplete[0]["parent_plan"]["filters"]
        with self.assertRaisesRegex(ValidationError, "missing filters"):
            validate_shard_manifests(incomplete)


if __name__ == "__main__":
    unittest.main()
