from __future__ import annotations

from pathlib import Path
import tempfile
import threading
import unittest
from unittest.mock import patch

from ads_benchmark.framework.errors import StorageError
from ads_benchmark.framework.storage import ResultStore


class ShardStorageTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.store = ResultStore(self.root)
        self.run_id = "parent-plan"
        self.store.create_run(
            self.run_id,
            {
                "schema_version": 1,
                "kind": "ads-benchmark-plan",
                "run_id": self.run_id,
            },
        )

    @staticmethod
    def document(index: int) -> dict[str, object]:
        return {
            "schema_version": 1,
            "kind": "ads-benchmark-plan-shard",
            "shard_index": index,
        }

    def test_indexed_shards_round_trip_without_overwrite(self) -> None:
        first = self.document(0)
        second = self.document(1)
        first_path = self.store.write_shard_manifest(self.run_id, 0, first)
        second_path = self.store.write_shard_manifest(self.run_id, 1, second)
        self.assertEqual(first_path.name, "shard-000000.json")
        self.assertEqual(second_path.name, "shard-000001.json")
        self.assertEqual(self.store.read_shard_manifest(self.run_id, 0), first)
        self.assertEqual(
            self.store.read_all_shard_manifests(self.run_id),
            (first, second),
        )
        with self.assertRaisesRegex(StorageError, "refusing overwrite"):
            self.store.write_shard_manifest(self.run_id, 0, first)

    def test_identity_bounds_and_symlinks_are_rejected(self) -> None:
        with self.assertRaisesRegex(StorageError, "identity"):
            self.store.write_shard_manifest(self.run_id, 0, self.document(1))
        with self.assertRaisesRegex(StorageError, "between 0 and 999999"):
            self.store.read_shard_manifest(self.run_id, -1)

        shards = self.root / "benchmarks" / self.run_id / "shards"
        shards.mkdir()
        outside = self.root / "outside.json"
        outside.write_text("{}\n", encoding="utf-8")
        (shards / "shard-000000.json").symlink_to(outside)
        with self.assertRaises(StorageError):
            self.store.read_all_shard_manifests(self.run_id)

    def test_bulk_generation_lock_refuses_a_concurrent_writer(self) -> None:
        documents = (self.document(0), self.document(1))
        entered_first_write = threading.Event()
        release_first_write = threading.Event()
        worker_errors: list[BaseException] = []
        original_write = self.store._write_shard_manifest_at

        def delayed_first_write(run_fd, run_id, index, document):
            if index == 0:
                entered_first_write.set()
                if not release_first_write.wait(timeout=5):
                    raise RuntimeError("test timed out waiting to release writer")
            return original_write(run_fd, run_id, index, document)

        def write_complete_set() -> None:
            try:
                self.store.write_shard_manifests(self.run_id, documents)
            except BaseException as error:  # pragma: no cover - asserted below
                worker_errors.append(error)

        with patch.object(
            self.store,
            "_write_shard_manifest_at",
            side_effect=delayed_first_write,
        ):
            worker = threading.Thread(target=write_complete_set)
            worker.start()
            self.assertTrue(entered_first_write.wait(timeout=5))
            try:
                with self.assertRaisesRegex(StorageError, "currently in use"):
                    self.store.write_shard_manifests(self.run_id, documents)
            finally:
                release_first_write.set()
                worker.join(timeout=5)

        self.assertFalse(worker.is_alive())
        self.assertEqual(worker_errors, [])
        self.assertEqual(
            self.store.read_all_shard_manifests(self.run_id), documents
        )


if __name__ == "__main__":
    unittest.main()
