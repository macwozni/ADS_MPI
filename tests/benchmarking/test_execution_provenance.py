from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

from ads_benchmark.components.planning import MpiLauncher, describe_launcher
from ads_benchmark.framework.errors import ProvenanceError, StorageError
from ads_benchmark.framework.provenance import (
    EXECUTION_KIND,
    EXECUTION_SCHEMA_VERSION,
    collect_execution_provenance,
    comparison_environment_identity,
    manifest_digest,
    validate_execution_provenance,
)
from ads_benchmark.framework.storage import ResultStore


def manifest(run_id: str = "execution-contract") -> dict[str, object]:
    return {
        "schema_version": 1,
        "kind": "ads-benchmark-plan",
        "run_id": run_id,
        "created_at": "2026-01-01T00:00:00+00:00",
        "profile": {"name": "fake", "description": "fake profile"},
        "filters": {},
        "config_hash": "sha256:" + "a" * 64,
        "repository": {"commit": "0" * 40, "dirty": False},
        "case_count": 2,
        "cases": [
            {
                "case_id": "fake-a",
                "configuration": {
                    "problem": "fake",
                    "build": {"profile": "release"},
                    "mpi": {"ranks": 1, "process_grid": [1, 1, 1]},
                    "openmp": {"threads": 1},
                },
            },
            {
                "case_id": "fake-b",
                "configuration": {
                    "problem": "fake",
                    "build": {"profile": "release"},
                    "mpi": {"ranks": 2, "process_grid": [2, 1, 1]},
                    "openmp": {
                        "threads": 4,
                        "dynamic": False,
                        "proc_bind": "spread",
                        "places": "cores",
                    },
                },
            },
        ],
        "result_layout": {},
    }


class ExecutionProvenanceTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory(prefix="ads-execution-provenance-")
        self.addCleanup(temporary.cleanup)
        self.repository = Path(temporary.name) / "repository"
        self.repository.mkdir()
        profile = self.repository / "benchmarking" / "build" / "release"
        self.profile = profile
        self.compile_flags = (
            f"-J{profile / 'benchmark_OBJ' / 'fake'} "
            f"-I{profile / 'benchmark_OBJ' / 'fake'} -I{profile}"
        )
        executable = profile / "EXEC" / "fake_manufactured"
        executable.parent.mkdir(parents=True)
        executable.write_bytes(b"fake executable version one\n")
        self.executable = executable

        mumps = self.repository / "dependencies" / "mumps"
        mumps_library = mumps / "lib" / "libdmumps.a"
        mumps_library.parent.mkdir(parents=True)
        mumps_library.write_bytes(b"fake mumps archive\n")
        header = mumps / "include" / "dmumps_c.h"
        header.parent.mkdir(parents=True)
        header.write_text('#define MUMPS_VERSION "5.8.test"\n', encoding="utf-8")
        self.mumps_library = mumps_library
        self.mumps = mumps

        core_stamp = profile / "_OBJ" / ".ads-core-build"
        core_stamp.parent.mkdir(parents=True)
        core_stamp.write_text(
            f"FF={sys.executable} -O3\nFFLAGS=-fopenmp\nAR=ar\nARFLAGS=rcs\n",
            encoding="utf-8",
        )
        adapter_stamp = (
            profile
            / "benchmark_OBJ"
            / "fake"
            / ".ads-benchmark-adapter-build"
        )
        adapter_stamp.parent.mkdir(parents=True)
        adapter_stamp.write_text(
            "\n".join(
                (
                    "CONFIG=/tmp/fake-config",
                    "BUILD=release",
                    "SOURCES=fake.F90",
                    f"FF={sys.executable} -O3 -fopenmp -I{mumps / 'include'}",
                    f"FFLAGS={self.compile_flags}",
                    f"USER_LIB={mumps_library}",
                    "",
                )
            ),
            encoding="utf-8",
        )
        self.adapter_stamp = adapter_stamp
        launcher = MpiLauncher(
            executable=sys.executable,
            rank_flag="-n",
            available_mpi_slots=2,
            available_cpu_slots=8,
            extra_args=("-B",),
        )
        self.launcher_description = describe_launcher(launcher)
        self.environment = {
            "PATH": os.environ.get("PATH", ""),
            "MUMPS_DIR": str(mumps),
            "OMP_WAIT_POLICY": "PASSIVE",
            "KMP_BLOCKTIME": "0",
        }
        self.instant = datetime(2026, 10, 5, 12, 30, tzinfo=timezone.utc)

    def collect(self, source: dict[str, object] | None = None) -> dict[str, object]:
        return collect_execution_provenance(
            self.repository,
            source or manifest(),
            launcher_descriptions=[self.launcher_description],
            environment=self.environment,
            recorded_at=self.instant,
        )

    def test_collection_is_complete_explicit_and_self_validating(self) -> None:
        source = manifest()
        record = self.collect(source)
        self.assertEqual(record["schema_version"], EXECUTION_SCHEMA_VERSION)
        self.assertEqual(record["kind"], EXECUTION_KIND)
        self.assertEqual(record["manifest_hash"], manifest_digest(source))
        self.assertEqual(record["config_hash"], source["config_hash"])
        self.assertRegex(record["build_identity_hash"], r"^sha256:[0-9a-f]{64}$")
        self.assertRegex(
            record["execution_compatibility_hash"], r"^sha256:[0-9a-f]{64}$"
        )
        self.assertRegex(record["record_hash"], r"^sha256:[0-9a-f]{64}$")

        build = record["build"]
        self.assertEqual(build["profiles"], ["release"])
        unit = build["units"][0]
        self.assertEqual(unit["profile"], "release")
        self.assertEqual(unit["problem"], "fake")
        self.assertEqual(unit["executable"]["size_bytes"], self.executable.stat().st_size)
        self.assertRegex(unit["executable"]["sha256"], r"^sha256:[0-9a-f]{64}$")
        self.assertIn("-O3", unit["compiler_command"]["value"])
        self.assertEqual(unit["compile_flags"]["value"], self.compile_flags)
        self.assertEqual(unit["link_flags"]["value"], str(self.mumps_library))
        self.assertTrue(build["compiler_versions"]["value"])
        self.assertIn(
            "Python",
            build["compiler_versions"]["value"][0]["version"]["value"],
        )
        self.assertEqual(
            build["mumps_version"]["value"]["version"], "5.8.test"
        )
        self.assertEqual(
            build["configured_libraries"]["mumps"]["value"],
            sorted({str(self.mumps), str(self.mumps_library)}),
        )
        self.assertIsNone(build["configured_libraries"]["lapack"]["value"])
        self.assertTrue(build["configured_libraries"]["lapack"]["reason"])

        self.assertEqual(
            record["topology"]["rank_grids"],
            [
                {"ranks": 1, "process_grid": [1, 1, 1]},
                {"ranks": 2, "process_grid": [2, 1, 1]},
            ],
        )
        self.assertEqual(record["environment"]["OMP_WAIT_POLICY"]["value"], "PASSIVE")
        self.assertIsNone(record["environment"]["OMP_PLACES"]["value"])
        self.assertTrue(record["environment"]["OMP_PLACES"]["reason"])
        for field in (
            "hostname",
            "operating_system",
            "kernel",
            "architecture",
            "cpu_model",
            "physical_cores",
            "logical_cores",
            "memory_bytes",
        ):
            observation = record["system"][field]
            self.assertEqual(set(observation), {"value", "reason"})
            self.assertNotEqual(observation["value"] is None, observation["reason"] is None)
        validate_execution_provenance(record, manifest=source)

    def test_build_and_compatibility_hashes_track_binary_and_flags_not_time(self) -> None:
        first = self.collect()
        later = collect_execution_provenance(
            self.repository,
            manifest(),
            launcher_descriptions=[self.launcher_description],
            environment=self.environment,
            recorded_at=datetime(2026, 10, 6, 12, 30, tzinfo=timezone.utc),
        )
        self.assertEqual(first["build_identity_hash"], later["build_identity_hash"])
        self.assertEqual(
            first["execution_compatibility_hash"],
            later["execution_compatibility_hash"],
        )
        self.assertNotEqual(first["record_hash"], later["record_hash"])

        self.executable.write_bytes(b"fake executable version two\n")
        changed_binary = self.collect()
        self.assertNotEqual(
            first["build_identity_hash"], changed_binary["build_identity_hash"]
        )
        self.assertNotEqual(
            first["execution_compatibility_hash"],
            changed_binary["execution_compatibility_hash"],
        )

        self.adapter_stamp.write_text(
            self.adapter_stamp.read_text(encoding="utf-8").replace(
                f"FFLAGS={self.compile_flags}",
                f"FFLAGS={self.compile_flags} -fcheck=all",
            ),
            encoding="utf-8",
        )
        changed_flags = self.collect()
        self.assertNotEqual(
            changed_binary["build_identity_hash"],
            changed_flags["build_identity_hash"],
        )

    def test_missing_values_are_explicit_and_tampering_is_rejected(self) -> None:
        record = collect_execution_provenance(
            self.repository,
            manifest(),
            environment={"PATH": os.environ.get("PATH", "")},
            recorded_at=self.instant,
        )
        self.assertIsNone(record["launchers"]["value"])
        self.assertTrue(record["launchers"]["reason"])
        tampered = deepcopy(record)
        tampered["build"]["units"][0]["compile_flags"]["value"] = "forged"
        with self.assertRaisesRegex(ProvenanceError, "build_identity_hash"):
            validate_execution_provenance(tampered, manifest=manifest())
        unknown = deepcopy(record)
        unknown["unexpected"] = True
        with self.assertRaisesRegex(ProvenanceError, "unknown unexpected"):
            validate_execution_provenance(unknown, manifest=manifest())

    def test_result_store_exclusively_round_trips_and_rejects_tamper(self) -> None:
        source = manifest()
        record = self.collect(source)
        store = ResultStore(self.repository)
        run_path = store.create_run(source["run_id"], source)
        execution_path = store.write_execution_record(source["run_id"], record)
        self.assertEqual(execution_path, run_path / "execution.json")
        self.assertEqual(store.read_execution_record(source["run_id"]), record)
        with self.assertRaisesRegex(StorageError, "already exists"):
            store.write_execution_record(source["run_id"], record)

        tampered = deepcopy(record)
        tampered["recorded_at"] = "2026-10-07T00:00:00+00:00"
        execution_path.write_text(json.dumps(tampered), encoding="utf-8")
        with self.assertRaisesRegex(StorageError, "record_hash"):
            store.read_execution_record(source["run_id"])

    def test_comparison_identity_is_worktree_neutral_but_semantic(self) -> None:
        external = self.repository.parent / "external" / "libblas.a"
        external.parent.mkdir()
        external.write_bytes(b"external blas version one\n")
        original_stamp = self.adapter_stamp.read_text(encoding="utf-8")
        self.adapter_stamp.write_text(
            original_stamp.replace(
                f"USER_LIB={self.mumps_library}",
                f"USER_LIB={self.mumps_library} {external}",
            ),
            encoding="utf-8",
        )
        first_record = self.collect()
        first_identity = comparison_environment_identity(first_record)

        second_root = self.repository.parent / "repository-copy"
        shutil.copytree(self.repository, second_root)
        second_executable = (
            second_root
            / "benchmarking"
            / "build"
            / "release"
            / "EXEC"
            / "fake_manufactured"
        )
        second_executable.write_bytes(b"deliberately different executable bytes\n")
        second_stamp = (
            second_root
            / "benchmarking"
            / "build"
            / "release"
            / "benchmark_OBJ"
            / "fake"
            / ".ads-benchmark-adapter-build"
        )
        second_stamp.write_text(
            second_stamp.read_text(encoding="utf-8").replace(
                str(self.repository), str(second_root)
            ),
            encoding="utf-8",
        )
        second_environment = dict(self.environment)
        second_environment["MUMPS_DIR"] = str(
            second_root / "dependencies" / "mumps"
        )

        def second_record() -> dict[str, object]:
            return collect_execution_provenance(
                second_root,
                manifest(),
                launcher_descriptions=[self.launcher_description],
                environment=second_environment,
                recorded_at=self.instant,
            )

        equivalent = comparison_environment_identity(second_record())
        self.assertEqual(first_identity, equivalent)
        rendered = json.dumps(equivalent, sort_keys=True)
        self.assertNotIn(str(self.repository), rendered)
        self.assertNotIn(str(second_root), rendered)
        self.assertIn("<BENCHMARK_BUILD_ROOT>", rendered)
        self.assertIn(str(external), rendered)
        external_identity = next(
            library
            for library in equivalent["build"]["libraries"]
            if library["path"] == str(external)
        )
        self.assertRegex(external_identity["sha256"], r"^sha256:[0-9a-f]{64}$")

        unchanged_stamp = second_stamp.read_text(encoding="utf-8")
        second_stamp.write_text(
            unchanged_stamp.replace(
                str(self.profile),
                str(second_root / "benchmarking" / "build" / "release"),
            ).replace("FFLAGS=", "FFLAGS=-g ", 1),
            encoding="utf-8",
        )
        self.assertNotEqual(
            first_identity,
            comparison_environment_identity(second_record()),
        )
        second_stamp.write_text(unchanged_stamp, encoding="utf-8")

        second_stamp.write_text(
            unchanged_stamp.replace(f"FF={sys.executable}", "FF=/missing/compiler"),
            encoding="utf-8",
        )
        self.assertNotEqual(
            first_identity,
            comparison_environment_identity(second_record()),
        )
        second_stamp.write_text(unchanged_stamp, encoding="utf-8")

        external.write_bytes(b"external blas version two\n")
        self.assertNotEqual(
            first_identity,
            comparison_environment_identity(second_record()),
        )
        external.write_bytes(b"external blas version one\n")

        with patch(
            "ads_benchmark.framework.provenance.socket.gethostname",
            return_value="different-host",
        ):
            different_system = comparison_environment_identity(second_record())
        self.assertNotEqual(first_identity, different_system)

    def test_manifest_and_execution_record_can_be_created_as_one_owned_run(self) -> None:
        source = manifest("atomic-execution")
        record = self.collect(source)
        store = ResultStore(self.repository)
        run_path = store.create_run(
            source["run_id"], source, execution_record=record
        )
        self.assertEqual(store.read_manifest(source["run_id"]), source)
        self.assertEqual(store.read_execution_record(source["run_id"]), record)
        self.assertTrue((run_path / "manifest.json").is_file())
        self.assertTrue((run_path / "execution.json").is_file())


if __name__ == "__main__":
    unittest.main(verbosity=2)
