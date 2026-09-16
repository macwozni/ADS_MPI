"""Command-line entry point for deterministic ADS benchmark planning."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from .catalog import build_catalog
from .framework.config import load_profiles
from .framework.errors import BenchmarkError, ValidationError
from .framework.filtering import (
    CaseFilters,
    parse_degree_pair,
    parse_vector3,
    positive_integer_set,
)
from .framework.planner import Planner
from .framework.provenance import inspect_repository
from .framework.registry import Catalog
from .framework.storage import ResultStore


DEFAULT_REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    plan = subparsers.add_parser("plan", help="expand and validate a profile")
    plan.add_argument("--profile", default="smoke")
    plan.add_argument("--repository-root", type=Path, default=DEFAULT_REPOSITORY_ROOT)
    plan.add_argument("--config-dir", type=Path)
    plan.add_argument("--problem", action="append", default=[])
    plan.add_argument("--scheme", action="append", default=[])
    plan.add_argument("--degree-pair", action="append", default=[])
    plan.add_argument("--mesh", action="append", default=[])
    plan.add_argument("--mpi-grid", action="append", default=[])
    plan.add_argument("--mpi-ranks", action="append", type=int, default=[])
    plan.add_argument("--omp", action="append", type=int, default=[])
    mode = plan.add_mutually_exclusive_group()
    mode.add_argument(
        "--dry-run",
        action="store_true",
        help="validate and print a summary without filesystem writes (default)",
    )
    mode.add_argument(
        "--write-manifest",
        action="store_true",
        help="create benchmarks/<run-id>/manifest.json exclusively",
    )
    plan.add_argument("--run-id")
    plan.add_argument(
        "--json",
        action="store_true",
        help="print the complete plan manifest to stdout",
    )
    return parser


def _filters(options: argparse.Namespace, catalog: Catalog) -> CaseFilters:
    problems = frozenset(options.problem)
    schemes = frozenset(value.lower() for value in options.scheme)
    for problem in problems:
        if not catalog.adapters.contains(problem):
            raise ValidationError(f"unknown problem filter: {problem}")
    for scheme in schemes:
        if not catalog.schemes.contains(scheme):
            raise ValidationError(f"unknown scheme filter: {scheme}")

    degree_pairs = frozenset(parse_degree_pair(value) for value in options.degree_pair)
    for test_degree, trial_degree in degree_pairs:
        for axis, (test_value, trial_value) in enumerate(
            zip(test_degree, trial_degree, strict=True), start=1
        ):
            if test_value > 9 or trial_value > 9:
                raise ValidationError(
                    f"degree-pair filter axis {axis} exceeds maximum 9"
                )
            if test_value <= trial_value:
                raise ValidationError(
                    f"degree-pair filter requires test > trial on axis {axis}"
                )

    return CaseFilters(
        problems=problems,
        schemes=schemes,
        degree_pairs=degree_pairs,
        meshes=frozenset(parse_vector3(value, "mesh") for value in options.mesh),
        mpi_grids=frozenset(
            parse_vector3(value, "MPI grid") for value in options.mpi_grid
        ),
        mpi_ranks=positive_integer_set(options.mpi_ranks, "MPI ranks"),
        openmp_threads=positive_integer_set(options.omp, "OpenMP"),
    )


def _plan(options: argparse.Namespace) -> int:
    repository_root = options.repository_root.resolve()
    config_directory = (
        options.config_dir.resolve()
        if options.config_dir
        else repository_root / "benchmarking" / "configs"
    )
    catalog = build_catalog()
    profiles = load_profiles(config_directory)
    planner = Planner(profiles, catalog)
    plan = planner.plan(options.profile, _filters(options, catalog))
    repository = inspect_repository(repository_root)

    if options.run_id and not options.write_manifest:
        raise ValidationError("--run-id is meaningful only with --write-manifest")
    if options.write_manifest and not options.run_id:
        raise ValidationError("--write-manifest requires --run-id")

    run_id = options.run_id or "dry-run"
    manifest = plan.manifest(run_id=run_id, repository=repository)
    written_path: Path | None = None
    if options.write_manifest:
        store = ResultStore(repository_root)
        written_path = store.create_run(options.run_id, manifest) / "manifest.json"

    if options.json:
        print(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False))
    else:
        first_case = plan.cases[0].case_id
        last_case = plan.cases[-1].case_id
        print(f"profile:      {plan.profile.name}")
        print(f"cases:        {len(plan.cases)}")
        print(f"config hash:  sha256:{plan.config_hash}")
        print(f"repository:   {repository.commit} (dirty={str(repository.dirty).lower()})")
        print(f"first case:   {first_case}")
        print(f"last case:    {last_case}")
        if written_path is None:
            print("mode:         dry-run (no files written)")
        else:
            print(f"manifest:     {written_path}")
    return 0


def main(arguments: list[str] | None = None) -> int:
    options = _parser().parse_args(arguments)
    try:
        if options.command == "plan":
            return _plan(options)
    except BenchmarkError as error:
        print(f"benchmark error: {error}", file=sys.stderr)
        return 2
    raise AssertionError(f"unhandled command: {options.command}")


if __name__ == "__main__":
    raise SystemExit(main())
