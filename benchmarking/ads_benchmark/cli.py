"""Command-line entry point for deterministic ADS benchmark runs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from .analysis import AnalysisPipeline, load_run_results, write_run_report
from .catalog import build_catalog
from .framework.config import load_profiles
from .framework.errors import BenchmarkError, ValidationError
from .framework.executor import ExecutionOutcome, Executor
from .framework.filtering import (
    CaseFilters,
    parse_degree_pair,
    parse_vector3,
    positive_integer_set,
)
from .framework.model import RepositoryState
from .framework.planner import (
    Plan,
    Planner,
    frozen_plan_from_manifest,
    validate_resume_request,
)
from .framework.provenance import inspect_repository
from .framework.registry import Catalog
from .framework.storage import ResultStore


DEFAULT_REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def _positive_integer(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be an integer") from error
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return parsed


def _add_selection_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--profile", default="smoke")
    parser.add_argument(
        "--repository-root", type=Path, default=DEFAULT_REPOSITORY_ROOT
    )
    parser.add_argument("--config-dir", type=Path)
    parser.add_argument("--problem", action="append", default=[])
    parser.add_argument("--scheme", action="append", default=[])
    parser.add_argument("--degree-pair", action="append", default=[])
    parser.add_argument("--mesh", action="append", default=[])
    parser.add_argument("--mpi-grid", action="append", default=[])
    parser.add_argument("--mpi-ranks", action="append", type=int, default=[])
    parser.add_argument("--omp", action="append", type=int, default=[])
    parser.add_argument("--steps", action="append", type=int, default=[])
    parser.add_argument(
        "--available-mpi-slots",
        type=_positive_integer,
        help="reject selected cases requiring more MPI ranks than this allocation",
    )
    parser.add_argument(
        "--available-cpu-slots",
        type=_positive_integer,
        help=(
            "reject selected cases whose MPI-rank times OpenMP-thread product "
            "exceeds this allocation"
        ),
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    plan = subparsers.add_parser("plan", help="expand and validate a profile")
    _add_selection_arguments(plan)
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
    run = subparsers.add_parser(
        "run", help="create a new run and execute every selected case"
    )
    _add_selection_arguments(run)
    run.add_argument(
        "--run-id",
        required=True,
        help="new directory name below benchmarks/ (existing runs are refused)",
    )
    run.add_argument(
        "--resume",
        action="store_true",
        help="resume this run ID instead of creating a new run",
    )
    resume = subparsers.add_parser(
        "resume", help="resume an existing compatible frozen run manifest"
    )
    _add_selection_arguments(resume)
    resume.add_argument(
        "--run-id",
        required=True,
        help="existing directory name below benchmarks/",
    )
    analyze = subparsers.add_parser(
        "analyze", help="analyze a complete frozen benchmark run"
    )
    analyze.add_argument(
        "--repository-root", type=Path, default=DEFAULT_REPOSITORY_ROOT
    )
    analyze.add_argument(
        "--run-id", required=True, help="existing directory name below benchmarks/"
    )
    analyze.add_argument(
        "--plot",
        action="store_true",
        help="also render a log-log PNG when matplotlib is available",
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
        steps=positive_integer_set(options.steps, "time steps"),
    )


def _planned_cases(
    options: argparse.Namespace,
) -> tuple[Path, Catalog, Plan, RepositoryState]:
    repository_root = options.repository_root.resolve()
    config_directory = (
        options.config_dir.resolve()
        if options.config_dir
        else repository_root / "benchmarking" / "configs"
    )
    catalog = build_catalog(
        available_mpi_slots=options.available_mpi_slots,
        available_cpu_slots=options.available_cpu_slots,
    )
    profiles = load_profiles(config_directory)
    planner = Planner(profiles, catalog)
    plan = planner.plan(options.profile, _filters(options, catalog))
    repository = inspect_repository(repository_root)
    return repository_root, catalog, plan, repository


def _require_validation_execution_resources(
    options: argparse.Namespace, plan: Plan
) -> None:
    """Require explicit scheduler/local capacity for real validation runs."""

    if not any(case.spec.family == "validation" for case in plan.cases):
        return
    missing = []
    if options.available_mpi_slots is None:
        missing.append("--available-mpi-slots")
    if options.available_cpu_slots is None:
        missing.append("--available-cpu-slots")
    if missing:
        raise ValidationError(
            "validation execution requires an explicit resource allocation: "
            + " and ".join(missing)
            + "; use plan/dry-run for structural validation without an allocation"
        )


def _plan(options: argparse.Namespace) -> int:
    repository_root, _, plan, repository = _planned_cases(options)

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


def _run(options: argparse.Namespace, *, resume: bool = False) -> int:
    repository_root, catalog, plan, repository = _planned_cases(options)
    _require_validation_execution_resources(options, plan)
    store = ResultStore(repository_root)
    executor = Executor(catalog, store, repository_root)
    if resume:
        run_path = store.results_root / options.run_id
    else:
        # Availability failures must not leave a run manifest or partial case
        # tree behind.  Planning/dry-run deliberately does not perform this
        # execution-only check.
        executor.preflight(plan.cases)
        manifest = plan.manifest(run_id=options.run_id, repository=repository)
        run_path = store.create_run(options.run_id, manifest)

    # Execution always uses a strict round-trip through the persisted manifest;
    # the in-memory planner result is used only as a compatibility assertion.
    frozen = frozen_plan_from_manifest(store.read_manifest(options.run_id), catalog)
    if frozen.run_id != options.run_id:
        raise ValidationError("run ID does not match its frozen manifest")
    validate_resume_request(frozen, plan, repository)
    if resume:
        executor.preflight(frozen.cases)

    case_count = len(frozen.cases)
    print(f"run:          {options.run_id}")
    print(f"profile:      {frozen.profile_name}")
    print(f"cases:        {case_count}")
    print(f"results:      {run_path}")
    print(f"mode:         {'resume' if resume else 'new'}")

    def report(position: int, total: int, outcome: ExecutionOutcome) -> None:
        case_id = outcome.case.case_id
        state = outcome.status.get("state")
        if outcome.skipped:
            print(f"[{position}/{total}] {case_id}: passed (verified, skipped)")
        elif state == "passed":
            print(f"[{position}/{total}] {case_id}: passed")
        else:
            detail = outcome.status.get("error")
            suffix = f" ({detail})" if detail else ""
            print(
                f"[{position}/{total}] {case_id}: {state}{suffix}",
                file=sys.stderr,
            )

    summary = executor.execute_frozen(
        options.run_id,
        frozen.cases,
        resume=resume,
        observer=report,
        preflight=False,
    )
    successful = summary.passed + summary.skipped
    print(
        f"summary:      passed={successful} failed={summary.failed} "
        f"total={case_count}"
    )
    if resume:
        print(
            f"resume:       skipped={summary.skipped} retried={case_count - summary.skipped}"
        )
    return 0 if summary.failed == 0 else 1


def _analyze(options: argparse.Namespace) -> int:
    repository_root = options.repository_root.resolve()
    catalog = build_catalog()
    store = ResultStore(repository_root)
    with store.execution_lock(options.run_id):
        frozen = frozen_plan_from_manifest(
            store.read_manifest(options.run_id), catalog
        )
        if frozen.run_id != options.run_id:
            raise ValidationError("run ID does not match its frozen manifest")
        families = {case.spec.family for case in frozen.cases}
        if len(families) != 1:
            raise ValidationError(
                "analysis requires exactly one experiment family per run"
            )
        family = catalog.families.get(next(iter(families)))
        if family.analyzer is None:
            raise ValidationError(
                f"experiment family {family.name} has no registered analyzer"
            )
        executor = Executor(catalog, store, repository_root)
        results = load_run_results(executor, options.run_id, frozen.cases)
        pipeline = AnalysisPipeline(catalog.analyzers)
        report = pipeline.analyze(
            family.analyzer,
            results,
            planned_cases=frozen.cases,
            source_run=options.run_id,
        )
        written = write_run_report(
            report, store, options.run_id, plot=options.plot
        )

    print(f"run:          {options.run_id}")
    print(f"analyzer:     {report.analyzer}")
    print(f"status:       {report.status}")
    print(
        f"series:       passed={report.summary['passed_series']} "
        f"failed={report.summary['failed_series']} "
        f"total={report.summary['series_count']}"
    )
    if report.family == "validation":
        print(
            f"field checks: parallel={report.summary['parallel_comparison_count']} "
            f"invalid-timings={report.summary['invalid_timing_count']} "
            f"analytic-failures={report.summary['analytic_final_failure_count']} "
            f"scheme-failures={report.summary['scheme_pair_failure_count']}"
        )
    print(f"json:         {written.json_path}")
    print(f"csv:          {written.csv_path}")
    if written.plot_path is not None:
        print(f"plot:         {written.plot_path}")
    elif written.plot_message is not None:
        print(f"plot:         {written.plot_message}")
    if report.diagnostics:
        print("diagnostics:")
        for diagnostic in report.diagnostics[:20]:
            print(f"  - {diagnostic}")
        if len(report.diagnostics) > 20:
            print(f"  - ... {len(report.diagnostics) - 20} more")
    return 0 if report.status == "passed" else 1


def main(arguments: list[str] | None = None) -> int:
    options = _parser().parse_args(arguments)
    try:
        if options.command == "plan":
            return _plan(options)
        if options.command == "run":
            return _run(options, resume=options.resume)
        if options.command == "resume":
            return _run(options, resume=True)
        if options.command == "analyze":
            return _analyze(options)
    except BenchmarkError as error:
        print(f"benchmark error: {error}", file=sys.stderr)
        return 2
    raise AssertionError(f"unhandled command: {options.command}")


if __name__ == "__main__":
    raise SystemExit(main())
