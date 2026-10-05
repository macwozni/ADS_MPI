"""Command-line entry point for deterministic ADS benchmark runs."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from contextlib import ExitStack
import json
from pathlib import Path
import shlex
import sys

from .analysis import (
    AnalysisPipeline,
    AnalysisReport,
    WrittenAnalysis,
    load_run_results,
    write_run_report,
)
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
from .framework.model import ExecutionContext, PlannedCase, RepositoryState
from .framework.planner import (
    FrozenPlan,
    Plan,
    Planner,
    frozen_plan_from_manifest,
    validate_resume_request,
)
from .framework.provenance import inspect_repository
from .framework.registry import Catalog
from .framework.sharding import (
    ValidatedShardSet,
    merge_shard_manifests,
    shard_manifest,
    shard_subset_manifest,
    validate_shard_manifests,
)
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


def _nonnegative_integer(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be an integer") from error
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be nonnegative")
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
    parser.add_argument(
        "--launcher-template",
        help=(
            "shell-free scheduler argv template; supports {ranks}, {threads}, "
            "{procx}, {procy}, {procz}, and one final {payload}"
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
    plan.add_argument(
        "--show-commands",
        action="store_true",
        help="print the final shell-quoted launcher argv for every planned case",
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
        help="also render the family-specific PNG when matplotlib is available",
    )
    shard = subparsers.add_parser(
        "shard", help="create deterministic shard manifests for one frozen plan"
    )
    _add_selection_arguments(shard)
    shard.add_argument("--run-id", required=True, help="new parent plan run ID")
    shard.add_argument(
        "--strategy",
        choices=("index", "group", "problem-scheme-degree"),
        default="problem-scheme-degree",
    )
    shard.add_argument(
        "--shard-count",
        type=_positive_integer,
        help="required only for index sharding",
    )
    run_shard = subparsers.add_parser(
        "run-shard", help="execute one generated shard as an isolated run"
    )
    run_shard.add_argument(
        "--repository-root", type=Path, default=DEFAULT_REPOSITORY_ROOT
    )
    run_shard.add_argument("--parent-run-id", required=True)
    run_shard.add_argument("--shard-index", required=True, type=_nonnegative_integer)
    run_shard.add_argument("--run-id", required=True)
    run_shard.add_argument("--resume", action="store_true")
    run_shard.add_argument("--available-mpi-slots", type=_positive_integer)
    run_shard.add_argument("--available-cpu-slots", type=_positive_integer)
    run_shard.add_argument("--launcher-template")
    merge = subparsers.add_parser(
        "merge-shards",
        help="verify a complete shard-run set and analyze it as the parent plan",
    )
    merge.add_argument(
        "--repository-root", type=Path, default=DEFAULT_REPOSITORY_ROOT
    )
    merge.add_argument("--parent-run-id", required=True)
    merge.add_argument(
        "--shard-run",
        action="append",
        required=True,
        help="completed isolated shard run ID; repeat for every shard",
    )
    merge.add_argument("--plot", action="store_true")
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
    catalog = _catalog(options)
    profiles = load_profiles(config_directory)
    planner = Planner(profiles, catalog)
    plan = planner.plan(options.profile, _filters(options, catalog))
    repository = inspect_repository(repository_root)
    return repository_root, catalog, plan, repository


def _catalog(options: argparse.Namespace) -> Catalog:
    """Compose launchers from explicit CLI data, with user-facing errors."""

    try:
        return build_catalog(
            available_mpi_slots=getattr(options, "available_mpi_slots", None),
            available_cpu_slots=getattr(options, "available_cpu_slots", None),
            launcher_template=getattr(options, "launcher_template", None),
        )
    except ValueError as error:
        raise ValidationError(f"invalid launcher configuration: {error}") from error


def _require_validation_execution_resources(
    options: argparse.Namespace,
    plan_or_cases: Plan | Sequence[PlannedCase],
) -> None:
    """Require explicit scheduler/local capacity for parallel correctness runs."""

    cases = tuple(getattr(plan_or_cases, "cases", plan_or_cases))
    protected_families = {
        case.spec.family
        for case in cases
        if case.spec.family in {"validation", "strong", "weak"}
    }
    if not protected_families:
        return
    missing = []
    if options.available_mpi_slots is None:
        missing.append("--available-mpi-slots")
    if options.available_cpu_slots is None:
        missing.append("--available-cpu-slots")
    if missing:
        labels = {
            "validation": "validation",
            "strong": "strong-scaling",
            "weak": "weak-scaling",
        }
        family_label = (
            labels[next(iter(protected_families))]
            if len(protected_families) == 1
            else "parallel validation/scaling"
        )
        raise ValidationError(
            f"{family_label} execution requires an explicit resource allocation: "
            + " and ".join(missing)
            + "; use plan/dry-run for structural validation without an allocation"
        )


def _preview_command(
    repository_root: Path,
    catalog: Catalog,
    case: PlannedCase,
    run_id: str,
) -> str:
    """Render the exact shell-quoted argv without invoking a shell or process."""

    context = ExecutionContext(
        repository_root=repository_root,
        case_directory=(
            repository_root / "benchmarks" / run_id / "cases" / case.case_id
        ),
    )
    try:
        adapter = catalog.adapters.get(case.spec.problem)
        payload = tuple(adapter.build_payload_command(case.spec, context))
        launcher = catalog.launchers.get(case.spec.launcher)
        command = tuple(launcher.command(payload, case.spec))
    except BenchmarkError:
        raise
    except Exception as error:
        raise ValidationError(
            f"cannot construct command for {case.case_id}: {error}"
        ) from error
    if not command or any(
        not isinstance(argument, str) or not argument or "\0" in argument
        for argument in command
    ):
        raise ValidationError(
            f"launcher produced invalid argv for {case.case_id}"
        )
    return shlex.join(command)


def _plan(options: argparse.Namespace) -> int:
    repository_root, catalog, plan, repository = _planned_cases(options)

    if options.run_id and not options.write_manifest:
        raise ValidationError("--run-id is meaningful only with --write-manifest")
    if options.write_manifest and not options.run_id:
        raise ValidationError("--write-manifest requires --run-id")
    if options.json and options.show_commands:
        raise ValidationError("--json and --show-commands cannot be combined")

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
        if options.show_commands:
            print("commands:")
            for case in plan.cases:
                command = _preview_command(
                    repository_root, catalog, case, run_id
                )
                print(f"  {case.case_id}: {command}")
    return 0


def _execute_frozen_run(
    executor: Executor,
    frozen: FrozenPlan,
    *,
    run_path: Path,
    resume: bool,
    parent_run_id: str | None = None,
    shard_index: int | None = None,
) -> int:
    """Execute one already persisted and round-tripped immutable manifest."""

    case_count = len(frozen.cases)
    print(f"run:          {frozen.run_id}")
    print(f"profile:      {frozen.profile_name}")
    print(f"cases:        {case_count}")
    print(f"results:      {run_path}")
    print(f"mode:         {'resume' if resume else 'new'}")
    if parent_run_id is not None:
        print(f"parent:       {parent_run_id}")
    if shard_index is not None:
        print(f"shard:        {shard_index}")

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
        frozen.run_id,
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

    return _execute_frozen_run(
        executor,
        frozen,
        run_path=run_path,
        resume=resume,
    )


def _analyze_results(
    catalog: Catalog,
    frozen: FrozenPlan,
    results: Sequence[Mapping[str, object]],
    *,
    source_run: str,
) -> AnalysisReport:
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
    pipeline = AnalysisPipeline(catalog.analyzers)
    return pipeline.analyze(
        family.analyzer,
        results,
        planned_cases=frozen.cases,
        source_run=source_run,
    )


def _print_analysis(
    run_id: str,
    report: AnalysisReport,
    written: WrittenAnalysis,
) -> int:
    print(f"run:          {run_id}")
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
    elif report.family in {"strong", "weak"}:
        print(
            f"scaling checks: invalid-timings={report.summary['invalid_timing_count']} "
            f"unreliable={report.summary['unreliable_measurement_count']} "
            f"field-failures={report.summary['field_failure_count']}"
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


def _analyze(options: argparse.Namespace) -> int:
    repository_root = options.repository_root.resolve()
    catalog = _catalog(options)
    store = ResultStore(repository_root)
    with store.execution_lock(options.run_id):
        frozen = frozen_plan_from_manifest(
            store.read_manifest(options.run_id), catalog
        )
        if frozen.run_id != options.run_id:
            raise ValidationError("run ID does not match its frozen manifest")
        executor = Executor(catalog, store, repository_root)
        results = load_run_results(executor, options.run_id, frozen.cases)
        report = _analyze_results(
            catalog, frozen, results, source_run=options.run_id
        )
        written = write_run_report(
            report, store, options.run_id, plot=options.plot
        )
    return _print_analysis(options.run_id, report, written)


def _shard(options: argparse.Namespace) -> int:
    """Freeze one complete plan, then write a self-validating shard set."""

    repository_root, _, plan, repository = _planned_cases(options)
    manifest = plan.manifest(run_id=options.run_id, repository=repository)
    shards = shard_manifest(
        manifest,
        strategy=options.strategy,
        shard_count=options.shard_count,
    )
    reconstructed = merge_shard_manifests(shards)
    if reconstructed != manifest:
        raise ValidationError("generated shards do not reconstruct the parent plan")

    store = ResultStore(repository_root)
    run_path = store.create_run(options.run_id, manifest)
    store.write_shard_manifests(options.run_id, shards)

    persisted = store.read_all_shard_manifests(options.run_id)
    if merge_shard_manifests(persisted) != manifest:
        raise ValidationError("persisted shards do not reconstruct the parent plan")
    print(f"parent:       {options.run_id}")
    print(f"profile:      {plan.profile.name}")
    print(f"strategy:     {shards[0]['strategy']}")
    print(f"shards:       {len(shards)}")
    print(f"cases:        {len(plan.cases)}")
    print(f"manifests:    {run_path / 'shards'}")
    return 0


def _require_repository_identity(
    expected: RepositoryState,
    current: RepositoryState,
) -> None:
    if expected.commit != current.commit:
        raise ValidationError(
            "shard repository SHA mismatch: "
            f"manifest={expected.commit} current={current.commit}"
        )
    if expected.dirty != current.dirty:
        raise ValidationError(
            "shard repository dirty-state mismatch with parent manifest"
        )
    if expected.worktree_fingerprint != current.worktree_fingerprint:
        raise ValidationError(
            "shard repository worktree fingerprint mismatch with parent manifest"
        )


def _validated_parent_shards(
    store: ResultStore,
    parent_run_id: str,
) -> tuple[tuple[dict[str, object], ...], ValidatedShardSet]:
    documents = store.read_all_shard_manifests(parent_run_id)
    validated = validate_shard_manifests(documents)
    parent_manifest = store.read_manifest(parent_run_id)
    if validated.parent_manifest != parent_manifest:
        raise ValidationError(
            "shard manifests do not match their persisted parent plan"
        )
    return documents, validated


def _run_shard(options: argparse.Namespace) -> int:
    repository_root = options.repository_root.resolve()
    if options.run_id == options.parent_run_id:
        raise ValidationError("shard run ID must differ from its parent run ID")
    catalog = _catalog(options)
    store = ResultStore(repository_root)
    documents, validated = _validated_parent_shards(
        store, options.parent_run_id
    )
    shard_count = validated.shard_count
    if options.shard_index >= shard_count:
        raise ValidationError(
            f"shard index {options.shard_index} is outside 0..{shard_count - 1}"
        )
    shard = next(
        document
        for document in documents
        if document.get("shard_index") == options.shard_index
    )
    expected_manifest = shard_subset_manifest(shard, run_id=options.run_id)
    frozen = frozen_plan_from_manifest(expected_manifest, catalog)
    current_repository = inspect_repository(repository_root)
    _require_repository_identity(frozen.repository, current_repository)
    _require_validation_execution_resources(options, frozen.cases)

    executor = Executor(catalog, store, repository_root)
    if options.resume:
        persisted = store.read_manifest(options.run_id)
        if persisted != expected_manifest:
            raise ValidationError(
                "shard resume manifest does not exactly match the generated subset"
            )
        frozen = frozen_plan_from_manifest(persisted, catalog)
        run_path = store.results_root / options.run_id
    else:
        # Resource/binary failures must not leave an isolated child run behind.
        executor.preflight(frozen.cases)
        run_path = store.create_run(options.run_id, expected_manifest)
        frozen = frozen_plan_from_manifest(
            store.read_manifest(options.run_id), catalog
        )
    if frozen.run_id != options.run_id:
        raise ValidationError("shard run ID does not match its frozen manifest")
    if options.resume:
        executor.preflight(frozen.cases)
    return _execute_frozen_run(
        executor,
        frozen,
        run_path=run_path,
        resume=options.resume,
        parent_run_id=options.parent_run_id,
        shard_index=options.shard_index,
    )


def _merge_shards(options: argparse.Namespace) -> int:
    """Verify every isolated child, then analyze exactly the parent case set."""

    repository_root = options.repository_root.resolve()
    catalog = _catalog(options)
    store = ResultStore(repository_root)
    documents, validated = _validated_parent_shards(
        store, options.parent_run_id
    )
    parent_frozen = frozen_plan_from_manifest(
        validated.parent_manifest, catalog
    )
    if parent_frozen.run_id != options.parent_run_id:
        raise ValidationError("parent run ID does not match its frozen manifest")

    shard_by_index = {
        int(document["shard_index"]): document for document in documents
    }
    expected_by_case_ids: dict[frozenset[str], int] = {
        frozenset(case_ids): index
        for index, case_ids in enumerate(validated.case_ids_by_shard)
    }
    if len(set(options.shard_run)) != len(options.shard_run):
        raise ValidationError("duplicate --shard-run identifiers are not allowed")
    child_by_index: dict[int, tuple[str, FrozenPlan]] = {}
    for child_run_id in options.shard_run:
        if child_run_id == options.parent_run_id:
            raise ValidationError("parent plan cannot also be a shard result run")
        child_manifest = store.read_manifest(child_run_id)
        child_frozen = frozen_plan_from_manifest(child_manifest, catalog)
        case_ids = frozenset(case.case_id for case in child_frozen.cases)
        shard_index = expected_by_case_ids.get(case_ids)
        if shard_index is None:
            raise ValidationError(
                f"run {child_run_id} does not match any expected shard case set"
            )
        if shard_index in child_by_index:
            previous = child_by_index[shard_index][0]
            raise ValidationError(
                f"runs {previous} and {child_run_id} both claim shard {shard_index}"
            )
        expected = shard_subset_manifest(
            shard_by_index[shard_index], run_id=child_run_id
        )
        if child_manifest != expected:
            raise ValidationError(
                f"run {child_run_id} conflicts with shard {shard_index} manifest"
            )
        child_by_index[shard_index] = (child_run_id, child_frozen)

    missing = sorted(set(range(validated.shard_count)) - set(child_by_index))
    if missing:
        raise ValidationError(
            "missing completed shard runs for indices: "
            + ", ".join(str(index) for index in missing)
        )

    executor = Executor(catalog, store, repository_root)
    results: list[Mapping[str, object]] = []
    locked_runs = sorted(
        {options.parent_run_id, *(run_id for run_id, _ in child_by_index.values())}
    )
    with ExitStack() as stack:
        for run_id in locked_runs:
            stack.enter_context(store.execution_lock(run_id))
        for index in range(validated.shard_count):
            child_run_id, child_frozen = child_by_index[index]
            results.extend(
                load_run_results(executor, child_run_id, child_frozen.cases)
            )
        if len(results) != len(parent_frozen.cases):
            raise ValidationError(
                "merged result count does not match the frozen parent plan"
            )
        report = _analyze_results(
            catalog,
            parent_frozen,
            results,
            source_run=options.parent_run_id,
        )
        written = write_run_report(
            report,
            store,
            options.parent_run_id,
            plot=options.plot,
        )

    print(f"shards:       {validated.shard_count}")
    print(f"child runs:   {len(child_by_index)}")
    return _print_analysis(options.parent_run_id, report, written)


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
        if options.command == "shard":
            return _shard(options)
        if options.command == "run-shard":
            return _run_shard(options)
        if options.command == "merge-shards":
            return _merge_shards(options)
    except BenchmarkError as error:
        print(f"benchmark error: {error}", file=sys.stderr)
        return 2
    raise AssertionError(f"unhandled command: {options.command}")


if __name__ == "__main__":
    raise SystemExit(main())
