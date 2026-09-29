"""Stage-two oracle for the paired manufactured smoke runs."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
import sys
from typing import Any

from .framework.errors import BenchmarkError
from .framework.storage import validate_identifier


INITIAL_ERROR_LIMIT = 1.0e-10
PROBLEMS = ("igrm_l2", "igrm_heat", "pure_diffusion_igrm")
SCHEMES = ("dg", "pr", "be")
EXPECTED_CASES = frozenset(
    (problem, scheme) for problem in PROBLEMS for scheme in SCHEMES
)


class SmokeRefinementError(Exception):
    """A paired smoke result is incomplete or violates its numeric oracle."""


@dataclass(frozen=True)
class SmokeResult:
    path: Path
    problem: str
    scheme: str
    initial_l2_error: float
    initial_linf_error: float
    l2_error: float


def _load_json(path: Path) -> dict[str, Any]:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise SmokeRefinementError(f"cannot read result {path}: {error}") from error
    if not isinstance(document, dict):
        raise SmokeRefinementError(f"result is not a JSON object: {path}")
    return document


def _finite_nonnegative(document: dict[str, Any], field: str, path: Path) -> float:
    value = document.get(field)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SmokeRefinementError(f"{path}: {field} must be a finite number")
    normalized = float(value)
    if not math.isfinite(normalized) or normalized < 0.0:
        raise SmokeRefinementError(
            f"{path}: {field} must be finite and nonnegative"
        )
    return normalized


def _load_run(
    run_directory: Path, expected_steps: int
) -> dict[tuple[str, str], SmokeResult]:
    if not run_directory.is_dir():
        raise SmokeRefinementError(f"smoke run directory is missing: {run_directory}")
    paths = sorted(run_directory.glob("cases/*/result.json"))
    results: dict[tuple[str, str], SmokeResult] = {}
    for path in paths:
        document = _load_json(path)
        if document.get("status") != "passed":
            raise SmokeRefinementError(f"result is not marked passed: {path}")
        domain = document.get("domain_result")
        if not isinstance(domain, dict):
            raise SmokeRefinementError(f"result lacks domain_result: {path}")
        problem = domain.get("problem")
        scheme = domain.get("scheme")
        if not isinstance(problem, str) or not isinstance(scheme, str):
            raise SmokeRefinementError(f"result lacks problem/scheme strings: {path}")
        key = (problem, scheme)
        if key not in EXPECTED_CASES:
            raise SmokeRefinementError(
                f"unexpected smoke result {problem}/{scheme}: {path}"
            )
        if key in results:
            raise SmokeRefinementError(
                f"duplicate smoke result for {problem}/{scheme}: {path}"
            )
        steps = domain.get("steps")
        if type(steps) is not int or steps != expected_steps:
            raise SmokeRefinementError(
                f"{path}: expected {expected_steps} physical steps, got {steps!r}"
            )
        result = SmokeResult(
            path=path,
            problem=problem,
            scheme=scheme,
            initial_l2_error=_finite_nonnegative(
                domain, "initial_l2_error", path
            ),
            initial_linf_error=_finite_nonnegative(
                domain, "initial_linf_error", path
            ),
            l2_error=_finite_nonnegative(domain, "l2_error", path),
        )
        if (
            result.initial_l2_error > INITIAL_ERROR_LIMIT
            or result.initial_linf_error > INITIAL_ERROR_LIMIT
        ):
            raise SmokeRefinementError(
                f"{problem}/{scheme} initial projection exceeds "
                f"{INITIAL_ERROR_LIMIT:g} "
                f"in {path.parent.parent.parent.name}: "
                f"L2={result.initial_l2_error:.17g}, "
                f"Linf={result.initial_linf_error:.17g}"
            )
        results[key] = result

    missing = sorted(EXPECTED_CASES - set(results))
    if missing:
        labels = ", ".join(f"{problem}/{scheme}" for problem, scheme in missing)
        raise SmokeRefinementError(
            f"run {run_directory.name} is missing smoke results: {labels}"
        )
    return results


def verify_refinement(
    coarse_run: Path, refined_run: Path
) -> tuple[tuple[str, str, float, float], ...]:
    """Verify the exact 3x3 matrix and strict L2 improvement from N=4 to N=8."""

    coarse = _load_run(coarse_run, expected_steps=4)
    refined = _load_run(refined_run, expected_steps=8)
    rows: list[tuple[str, str, float, float]] = []
    for problem in PROBLEMS:
        for scheme in SCHEMES:
            key = (problem, scheme)
            coarse_error = coarse[key].l2_error
            refined_error = refined[key].l2_error
            if not refined_error < coarse_error:
                raise SmokeRefinementError(
                    f"{problem}/{scheme} did not improve under time refinement: "
                    f"N=4 L2={coarse_error:.17g}, N=8 L2={refined_error:.17g}"
                )
            rows.append((problem, scheme, coarse_error, refined_error))
    return tuple(rows)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--coarse-run-id", required=True)
    parser.add_argument("--refined-run-id", required=True)
    return parser


def main(arguments: list[str] | None = None) -> int:
    options = _parser().parse_args(arguments)
    try:
        coarse_id = validate_identifier(options.coarse_run_id, "coarse run_id")
        refined_id = validate_identifier(options.refined_run_id, "refined run_id")
        results_root = options.repository_root.resolve() / "benchmarks"
        rows = verify_refinement(
            results_root / coarse_id,
            results_root / refined_id,
        )
    except (BenchmarkError, SmokeRefinementError) as error:
        print(f"smoke refinement failed: {error}", file=sys.stderr)
        return 1

    for problem, scheme, coarse_error, refined_error in rows:
        print(
            f"{problem:22s} {scheme:2s} "
            f"N=4 L2={coarse_error:.9e} N=8 L2={refined_error:.9e}"
        )
    print("smoke refinement: passed 9/9")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
