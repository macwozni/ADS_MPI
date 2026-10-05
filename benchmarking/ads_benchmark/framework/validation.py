"""Validation of fully expanded benchmark cases."""

from __future__ import annotations

from decimal import Decimal, InvalidOperation
from fractions import Fraction
import math

from .errors import BenchmarkError, RegistryError, ValidationError
from .model import (
    MAX_SAMPLE_POINTS_PER_AXIS,
    MAX_WRITTEN_SAMPLE_POINTS_PER_AXIS,
    CaseSpec,
    Vector3,
)
from .registry import Catalog


MAX_FORTRAN_INTEGER = 2_147_483_647
MAX_PROCESS_GRID_AXIS = 128
MAX_TIMEOUT_SECONDS = Decimal("2592000")


def decimal_value(text: str, field: str) -> Decimal:
    try:
        value = Decimal(text)
    except (InvalidOperation, TypeError, ValueError) as error:
        raise ValidationError(f"{field} is not a decimal number: {text!r}") from error
    if not value.is_finite():
        raise ValidationError(f"{field} must be finite")
    return value


def rational_value(text: str, field: str) -> Fraction:
    try:
        value = Fraction(text)
    except (ValueError, ZeroDivisionError, TypeError) as error:
        raise ValidationError(f"{field} is not an exact number: {text!r}") from error
    return value


def _positive_vector(value: Vector3, field: str) -> None:
    if len(value) != 3 or any(type(item) is not int or item <= 0 for item in value):
        raise ValidationError(f"{field} must contain three positive integers")


def _safe_product(values: Vector3, field: str) -> int:
    product = math.prod(values)
    if product > MAX_FORTRAN_INTEGER:
        raise ValidationError(f"{field} product exceeds a 32-bit Fortran integer")
    return product


def validate_case(case: CaseSpec, catalog: Catalog) -> None:
    """Reject every unsupported or internally inconsistent case."""

    lookups = (
        (catalog.families, case.family),
        (catalog.adapters, case.problem),
        (catalog.schemes, case.scheme),
        (catalog.exact_cases, case.exact_case),
        (catalog.build_profiles, case.build_profile),
        (catalog.launchers, case.launcher),
    )
    for registry, name in lookups:
        try:
            registry.get(name)
        except RegistryError as error:
            raise ValidationError(str(error)) from error

    final_time = rational_value(case.time.final_time, "final_time")
    time_step = rational_value(case.time.time_step, "time_step")
    if final_time <= 0:
        raise ValidationError("final_time must be positive")
    if time_step <= 0:
        raise ValidationError("time_step must be positive")
    if type(case.time.steps) is not int or case.time.steps <= 0:
        raise ValidationError("steps must be a positive integer")
    if time_step * case.time.steps != final_time:
        raise ValidationError("final_time/time_step must be an integral step count")

    _positive_vector(case.mesh, "mesh")
    _positive_vector(case.test_degree, "test_degree")
    _positive_vector(case.trial_degree, "trial_degree")
    _positive_vector(case.mpi.process_grid, "process_grid")
    _safe_product(case.mesh, "mesh")
    process_count = _safe_product(case.mpi.process_grid, "process_grid")

    for axis, (test_degree, trial_degree) in enumerate(
        zip(case.test_degree, case.trial_degree, strict=True), start=1
    ):
        if test_degree > 9 or trial_degree > 9:
            raise ValidationError(f"degree axis {axis} exceeds supported maximum 9")
        if test_degree <= trial_degree:
            raise ValidationError(
                f"test degree must exceed trial degree on axis {axis}"
            )

    if type(case.mpi.ranks) is not int or case.mpi.ranks <= 0:
        raise ValidationError("MPI ranks must be a positive integer")
    if case.mpi.ranks != process_count:
        raise ValidationError("MPI ranks must equal procx*procy*procz")
    for axis, processes in enumerate(case.mpi.process_grid, start=1):
        if processes > MAX_PROCESS_GRID_AXIS:
            raise ValidationError(
                f"process grid axis {axis} exceeds supported maximum "
                f"{MAX_PROCESS_GRID_AXIS}"
            )
    for axis, (processes, elements, trial_degree) in enumerate(
        zip(
            case.mpi.process_grid,
            case.mesh,
            case.trial_degree,
            strict=True,
        ),
        start=1,
    ):
        if processes > elements + trial_degree:
            raise ValidationError(
                f"process grid axis {axis} exceeds the {elements + trial_degree} "
                "trial-space DOFs available on that axis"
            )

    if type(case.openmp_threads) is not int or case.openmp_threads <= 0:
        raise ValidationError("OpenMP threads must be a positive integer")
    openmp_policy = (
        case.openmp_dynamic,
        case.openmp_proc_bind,
        case.openmp_places,
    )
    if any(value is not None for value in openmp_policy) and any(
        value is None for value in openmp_policy
    ):
        raise ValidationError(
            "OpenMP dynamic, proc_bind and places must be specified together"
        )
    if case.openmp_dynamic is not None and type(case.openmp_dynamic) is not bool:
        raise ValidationError("OpenMP dynamic must be a boolean")
    for value, field in (
        (case.openmp_proc_bind, "OpenMP proc_bind"),
        (case.openmp_places, "OpenMP places"),
    ):
        if value is not None and (
            not isinstance(value, str) or not value or "\0" in value
        ):
            raise ValidationError(f"{field} must be a nonempty NUL-free string")
    if (
        type(case.sampling.points_per_axis) is not int
        or case.sampling.points_per_axis < 2
        or case.sampling.points_per_axis > MAX_SAMPLE_POINTS_PER_AXIS
    ):
        raise ValidationError(
            "sampling points_per_axis must be between 2 and "
            f"{MAX_SAMPLE_POINTS_PER_AXIS}"
        )
    if type(case.sampling.write_samples) is not bool:
        raise ValidationError("sampling write_samples must be a boolean")
    if (
        case.sampling.write_samples
        and case.sampling.points_per_axis > MAX_WRITTEN_SAMPLE_POINTS_PER_AXIS
    ):
        raise ValidationError(
            "sampling points_per_axis must not exceed "
            f"{MAX_WRITTEN_SAMPLE_POINTS_PER_AXIS} when write_samples is true"
        )
    if type(case.measurement.warmups) is not int or case.measurement.warmups < 0:
        raise ValidationError("warmups must be a nonnegative integer")
    if type(case.measurement.samples) is not int or case.measurement.samples <= 0:
        raise ValidationError("samples must be a positive integer")
    timeout = decimal_value(case.measurement.timeout_seconds, "timeout_seconds")
    if timeout <= 0:
        raise ValidationError("timeout_seconds must be positive")
    if timeout > MAX_TIMEOUT_SECONDS:
        raise ValidationError(
            f"timeout_seconds must not exceed {MAX_TIMEOUT_SECONDS} seconds"
        )
    try:
        runtime_timeout = float(timeout)
    except (OverflowError, ValueError) as error:
        raise ValidationError("timeout_seconds is outside runtime range") from error
    if not math.isfinite(runtime_timeout) or runtime_timeout <= 0:
        raise ValidationError("timeout_seconds is outside runtime range")
    minimum_sample_seconds: Decimal | None = None
    if case.measurement.minimum_sample_seconds is not None:
        minimum_sample_seconds = decimal_value(
            case.measurement.minimum_sample_seconds,
            "minimum_sample_seconds",
        )
        if minimum_sample_seconds <= 0:
            raise ValidationError("minimum_sample_seconds must be positive")
        if minimum_sample_seconds > timeout:
            raise ValidationError(
                "minimum_sample_seconds must not exceed timeout_seconds"
            )

    if case.family == "strong":
        if case.build_profile != "release":
            raise ValidationError("strong scaling requires the release build profile")
        if case.exact_case != "spatial-cosine":
            raise ValidationError("strong scaling requires exact case spatial-cosine")
        if not case.sampling.write_samples:
            raise ValidationError("strong scaling requires write_samples=true")
        if case.openmp_dynamic is not False:
            raise ValidationError("strong scaling requires OMP_DYNAMIC=FALSE")
        if case.openmp_proc_bind is None or case.openmp_places is None:
            raise ValidationError("strong scaling requires a complete OpenMP policy")
        if minimum_sample_seconds is None:
            raise ValidationError(
                "strong scaling requires minimum_sample_seconds"
            )

    if case.family == "weak":
        weak = case.weak_scaling
        if weak is None:
            raise ValidationError("weak scaling requires weak_scaling metadata")
        _positive_vector(weak.local_elements, "weak local_elements")
        _safe_product(weak.local_elements, "weak local_elements")
        if weak.workload_basis != "per-rank":
            raise ValidationError(
                "weak scaling workload_basis must be per-rank"
            )
        if weak.role not in {"measurement", "field-reference"}:
            raise ValidationError(
                "weak scaling role must be measurement or field-reference"
            )
        if case.build_profile != "release":
            raise ValidationError("weak scaling requires the release build profile")
        if case.exact_case != "spatial-cosine":
            raise ValidationError("weak scaling requires exact case spatial-cosine")
        if not case.sampling.write_samples:
            raise ValidationError("weak scaling requires write_samples=true")
        if case.openmp_dynamic is not False:
            raise ValidationError("weak scaling requires OMP_DYNAMIC=FALSE")
        if case.openmp_proc_bind is None or case.openmp_places is None:
            raise ValidationError("weak scaling requires a complete OpenMP policy")
        if minimum_sample_seconds is None:
            raise ValidationError("weak scaling requires minimum_sample_seconds")

        if weak.role == "measurement":
            expected_global = tuple(
                local * processes
                for local, processes in zip(
                    weak.local_elements,
                    case.mpi.process_grid,
                    strict=True,
                )
            )
            if case.mesh != expected_global:
                raise ValidationError(
                    "weak measurement global mesh must equal "
                    "local_elements*process_grid"
                )
        else:
            if (
                case.mpi.ranks != 1
                or case.mpi.process_grid != (1, 1, 1)
                or case.openmp_threads != 1
            ):
                raise ValidationError(
                    "weak field-reference must use MPI=1, grid=1x1x1, OMP=1"
                )
            if case.mesh != weak.local_elements:
                raise ValidationError(
                    "weak field-reference local_elements must equal global mesh"
                )
    elif case.weak_scaling is not None:
        raise ValidationError(
            "weak_scaling metadata is only valid for the weak family"
        )

    adapter = catalog.adapters.get(case.problem)
    try:
        adapter.validate_case(case)
    except BenchmarkError:
        raise
    except Exception as error:
        raise ValidationError(
            f"problem adapter {case.problem} rejected the case: {error}"
        ) from error

    launcher = catalog.launchers.get(case.launcher)
    try:
        launcher.validate_case(case)
    except BenchmarkError:
        raise
    except Exception as error:
        raise ValidationError(
            f"launcher {case.launcher} rejected the case: {error}"
        ) from error
