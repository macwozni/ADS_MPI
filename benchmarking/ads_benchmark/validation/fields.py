"""Strict, dependency-free comparison of sampled three-dimensional fields.

The on-disk interchange format is CSV with coordinates in the first three
columns and one or more named scalar components after them.  Rows use the
canonical order ``x`` fastest, then ``y``, then ``z``.  Parsing normalizes the
text into a compact regular-grid representation, so comparison is independent
of CSV spelling and of the MPI decomposition that produced the samples.
"""

from __future__ import annotations

from array import array
from collections.abc import Mapping, Sequence
import csv
from dataclasses import dataclass
import hashlib
import io
import math
import struct
import sys
from typing import Protocol


Coordinate3 = tuple[float, float, float]
GridIndex3 = tuple[int, int, int]
GridShape3 = tuple[int, int, int]

DEFAULT_MAXIMUM_POINTS = 257**3
DEFAULT_MAXIMUM_COMPONENTS = 64
_COORDINATE_ABSOLUTE_TOLERANCE = 2.0e-14
_COORDINATE_RELATIVE_TOLERANCE = 64.0 * sys.float_info.epsilon


class FieldComparisonError(ValueError):
    """A sampled field is malformed or two fields are not comparable."""


class _ReadableText(Protocol):
    def read_text(self) -> str:
        """Return the already safety-checked artifact contents."""


class _Digest(Protocol):
    def update(self, value: bytes) -> None:
        """Add bytes to a stable checksum."""


@dataclass(frozen=True)
class RegularGridField:
    """A compact, canonical regular-grid field returned by the CSV parser."""

    label: str
    coordinate_names: tuple[str, str, str]
    shape: GridShape3
    axes: tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...]]
    component_names: tuple[str, ...]
    _component_values: tuple[array, ...]

    @property
    def point_count(self) -> int:
        return self.shape[0] * self.shape[1] * self.shape[2]

    def component_values(self, name: str) -> memoryview:
        """Expose one component through a read-only, zero-copy view."""

        try:
            index = self.component_names.index(name)
        except ValueError as error:
            raise FieldComparisonError(
                f"{self.label}: unknown field component {name!r}"
            ) from error
        return memoryview(self._component_values[index]).toreadonly()

    def grid_index(self, flat_index: int) -> GridIndex3:
        if type(flat_index) is not int or not 0 <= flat_index < self.point_count:
            raise FieldComparisonError(
                f"{self.label}: sample index {flat_index!r} is out of range"
            )
        nx, ny, _ = self.shape
        iz, remainder = divmod(flat_index, nx * ny)
        iy, ix = divmod(remainder, nx)
        return ix, iy, iz

    def coordinates(self, flat_index: int) -> Coordinate3:
        ix, iy, iz = self.grid_index(flat_index)
        return self.axes[0][ix], self.axes[1][iy], self.axes[2][iz]


@dataclass(frozen=True)
class ComponentStatistics:
    name: str
    l2_norm: float
    linf_norm: float

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "l2_norm": self.l2_norm,
            "linf_norm": self.linf_norm,
        }


@dataclass(frozen=True)
class FieldStatistics:
    l2_norm: float
    linf_norm: float
    checksum: str
    components: tuple[ComponentStatistics, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "l2_norm": self.l2_norm,
            "linf_norm": self.linf_norm,
            "checksum": self.checksum,
            "components": [component.to_dict() for component in self.components],
        }


@dataclass(frozen=True)
class PointDifference:
    flat_index: int
    grid_index: GridIndex3
    coordinates: Coordinate3
    component: str
    reference: float
    candidate: float
    absolute_difference: float
    relative_difference: float
    allowed_difference: float

    def to_dict(self) -> dict[str, object]:
        return {
            "flat_index": self.flat_index,
            "grid_index": list(self.grid_index),
            "coordinates": list(self.coordinates),
            "component": self.component,
            "reference": self.reference,
            "candidate": self.candidate,
            "absolute_difference": self.absolute_difference,
            "relative_difference": self.relative_difference,
            "allowed_difference": self.allowed_difference,
        }


@dataclass(frozen=True)
class ComponentDifference:
    name: str
    l2_difference: float
    linf_difference: float
    mismatch_count: int

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "l2_difference": self.l2_difference,
            "linf_difference": self.linf_difference,
            "mismatch_count": self.mismatch_count,
        }


@dataclass(frozen=True)
class FieldComparison:
    passed: bool
    absolute_tolerance: float
    relative_tolerance: float
    point_count: int
    component_count: int
    mismatch_count: int
    l2_difference: float
    linf_difference: float
    reference: FieldStatistics
    candidate: FieldStatistics
    component_differences: tuple[ComponentDifference, ...]
    largest_difference_point: PointDifference
    worst_point: PointDifference

    @property
    def comparison_count(self) -> int:
        return self.point_count * self.component_count

    def diagnostic(self) -> str:
        """Return a concise human-readable explanation of the worst point."""

        worst = self.worst_point
        status = "passed" if self.passed else "FAILED"
        coordinate = ", ".join(f"{value:.17g}" for value in worst.coordinates)
        return (
            f"{status}: {self.mismatch_count}/{self.comparison_count} field values "
            f"outside tolerance; worst component={worst.component!r}, "
            f"coordinate=({coordinate}), grid_index={worst.grid_index}, "
            f"reference={worst.reference:.17g}, candidate={worst.candidate:.17g}, "
            f"abs_difference={worst.absolute_difference:.17g}, "
            f"allowed={worst.allowed_difference:.17g}, "
            f"L2(difference)={self.l2_difference:.17g}, "
            f"Linf(difference)={self.linf_difference:.17g}"
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "passed": self.passed,
            "absolute_tolerance": self.absolute_tolerance,
            "relative_tolerance": self.relative_tolerance,
            "point_count": self.point_count,
            "component_count": self.component_count,
            "comparison_count": self.comparison_count,
            "mismatch_count": self.mismatch_count,
            "l2_difference": self.l2_difference,
            "linf_difference": self.linf_difference,
            "reference": self.reference.to_dict(),
            "candidate": self.candidate.to_dict(),
            "component_differences": [
                component.to_dict() for component in self.component_differences
            ],
            "largest_difference_point": self.largest_difference_point.to_dict(),
            "worst_point": self.worst_point.to_dict(),
            "diagnostic": self.diagnostic(),
        }


class _ScaledSumSquares:
    """Streaming L2 accumulator with the overflow safety of LAPACK xLASSQ."""

    def __init__(self) -> None:
        self._scale = 0.0
        self._sum_squares = 1.0

    def add(self, value: float) -> None:
        magnitude = abs(value)
        if not math.isfinite(magnitude):
            raise FieldComparisonError("a derived weighted value is not finite")
        if magnitude == 0.0:
            return
        if self._scale < magnitude:
            ratio = self._scale / magnitude
            self._sum_squares = 1.0 + self._sum_squares * ratio * ratio
            self._scale = magnitude
        else:
            ratio = magnitude / self._scale
            self._sum_squares += ratio * ratio

    def result(self) -> float:
        if self._scale == 0.0:
            return 0.0
        result = self._scale * math.sqrt(self._sum_squares)
        if not math.isfinite(result):
            raise FieldComparisonError("a derived L2 value is not finite")
        return result


def _strict_positive_integer(value: object, name: str) -> int:
    if type(value) is not int or value <= 0:
        raise FieldComparisonError(f"{name} must be a positive integer")
    return value


def _normalize_shape(shape: int | Sequence[int]) -> GridShape3:
    if type(shape) is int:
        dimensions = (shape, shape, shape)
    elif isinstance(shape, Sequence) and not isinstance(shape, (str, bytes)):
        dimensions = tuple(shape)
    else:
        raise FieldComparisonError("shape must be an integer or three integers")
    if len(dimensions) != 3:
        raise FieldComparisonError("shape must contain exactly three dimensions")
    normalized = tuple(
        _strict_positive_integer(value, f"shape[{index}]")
        for index, value in enumerate(dimensions)
    )
    if any(value < 2 for value in normalized):
        raise FieldComparisonError("each regular-grid dimension must be at least 2")
    return normalized  # type: ignore[return-value]


def _finite_float(value: str, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise FieldComparisonError(f"{name} must be a finite number") from error
    if not math.isfinite(result):
        raise FieldComparisonError(f"{name} must be finite")
    return result


def _coordinate_matches(first: float, second: float) -> bool:
    return math.isclose(
        first,
        second,
        rel_tol=_COORDINATE_RELATIVE_TOLERANCE,
        abs_tol=_COORDINATE_ABSOLUTE_TOLERANCE,
    )


def _validate_regular_axis(axis: tuple[float, ...], name: str, label: str) -> None:
    if any(right <= left for left, right in zip(axis, axis[1:])):
        raise FieldComparisonError(
            f"{label}: coordinate axis {name!r} must be strictly increasing"
        )
    spacing = (axis[-1] - axis[0]) / (len(axis) - 1)
    if not math.isfinite(spacing) or spacing <= 0.0:
        raise FieldComparisonError(
            f"{label}: coordinate axis {name!r} has invalid spacing"
        )
    for index, actual in enumerate(axis):
        expected = axis[0] + index * spacing
        if not _coordinate_matches(actual, expected):
            raise FieldComparisonError(
                f"{label}: coordinate axis {name!r} is not regular at index "
                f"{index}: expected {expected:.17g}, found {actual:.17g}"
            )


def _artifact_text(source: str | _ReadableText, label: str) -> str:
    if isinstance(source, str):
        return source
    read_text = getattr(source, "read_text", None)
    if not callable(read_text):
        raise FieldComparisonError(
            f"{label}: field source must be CSV text or expose read_text()"
        )
    text = read_text()
    if not isinstance(text, str):
        raise FieldComparisonError(f"{label}: read_text() did not return text")
    return text


def parse_regular_grid_csv(
    source: str | _ReadableText,
    *,
    shape: int | Sequence[int],
    components: Sequence[str] | None = None,
    label: str = "field",
    maximum_points: int = DEFAULT_MAXIMUM_POINTS,
    maximum_components: int = DEFAULT_MAXIMUM_COMPONENTS,
) -> RegularGridField:
    """Parse strict ``x,y,z,<components...>`` samples in canonical grid order.

    Every numeric column is checked for finiteness, including columns not
    selected into the returned field.  ``source`` may be inline text or the
    lazy artifact reference used by the benchmark result store.
    """

    if not isinstance(label, str) or not label:
        raise FieldComparisonError("label must be a nonempty string")
    grid_shape = _normalize_shape(shape)
    maximum_points = _strict_positive_integer(maximum_points, "maximum_points")
    maximum_components = _strict_positive_integer(
        maximum_components, "maximum_components"
    )
    expected_count = grid_shape[0] * grid_shape[1] * grid_shape[2]
    if expected_count > maximum_points:
        raise FieldComparisonError(
            f"{label}: grid contains {expected_count} points, maximum is "
            f"{maximum_points}"
        )

    text = _artifact_text(source, label)
    reader = csv.reader(io.StringIO(text, newline=""), strict=True)
    try:
        header = next(reader)
    except StopIteration as error:
        raise FieldComparisonError(f"{label}: field CSV has no header") from error
    except csv.Error as error:
        raise FieldComparisonError(f"{label}: malformed field CSV: {error}") from error

    if len(header) < 4 or header[:3] != ["x", "y", "z"]:
        raise FieldComparisonError(
            f"{label}: field CSV header must start with x,y,z and include a "
            "component"
        )
    if any(not name for name in header) or len(set(header)) != len(header):
        raise FieldComparisonError(
            f"{label}: field CSV column names must be nonempty and unique"
        )
    available_components = tuple(header[3:])
    if len(available_components) > maximum_components:
        raise FieldComparisonError(
            f"{label}: field CSV has {len(available_components)} components, "
            f"maximum is {maximum_components}"
        )
    if components is None:
        selected_components = available_components
    else:
        if not isinstance(components, Sequence) or isinstance(
            components, (str, bytes)
        ):
            raise FieldComparisonError("components must be a sequence of names")
        selected_components = tuple(components)
        if not selected_components:
            raise FieldComparisonError("at least one component must be selected")
        if any(not isinstance(name, str) or not name for name in selected_components):
            raise FieldComparisonError(
                "selected component names must be nonempty strings"
            )
        if len(set(selected_components)) != len(selected_components):
            raise FieldComparisonError("selected component names must be unique")
        missing = [
            name for name in selected_components if name not in available_components
        ]
        if missing:
            raise FieldComparisonError(
                f"{label}: field CSV does not contain component(s): "
                + ", ".join(repr(name) for name in missing)
            )
    if len(selected_components) > maximum_components:
        raise FieldComparisonError(
            f"{label}: selected {len(selected_components)} components, maximum is "
            f"{maximum_components}"
        )

    selected_columns = tuple(header.index(name) for name in selected_components)
    values = tuple(array("d") for _ in selected_components)
    x_axis: list[float] = []
    y_axis: list[float] = []
    z_axis: list[float] = []
    nx, ny, _ = grid_shape
    row_count = 0
    try:
        for row_count, row in enumerate(reader, start=1):
            if row_count > expected_count:
                raise FieldComparisonError(
                    f"{label}: field CSV has more than {expected_count} data rows"
                )
            if len(row) != len(header):
                raise FieldComparisonError(
                    f"{label}: field CSV row {row_count + 1} has {len(row)} "
                    f"columns, expected {len(header)}"
                )
            parsed = tuple(
                _finite_float(value, f"{label}.row[{row_count}].{header[column]}")
                for column, value in enumerate(row)
            )
            flat_index = row_count - 1
            iz, remainder = divmod(flat_index, nx * ny)
            iy, ix = divmod(remainder, nx)
            x, y, z = parsed[:3]
            if iy == 0 and iz == 0:
                x_axis.append(x)
            if ix == 0 and iz == 0:
                y_axis.append(y)
            if ix == 0 and iy == 0:
                z_axis.append(z)
            expected_coordinates = (x_axis[ix], y_axis[iy], z_axis[iz])
            if any(
                not _coordinate_matches(actual, expected)
                for actual, expected in zip(
                    (x, y, z), expected_coordinates, strict=True
                )
            ):
                raise FieldComparisonError(
                    f"{label}: field CSV row {row_count + 1} is not in canonical "
                    "regular-grid order (x fastest, then y, then z)"
                )
            for target, column in zip(values, selected_columns, strict=True):
                target.append(parsed[column])
    except csv.Error as error:
        raise FieldComparisonError(f"{label}: malformed field CSV: {error}") from error

    if row_count != expected_count:
        raise FieldComparisonError(
            f"{label}: field CSV has {row_count} data rows, expected {expected_count}"
        )
    axes = (tuple(x_axis), tuple(y_axis), tuple(z_axis))
    for name, axis in zip(("x", "y", "z"), axes, strict=True):
        _validate_regular_axis(axis, name, label)
    return RegularGridField(
        label=label,
        coordinate_names=("x", "y", "z"),
        shape=grid_shape,
        axes=axes,
        component_names=selected_components,
        _component_values=values,
    )


def select_components(
    field: RegularGridField,
    names: Sequence[str] | Mapping[str, str],
    *,
    label: str | None = None,
) -> RegularGridField:
    """Return a zero-copy component view, optionally aliasing component names.

    A sequence preserves each selected source name.  A mapping associates
    ``source_name -> result_name``; this is useful for comparing differently
    named scalar columns such as ``numerical`` and ``exact`` under one shared
    alias.  Grid axes and component arrays are reused, never copied.
    """

    if not isinstance(field, RegularGridField):
        raise FieldComparisonError("field must be a parsed RegularGridField")
    selected_label = field.label if label is None else label
    if not isinstance(selected_label, str) or not selected_label:
        raise FieldComparisonError("label must be a nonempty string")
    if isinstance(names, Mapping):
        selections = tuple(names.items())
    elif isinstance(names, Sequence) and not isinstance(names, (str, bytes)):
        selections = tuple((name, name) for name in names)
    else:
        raise FieldComparisonError(
            "component selection must be a sequence or alias mapping"
        )
    if not selections:
        raise FieldComparisonError("at least one component must be selected")
    if any(
        not isinstance(source_name, str)
        or not source_name
        or not isinstance(result_name, str)
        or not result_name
        for source_name, result_name in selections
    ):
        raise FieldComparisonError(
            "selected source and result component names must be nonempty strings"
        )
    source_names = tuple(source_name for source_name, _ in selections)
    result_names = tuple(result_name for _, result_name in selections)
    if len(set(source_names)) != len(source_names):
        raise FieldComparisonError("selected source component names must be unique")
    if len(set(result_names)) != len(result_names):
        raise FieldComparisonError("selected result component names must be unique")
    missing = [name for name in source_names if name not in field.component_names]
    if missing:
        raise FieldComparisonError(
            f"{field.label}: unknown field component(s): "
            + ", ".join(repr(name) for name in missing)
        )
    values = tuple(
        field._component_values[field.component_names.index(source_name)]
        for source_name in source_names
    )
    return RegularGridField(
        label=selected_label,
        coordinate_names=field.coordinate_names,
        shape=field.shape,
        axes=field.axes,
        component_names=result_names,
        _component_values=values,
    )


def scalar_component(
    field: RegularGridField,
    source_name: str,
    *,
    name: str = "value",
    label: str | None = None,
) -> RegularGridField:
    """Return one zero-copy scalar component under a comparison-safe alias."""

    return select_components(field, {source_name: name}, label=label)


def _point_sqrt_weights(field: RegularGridField):
    nx, ny, nz = field.shape
    hx = (field.axes[0][-1] - field.axes[0][0]) / (nx - 1)
    hy = (field.axes[1][-1] - field.axes[1][0]) / (ny - 1)
    hz = (field.axes[2][-1] - field.axes[2][0]) / (nz - 1)
    cell_volume = hx * hy * hz
    if not math.isfinite(cell_volume) or cell_volume <= 0.0:
        raise FieldComparisonError(
            f"{field.label}: regular-grid cell volume is not finite and positive"
        )
    for iz in range(nz):
        z_weight = 0.5 if iz in (0, nz - 1) else 1.0
        for iy in range(ny):
            y_weight = 0.5 if iy in (0, ny - 1) else 1.0
            for ix in range(nx):
                x_weight = 0.5 if ix in (0, nx - 1) else 1.0
                yield math.sqrt(cell_volume * x_weight * y_weight * z_weight)


def _update_text(digest: _Digest, value: str) -> None:
    encoded = value.encode("utf-8", errors="strict")
    digest.update(struct.pack(">Q", len(encoded)))
    digest.update(encoded)


def _update_float(digest: _Digest, value: float) -> None:
    digest.update(struct.pack(">d", 0.0 if value == 0.0 else value))


def _field_checksum(field: RegularGridField) -> str:
    digest = hashlib.sha256()
    digest.update(b"ads-regular-grid-field-v1\0")
    digest.update(struct.pack(">QQQ", *field.shape))
    for name in field.coordinate_names:
        _update_text(digest, name)
    for axis in field.axes:
        for coordinate in axis:
            _update_float(digest, coordinate)
    components_by_name = {
        name: component
        for name, component in zip(
            field.component_names, field._component_values, strict=True
        )
    }
    for name in sorted(components_by_name):
        component = components_by_name[name]
        _update_text(digest, name)
        for value in component:
            _update_float(digest, value)
    return "sha256:" + digest.hexdigest()


def summarize_field(field: RegularGridField) -> FieldStatistics:
    """Calculate physical trapezoidal norms and a canonical SHA-256 checksum."""

    if not isinstance(field, RegularGridField):
        raise FieldComparisonError("field must be a parsed RegularGridField")
    aggregate = _ScaledSumSquares()
    component_accumulators = tuple(
        _ScaledSumSquares() for _ in field.component_names
    )
    component_linf = [0.0 for _ in field.component_names]
    for flat_index, sqrt_weight in enumerate(_point_sqrt_weights(field)):
        for component_index, component in enumerate(field._component_values):
            value = component[flat_index]
            weighted = value * sqrt_weight
            component_accumulators[component_index].add(weighted)
            aggregate.add(weighted)
            component_linf[component_index] = max(
                component_linf[component_index], abs(value)
            )
    component_statistics = tuple(
        ComponentStatistics(
            name=name,
            l2_norm=accumulator.result(),
            linf_norm=linf,
        )
        for name, accumulator, linf in zip(
            field.component_names,
            component_accumulators,
            component_linf,
            strict=True,
        )
    )
    return FieldStatistics(
        l2_norm=aggregate.result(),
        linf_norm=max(component_linf),
        checksum=_field_checksum(field),
        components=component_statistics,
    )


def _tolerance(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise FieldComparisonError(f"{name} must be a finite nonnegative number")
    result = float(value)
    if not math.isfinite(result) or result < 0.0:
        raise FieldComparisonError(f"{name} must be a finite nonnegative number")
    return result


def _allowed_difference(
    reference: float, candidate: float, absolute: float, relative: float
) -> float:
    scale = max(abs(reference), abs(candidate))
    relative_part = relative * scale
    if not math.isfinite(relative_part):
        return sys.float_info.max
    allowed = absolute + relative_part
    return min(allowed, sys.float_info.max)


def _validate_comparable_grids(
    reference: RegularGridField, candidate: RegularGridField
) -> None:
    if reference.shape != candidate.shape:
        raise FieldComparisonError(
            f"field shapes differ: {reference.shape!r} != {candidate.shape!r}"
        )
    if set(reference.component_names) != set(candidate.component_names):
        raise FieldComparisonError(
            "field component sets differ: "
            f"{reference.component_names!r} != {candidate.component_names!r}"
        )
    for axis_name, reference_axis, candidate_axis in zip(
        reference.coordinate_names,
        reference.axes,
        candidate.axes,
        strict=True,
    ):
        for index, (reference_coordinate, candidate_coordinate) in enumerate(
            zip(reference_axis, candidate_axis, strict=True)
        ):
            if not _coordinate_matches(reference_coordinate, candidate_coordinate):
                raise FieldComparisonError(
                    f"field coordinate grids differ at {axis_name}[{index}]: "
                    f"{reference_coordinate:.17g} != {candidate_coordinate:.17g}"
                )


def compare_fields(
    reference: RegularGridField,
    candidate: RegularGridField,
    *,
    absolute_tolerance: float,
    relative_tolerance: float,
) -> FieldComparison:
    """Compare every named component using ``abs + rel * max(|a|, |b|)``."""

    if not isinstance(reference, RegularGridField) or not isinstance(
        candidate, RegularGridField
    ):
        raise FieldComparisonError("both fields must be parsed RegularGridField values")
    absolute = _tolerance(absolute_tolerance, "absolute_tolerance")
    relative = _tolerance(relative_tolerance, "relative_tolerance")
    _validate_comparable_grids(reference, candidate)

    candidate_by_name = {
        name: values
        for name, values in zip(
            candidate.component_names, candidate._component_values, strict=True
        )
    }
    difference_accumulators = tuple(
        _ScaledSumSquares() for _ in reference.component_names
    )
    aggregate_difference = _ScaledSumSquares()
    component_linf = [0.0 for _ in reference.component_names]
    component_mismatches = [0 for _ in reference.component_names]
    mismatch_count = 0
    worst_absolute: PointDifference | None = None
    worst_mismatch: PointDifference | None = None
    worst_mismatch_score: tuple[float, float] | None = None
    for flat_index, sqrt_weight in enumerate(_point_sqrt_weights(reference)):
        coordinates = reference.coordinates(flat_index)
        grid_index = reference.grid_index(flat_index)
        for component_index, (name, reference_values) in enumerate(
            zip(
                reference.component_names,
                reference._component_values,
                strict=True,
            )
        ):
            reference_value = reference_values[flat_index]
            candidate_value = candidate_by_name[name][flat_index]
            difference = abs(candidate_value - reference_value)
            if not math.isfinite(difference):
                raise FieldComparisonError(
                    f"field difference overflow at component {name!r}, "
                    f"coordinate {coordinates!r}"
                )
            scale = max(abs(reference_value), abs(candidate_value))
            relative_difference = 0.0 if difference == 0.0 else difference / scale
            allowed = _allowed_difference(
                reference_value, candidate_value, absolute, relative
            )
            mismatch = difference > allowed
            if mismatch:
                mismatch_count += 1
                component_mismatches[component_index] += 1
            difference_accumulators[component_index].add(difference * sqrt_weight)
            aggregate_difference.add(difference * sqrt_weight)
            component_linf[component_index] = max(
                component_linf[component_index], difference
            )
            point = PointDifference(
                flat_index=flat_index,
                grid_index=grid_index,
                coordinates=coordinates,
                component=name,
                reference=reference_value,
                candidate=candidate_value,
                absolute_difference=difference,
                relative_difference=relative_difference,
                allowed_difference=allowed,
            )
            if (
                worst_absolute is None
                or difference > worst_absolute.absolute_difference
            ):
                worst_absolute = point
            if mismatch:
                normalized_violation = (
                    difference / allowed if allowed > 0.0 else math.inf
                )
                mismatch_score = (normalized_violation, difference)
                if (
                    worst_mismatch_score is None
                    or mismatch_score > worst_mismatch_score
                ):
                    worst_mismatch_score = mismatch_score
                    worst_mismatch = point

    if worst_absolute is None:  # Grid validation makes this unreachable.
        raise FieldComparisonError("cannot compare an empty field")
    worst = worst_mismatch if worst_mismatch is not None else worst_absolute
    component_differences = tuple(
        ComponentDifference(
            name=name,
            l2_difference=accumulator.result(),
            linf_difference=linf,
            mismatch_count=mismatches,
        )
        for name, accumulator, linf, mismatches in zip(
            reference.component_names,
            difference_accumulators,
            component_linf,
            component_mismatches,
            strict=True,
        )
    )
    return FieldComparison(
        passed=mismatch_count == 0,
        absolute_tolerance=absolute,
        relative_tolerance=relative,
        point_count=reference.point_count,
        component_count=len(reference.component_names),
        mismatch_count=mismatch_count,
        l2_difference=aggregate_difference.result(),
        linf_difference=max(component_linf),
        reference=summarize_field(reference),
        candidate=summarize_field(candidate),
        component_differences=component_differences,
        largest_difference_point=worst_absolute,
        worst_point=worst,
    )
