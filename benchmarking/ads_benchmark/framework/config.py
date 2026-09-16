"""Strict, dependency-free JSON profile loading."""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Iterable

from .errors import ConfigurationError, RegistryError
from .model import MeasurementSpec, MpiSpec, TimeSpec, Vector3
from .registry import Registry


PROFILE_SCHEMA_VERSION = 1
_PROFILE_KEYS = {
    "schema_version",
    "name",
    "description",
    "family",
    "exact_cases",
    "problems",
    "schemes",
    "time_discretizations",
    "meshes",
    "degree_pairs",
    "process_layouts",
    "thread_counts",
    "execution",
    "build_profiles",
    "launcher",
}


@dataclass(frozen=True)
class ProfileDefinition:
    name: str
    description: str
    family: str
    exact_cases: tuple[str, ...]
    problems: tuple[str, ...]
    schemes: tuple[str, ...]
    time_discretizations: tuple[TimeSpec, ...]
    meshes: tuple[Vector3, ...]
    degree_pairs: tuple[tuple[Vector3, Vector3], ...]
    process_layouts: tuple[MpiSpec, ...]
    thread_counts: tuple[int, ...]
    measurement: MeasurementSpec
    build_profiles: tuple[str, ...]
    launcher: str


def _reject_constant(value: str) -> None:
    raise ConfigurationError(f"non-finite JSON number is not allowed: {value}")


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ConfigurationError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def load_json(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as stream:
            document = json.load(
                stream,
                parse_float=Decimal,
                parse_constant=_reject_constant,
                object_pairs_hook=_unique_object,
            )
    except ConfigurationError:
        raise
    except (OSError, json.JSONDecodeError) as error:
        raise ConfigurationError(f"cannot read profile {path}: {error}") from error
    if not isinstance(document, dict):
        raise ConfigurationError(f"profile {path} must contain one JSON object")
    return document


def _exact_keys(value: Any, required: set[str], field: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ConfigurationError(f"{field} must be an object")
    missing = sorted(required - set(value))
    extra = sorted(set(value) - required)
    if missing or extra:
        details: list[str] = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if extra:
            details.append("unknown " + ", ".join(extra))
        raise ConfigurationError(f"{field} has invalid keys: {'; '.join(details)}")
    return value


def _string(value: Any, field: str) -> str:
    if not isinstance(value, str) or not value:
        raise ConfigurationError(f"{field} must be a nonempty string")
    return value


def _integer(value: Any, field: str, *, minimum: int) -> int:
    if type(value) is not int or value < minimum:
        qualifier = "nonnegative" if minimum == 0 else "positive"
        raise ConfigurationError(f"{field} must be a {qualifier} integer")
    return value


def _decimal(value: Any, field: str, *, positive: bool = True) -> Decimal:
    if isinstance(value, bool):
        raise ConfigurationError(f"{field} must be a decimal number")
    try:
        parsed = value if isinstance(value, Decimal) else Decimal(str(value))
    except (InvalidOperation, TypeError, ValueError) as error:
        raise ConfigurationError(f"{field} must be a decimal number") from error
    if not parsed.is_finite():
        raise ConfigurationError(f"{field} must be finite")
    if positive and parsed <= 0:
        raise ConfigurationError(f"{field} must be positive")
    return parsed


def canonical_decimal(value: Decimal) -> str:
    if value == 0:
        return "0"
    # Decimal.normalize() obeys the active precision context and can round
    # values with more than 28 significant digits.  Formatting the original
    # value is context-free and preserves every input digit.
    text = format(value, "f")
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return text


def canonical_fraction(value: Fraction) -> str:
    """Return an exact decimal when finite, otherwise a reduced fraction."""

    numerator = value.numerator
    denominator = value.denominator
    reduced = denominator
    power_two = 0
    power_five = 0
    while reduced % 2 == 0:
        reduced //= 2
        power_two += 1
    while reduced % 5 == 0:
        reduced //= 5
        power_five += 1
    if reduced != 1:
        return f"{numerator}/{denominator}"

    scale = max(power_two, power_five)
    scaled = abs(numerator) * (2 ** (scale - power_two)) * (
        5 ** (scale - power_five)
    )
    if scale == 0:
        text = str(scaled)
    else:
        digits = str(scaled).zfill(scale + 1)
        text = f"{digits[:-scale]}.{digits[-scale:]}".rstrip("0").rstrip(".")
    if numerator < 0:
        text = "-" + text
    return text


def _list(value: Any, field: str) -> list[Any]:
    if not isinstance(value, list) or not value:
        raise ConfigurationError(f"{field} must be a nonempty array")
    return value


def _string_list(value: Any, field: str, *, lower: bool = False) -> tuple[str, ...]:
    result = []
    for index, entry in enumerate(_list(value, field)):
        item = _string(entry, f"{field}[{index}]")
        result.append(item.lower() if lower else item)
    return tuple(result)


def _vector3(value: Any, field: str) -> Vector3:
    if not isinstance(value, list) or len(value) != 3:
        raise ConfigurationError(f"{field} must contain exactly three integers")
    return tuple(
        _integer(entry, f"{field}[{index}]", minimum=1)
        for index, entry in enumerate(value)
    )  # type: ignore[return-value]


def _time_spec(value: Any, field: str) -> TimeSpec:
    if not isinstance(value, dict):
        raise ConfigurationError(f"{field} must be an object")
    keys = set(value)
    allowed = {"final_time", "steps", "time_step"}
    if "final_time" not in keys or not ({"steps", "time_step"} & keys):
        raise ConfigurationError(
            f"{field} requires final_time and steps and/or time_step"
        )
    extra = sorted(keys - allowed)
    if extra:
        raise ConfigurationError(f"{field} has unknown keys: {', '.join(extra)}")

    final_decimal = _decimal(value["final_time"], f"{field}.final_time")
    final_time = Fraction(*final_decimal.as_integer_ratio())
    steps = (
        _integer(value["steps"], f"{field}.steps", minimum=1)
        if "steps" in value
        else None
    )
    time_step_decimal = (
        _decimal(value["time_step"], f"{field}.time_step")
        if "time_step" in value
        else None
    )
    time_step = (
        Fraction(*time_step_decimal.as_integer_ratio())
        if time_step_decimal is not None
        else None
    )

    if steps is None:
        assert time_step is not None
        quotient = final_time / time_step
        if quotient.denominator != 1 or quotient <= 0:
            raise ConfigurationError(
                f"{field} final_time/time_step must be a positive integer"
            )
        steps = quotient.numerator
    elif time_step is None:
        time_step = final_time / steps
    elif time_step * steps != final_time:
        raise ConfigurationError(
            f"{field} final_time, time_step and steps are inconsistent"
        )

    assert time_step is not None
    return TimeSpec(
        final_time=canonical_fraction(final_time),
        time_step=canonical_fraction(time_step),
        steps=steps,
    )


def parse_profile(document: dict[str, Any], source: Path) -> ProfileDefinition:
    _exact_keys(document, _PROFILE_KEYS, str(source))
    if (
        type(document["schema_version"]) is not int
        or document["schema_version"] != PROFILE_SCHEMA_VERSION
    ):
        raise ConfigurationError(
            f"{source} schema_version must be {PROFILE_SCHEMA_VERSION}"
        )

    times = tuple(
        _time_spec(value, f"time_discretizations[{index}]")
        for index, value in enumerate(
            _list(document["time_discretizations"], "time_discretizations")
        )
    )
    meshes = tuple(
        _vector3(value, f"meshes[{index}]")
        for index, value in enumerate(_list(document["meshes"], "meshes"))
    )

    degree_pairs = []
    for index, raw_pair in enumerate(
        _list(document["degree_pairs"], "degree_pairs")
    ):
        pair = _exact_keys(
            raw_pair, {"test", "trial"}, f"degree_pairs[{index}]"
        )
        degree_pairs.append(
            (
                _vector3(pair["test"], f"degree_pairs[{index}].test"),
                _vector3(pair["trial"], f"degree_pairs[{index}].trial"),
            )
        )

    process_layouts = []
    for index, raw_layout in enumerate(
        _list(document["process_layouts"], "process_layouts")
    ):
        layout = _exact_keys(
            raw_layout, {"ranks", "grid"}, f"process_layouts[{index}]"
        )
        process_layouts.append(
            MpiSpec(
                ranks=_integer(
                    layout["ranks"], f"process_layouts[{index}].ranks", minimum=1
                ),
                process_grid=_vector3(
                    layout["grid"], f"process_layouts[{index}].grid"
                ),
            )
        )

    execution = _exact_keys(
        document["execution"],
        {"warmups", "samples", "timeout_seconds"},
        "execution",
    )
    measurement = MeasurementSpec(
        warmups=_integer(execution["warmups"], "execution.warmups", minimum=0),
        samples=_integer(execution["samples"], "execution.samples", minimum=1),
        timeout_seconds=canonical_decimal(
            _decimal(execution["timeout_seconds"], "execution.timeout_seconds")
        ),
    )

    thread_counts = tuple(
        _integer(value, f"thread_counts[{index}]", minimum=1)
        for index, value in enumerate(
            _list(document["thread_counts"], "thread_counts")
        )
    )
    return ProfileDefinition(
        name=_string(document["name"], "name"),
        description=_string(document["description"], "description"),
        family=_string(document["family"], "family"),
        exact_cases=_string_list(document["exact_cases"], "exact_cases"),
        problems=_string_list(document["problems"], "problems"),
        schemes=_string_list(document["schemes"], "schemes", lower=True),
        time_discretizations=times,
        meshes=meshes,
        degree_pairs=tuple(degree_pairs),
        process_layouts=tuple(process_layouts),
        thread_counts=thread_counts,
        measurement=measurement,
        build_profiles=_string_list(document["build_profiles"], "build_profiles"),
        launcher=_string(document["launcher"], "launcher"),
    )


def load_profiles(directory: Path) -> Registry[ProfileDefinition]:
    if not directory.is_dir():
        raise ConfigurationError(f"profile directory does not exist: {directory}")
    paths = sorted(directory.glob("*.json"))
    if not paths:
        raise ConfigurationError(f"profile directory contains no JSON files: {directory}")

    registry: Registry[ProfileDefinition] = Registry("profile")
    for path in paths:
        profile = parse_profile(load_json(path), path)
        try:
            registry.register(profile.name, profile)
        except RegistryError as error:
            raise ConfigurationError(f"{path}: {error}") from error
    return registry


def profiles_to_names(profiles: Iterable[ProfileDefinition]) -> tuple[str, ...]:
    return tuple(sorted(profile.name for profile in profiles))
