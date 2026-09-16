"""Small deterministic registries used as the framework composition API."""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Generic, Iterable, TypeVar

from .errors import RegistryError
from .model import BuildProfileDefinition, ExactCaseDefinition, FamilyDefinition
from .protocols import Launcher, ProblemAdapter


T = TypeVar("T")
_NAME = re.compile(r"^[a-z][a-z0-9_-]*$")


class Registry(Generic[T]):
    """Name-indexed collection with duplicate and spelling protection."""

    def __init__(self, kind: str) -> None:
        self.kind = kind
        self._items: dict[str, T] = {}

    def register(self, name: str, item: T) -> None:
        if not isinstance(name, str) or not _NAME.fullmatch(name):
            raise RegistryError(f"invalid {self.kind} name: {name!r}")
        if name in self._items:
            raise RegistryError(f"duplicate {self.kind}: {name}")
        self._items[name] = item

    def get(self, name: str) -> T:
        try:
            return self._items[name]
        except KeyError as error:
            available = ", ".join(self.names()) or "<none>"
            raise RegistryError(
                f"unknown {self.kind}: {name}; available: {available}"
            ) from error

    def contains(self, name: str) -> bool:
        return name in self._items

    def names(self) -> tuple[str, ...]:
        return tuple(sorted(self._items))

    def values(self) -> Iterable[T]:
        for name in self.names():
            yield self._items[name]


@dataclass(frozen=True)
class Catalog:
    """All extension axes supplied to planners and executors by injection."""

    adapters: Registry[ProblemAdapter]
    families: Registry[FamilyDefinition]
    exact_cases: Registry[ExactCaseDefinition]
    build_profiles: Registry[BuildProfileDefinition]
    launchers: Registry[Launcher]
    schemes: Registry[str]
