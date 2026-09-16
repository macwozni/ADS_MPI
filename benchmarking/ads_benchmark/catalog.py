"""Composition root for built-in benchmark components.

This is intentionally the only production module that knows both the neutral
framework and concrete adapters.  Adding an adapter means implementing the
protocol and registering it here; planner and executor code stay unchanged.
"""

from __future__ import annotations

from .components.planning import (
    DirectLauncher,
    PlanningAdapter,
    default_mpi_launcher,
)
from .framework.model import (
    BuildProfileDefinition,
    ExactCaseDefinition,
    FamilyDefinition,
)
from .framework.registry import Catalog, Registry


def build_catalog() -> Catalog:
    adapters = Registry("problem")
    for name in ("igrm_l2", "igrm_heat", "pure_diffusion_igrm"):
        adapters.register(name, PlanningAdapter(name=name))

    families = Registry("experiment family")
    for name, description in (
        ("temporal", "fixed-space temporal convergence"),
        ("h", "mesh-size convergence"),
        ("p", "polynomial-degree convergence"),
        ("validation", "MPI/OpenMP full-field validation"),
        ("strong", "strong scaling"),
        ("weak", "weak scaling"),
    ):
        families.register(name, FamilyDefinition(name, description))

    exact_cases = Registry("exact case")
    exact_cases.register(
        "temporal-polynomial",
        ExactCaseDefinition(
            "temporal-polynomial",
            "exp(-t) times the exactly representable cubic Neumann field",
        ),
    )

    build_profiles = Registry("build profile")
    for name in ("debug", "release"):
        build_profiles.register(
            name, BuildProfileDefinition(name, f"ADS {name} build")
        )

    launchers = Registry("launcher")
    direct = DirectLauncher()
    mpi = default_mpi_launcher()
    launchers.register(direct.name, direct)
    launchers.register(mpi.name, mpi)

    schemes = Registry("scheme")
    for name in ("dg", "pr", "be"):
        schemes.register(name, name)

    return Catalog(
        adapters=adapters,
        families=families,
        exact_cases=exact_cases,
        build_profiles=build_profiles,
        launchers=launchers,
        schemes=schemes,
    )
