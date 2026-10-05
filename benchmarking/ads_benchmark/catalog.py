"""Composition root for built-in benchmark components.

This is intentionally the only production module that knows both the neutral
framework and concrete adapters.  Adding an adapter means implementing the
protocol and registering it here; planner and executor code stay unchanged.
"""

from __future__ import annotations

from .analysis.spatial import HConvergenceAnalyzer, PConvergenceAnalyzer
from .analysis.strong import StrongScalingAnalyzer
from .analysis.temporal import TemporalConvergenceAnalyzer
from .analysis.validation import FieldValidationAnalyzer
from .analysis.weak import WeakScalingAnalyzer
from .components.manufactured import ManufacturedTransientAdapter
from .components.planning import DirectLauncher, default_mpi_launcher
from .framework.model import (
    BuildProfileDefinition,
    ExactCaseDefinition,
    FamilyDefinition,
)
from .framework.registry import Catalog, Registry


def build_catalog(
    *,
    available_mpi_slots: int | None = None,
    available_cpu_slots: int | None = None,
    launcher_template: str | None = None,
) -> Catalog:
    adapters = Registry("problem")
    for name in ("igrm_l2", "igrm_heat", "pure_diffusion_igrm"):
        adapters.register(name, ManufacturedTransientAdapter(name=name))

    families = Registry("experiment family")
    for name, description, analyzer in (
        (
            "temporal",
            "fixed-space temporal convergence",
            "temporal-convergence",
        ),
        ("h", "mesh-size convergence", "h-convergence"),
        ("p", "polynomial-degree convergence", "p-convergence"),
        (
            "validation",
            "MPI/OpenMP full-field validation",
            "field-validation",
        ),
        ("strong", "strong scaling", "strong-scaling"),
        ("weak", "weak scaling", "weak-scaling"),
    ):
        families.register(name, FamilyDefinition(name, description, analyzer))

    exact_cases = Registry("exact case")
    exact_cases.register(
        "temporal-polynomial",
        ExactCaseDefinition(
            "temporal-polynomial",
            "exp(-t) times the exactly representable cubic Neumann field",
        ),
    )
    exact_cases.register(
        "spatial-cosine",
        ExactCaseDefinition(
            "spatial-cosine",
            "exp(-t) times a non-polynomial cosine Neumann field",
        ),
    )

    build_profiles = Registry("build profile")
    for name in ("debug", "release"):
        build_profiles.register(
            name, BuildProfileDefinition(name, f"ADS {name} build")
        )

    launchers = Registry("launcher")
    direct = DirectLauncher()
    mpi = default_mpi_launcher(
        available_mpi_slots=available_mpi_slots,
        available_cpu_slots=available_cpu_slots,
        command_template=launcher_template,
    )
    launchers.register(direct.name, direct)
    launchers.register(mpi.name, mpi)

    schemes = Registry("scheme")
    for name in ("dg", "pr", "be"):
        schemes.register(name, name)

    analyzers = Registry("analyzer")
    temporal_analyzer = TemporalConvergenceAnalyzer()
    analyzers.register(temporal_analyzer.name, temporal_analyzer)
    h_analyzer = HConvergenceAnalyzer()
    analyzers.register(h_analyzer.name, h_analyzer)
    p_analyzer = PConvergenceAnalyzer()
    analyzers.register(p_analyzer.name, p_analyzer)
    validation_analyzer = FieldValidationAnalyzer()
    analyzers.register(validation_analyzer.name, validation_analyzer)
    strong_analyzer = StrongScalingAnalyzer()
    analyzers.register(strong_analyzer.name, strong_analyzer)
    weak_analyzer = WeakScalingAnalyzer()
    analyzers.register(weak_analyzer.name, weak_analyzer)

    return Catalog(
        adapters=adapters,
        families=families,
        exact_cases=exact_cases,
        build_profiles=build_profiles,
        launchers=launchers,
        schemes=schemes,
        analyzers=analyzers,
    )
