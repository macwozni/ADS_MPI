"""Reusable numerical validation primitives for ADS benchmark results."""

from .fields import (
    ComponentDifference,
    ComponentStatistics,
    FieldComparison,
    FieldComparisonError,
    FieldStatistics,
    PointDifference,
    RegularGridField,
    compare_fields,
    parse_regular_grid_csv,
    scalar_component,
    select_components,
    summarize_field,
)

__all__ = [
    "ComponentDifference",
    "ComponentStatistics",
    "FieldComparison",
    "FieldComparisonError",
    "FieldStatistics",
    "PointDifference",
    "RegularGridField",
    "compare_fields",
    "parse_regular_grid_csv",
    "scalar_component",
    "select_components",
    "summarize_field",
]
