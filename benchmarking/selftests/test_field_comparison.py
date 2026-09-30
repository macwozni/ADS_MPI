from __future__ import annotations

import csv
import io
import math
import unittest

from ads_benchmark.validation.fields import (
    FieldComparisonError,
    compare_fields,
    parse_regular_grid_csv,
    scalar_component,
    select_components,
    summarize_field,
)


def _field_csv(
    axes: tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...]],
    component_names: tuple[str, ...],
    value_at,
) -> str:
    rows: list[list[object]] = [["x", "y", "z", *component_names]]
    for iz, z in enumerate(axes[2]):
        for iy, y in enumerate(axes[1]):
            for ix, x in enumerate(axes[0]):
                values = value_at(ix, iy, iz, x, y, z)
                rows.append([x, y, z, *values])
    stream = io.StringIO(newline="")
    csv.writer(stream, lineterminator="\n").writerows(rows)
    return stream.getvalue()


UNIT_AXES = ((0.0, 0.5, 1.0),) * 3


class _Artifact:
    def __init__(self, text: str) -> None:
        self.text = text
        self.read_count = 0

    def read_text(self) -> str:
        self.read_count += 1
        return self.text


class FieldComparisonTests(unittest.TestCase):
    def test_zero_copy_selection_and_aliases_compare_numerical_to_exact(
        self,
    ) -> None:
        text = _field_csv(
            UNIT_AXES,
            ("numerical", "exact", "error"),
            lambda _ix, _iy, _iz, x, y, z: (
                x + y + z,
                x + y + z,
                0.0,
            ),
        )
        field = parse_regular_grid_csv(text, shape=3)
        selected = select_components(
            field,
            ("numerical", "error"),
            label="selected",
        )
        numerical = scalar_component(field, "numerical", label="numerical")
        exact = select_components(field, {"exact": "value"}, label="exact")

        self.assertEqual(field.component_names, ("numerical", "exact", "error"))
        self.assertEqual(selected.component_names, ("numerical", "error"))
        self.assertEqual(selected.label, "selected")
        self.assertIs(selected.axes, field.axes)
        self.assertIs(selected._component_values[0], field._component_values[0])
        self.assertEqual(numerical.component_names, ("value",))
        self.assertTrue(
            compare_fields(
                numerical,
                exact,
                absolute_tolerance=0.0,
                relative_tolerance=0.0,
            ).passed
        )
        with self.assertRaisesRegex(FieldComparisonError, "result.*unique"):
            select_components(field, {"numerical": "value", "exact": "value"})

    def test_artifact_text_is_parsed_and_component_order_is_normalized(self) -> None:
        first_text = _field_csv(
            UNIT_AXES,
            ("u", "tau"),
            lambda _ix, _iy, _iz, x, y, z: (x + 2.0 * y + z, x - z),
        )
        second_text = _field_csv(
            UNIT_AXES,
            ("tau", "u"),
            lambda _ix, _iy, _iz, x, y, z: (x - z, x + 2.0 * y + z),
        )
        artifact = _Artifact(first_text)
        reference = parse_regular_grid_csv(
            artifact, shape=3, label="reference"
        )
        candidate = parse_regular_grid_csv(
            second_text, shape=(3, 3, 3), label="candidate"
        )

        self.assertEqual(artifact.read_count, 1)
        self.assertEqual(reference.point_count, 27)
        self.assertEqual(reference.coordinates(5), (1.0, 0.5, 0.0))
        self.assertTrue(reference.component_values("u").readonly)
        comparison = compare_fields(
            reference,
            candidate,
            absolute_tolerance=0.0,
            relative_tolerance=0.0,
        )
        self.assertTrue(comparison.passed)
        self.assertEqual(comparison.mismatch_count, 0)
        self.assertEqual(comparison.l2_difference, 0.0)
        self.assertEqual(comparison.linf_difference, 0.0)
        self.assertEqual(
            comparison.reference.checksum, comparison.candidate.checksum
        )
        self.assertEqual(comparison.to_dict()["comparison_count"], 54)

    def test_absolute_plus_relative_tolerance_is_applied_pointwise(self) -> None:
        reference_text = _field_csv(
            UNIT_AXES,
            ("u",),
            lambda ix, iy, iz, _x, _y, _z: (1000.0 + ix + iy + iz,),
        )
        candidate_text = _field_csv(
            UNIT_AXES,
            ("u",),
            lambda ix, iy, iz, _x, _y, _z: (
                1000.0
                + ix
                + iy
                + iz
                + (5.0e-4 if (ix, iy, iz) == (1, 1, 1) else 0.0),
            ),
        )
        reference = parse_regular_grid_csv(reference_text, shape=3)
        candidate = parse_regular_grid_csv(candidate_text, shape=3)

        accepted = compare_fields(
            reference,
            candidate,
            absolute_tolerance=1.0e-8,
            relative_tolerance=1.0e-6,
        )
        rejected = compare_fields(
            reference,
            candidate,
            absolute_tolerance=1.0e-8,
            relative_tolerance=1.0e-8,
        )

        self.assertTrue(accepted.passed)
        self.assertFalse(rejected.passed)
        self.assertEqual(rejected.mismatch_count, 1)
        self.assertEqual(rejected.worst_point.grid_index, (1, 1, 1))
        self.assertEqual(rejected.worst_point.coordinates, (0.5, 0.5, 0.5))
        self.assertEqual(rejected.worst_point.component, "u")
        self.assertIn("component='u'", rejected.diagnostic())
        self.assertIn("coordinate=(0.5, 0.5, 0.5)", rejected.diagnostic())

    def test_local_permutation_is_detected_despite_same_sum_and_l2_norm(self) -> None:
        def reference_value(ix, iy, iz, _x, _y, _z):
            if (ix, iy, iz) == (1, 1, 0):
                return (1.0,)
            if (ix, iy, iz) == (1, 1, 2):
                return (2.0,)
            return (0.0,)

        def candidate_value(ix, iy, iz, _x, _y, _z):
            if (ix, iy, iz) == (1, 1, 0):
                return (2.0,)
            if (ix, iy, iz) == (1, 1, 2):
                return (1.0,)
            return (0.0,)

        reference = parse_regular_grid_csv(
            _field_csv(UNIT_AXES, ("u",), reference_value), shape=3
        )
        candidate = parse_regular_grid_csv(
            _field_csv(UNIT_AXES, ("u",), candidate_value), shape=3
        )
        reference_statistics = summarize_field(reference)
        candidate_statistics = summarize_field(candidate)

        # A simple sum/checksum and the global L2 norm cannot see the swap.
        self.assertEqual(
            sum(reference.component_values("u")),
            sum(candidate.component_values("u")),
        )
        self.assertEqual(
            reference_statistics.l2_norm, candidate_statistics.l2_norm
        )
        comparison = compare_fields(
            reference,
            candidate,
            absolute_tolerance=1.0e-14,
            relative_tolerance=1.0e-14,
        )
        self.assertFalse(comparison.passed)
        self.assertEqual(comparison.mismatch_count, 2)
        self.assertEqual(comparison.linf_difference, 1.0)
        self.assertGreater(comparison.l2_difference, 0.0)
        self.assertNotEqual(
            comparison.reference.checksum, comparison.candidate.checksum
        )

    def test_multiple_components_report_the_actual_worst_component(self) -> None:
        reference = parse_regular_grid_csv(
            _field_csv(
                UNIT_AXES,
                ("u", "tau"),
                lambda _ix, _iy, _iz, _x, _y, _z: (1.0, 2.0),
            ),
            shape=3,
        )
        candidate = parse_regular_grid_csv(
            _field_csv(
                UNIT_AXES,
                ("u", "tau"),
                lambda ix, iy, iz, _x, _y, _z: (
                    1.0 + (0.1 if (ix, iy, iz) == (0, 0, 0) else 0.0),
                    2.0 + (0.75 if (ix, iy, iz) == (2, 1, 0) else 0.0),
                ),
            ),
            shape=3,
        )
        comparison = compare_fields(
            reference,
            candidate,
            absolute_tolerance=1.0e-12,
            relative_tolerance=0.0,
        )

        self.assertFalse(comparison.passed)
        self.assertEqual(comparison.mismatch_count, 2)
        self.assertEqual(comparison.worst_point.component, "tau")
        self.assertEqual(comparison.worst_point.grid_index, (2, 1, 0))
        by_name = {item.name: item for item in comparison.component_differences}
        self.assertAlmostEqual(by_name["u"].linf_difference, 0.1)
        self.assertAlmostEqual(by_name["tau"].linf_difference, 0.75)

    def test_failure_reports_largest_normalized_tolerance_violation(self) -> None:
        def reference_value(ix, iy, iz, _x, _y, _z):
            if (ix, iy, iz) == (1, 0, 0):
                return (1000.0,)
            return (0.0,)

        def candidate_value(ix, iy, iz, _x, _y, _z):
            if (ix, iy, iz) == (0, 0, 0):
                return (1.5e-8,)
            if (ix, iy, iz) == (1, 0, 0):
                return (1000.0005,)
            if (ix, iy, iz) == (2, 0, 0):
                return (3.0e-8,)
            return (0.0,)

        reference = parse_regular_grid_csv(
            _field_csv(UNIT_AXES, ("u",), reference_value), shape=3
        )
        candidate = parse_regular_grid_csv(
            _field_csv(UNIT_AXES, ("u",), candidate_value), shape=3
        )

        rejected = compare_fields(
            reference,
            candidate,
            absolute_tolerance=1.0e-8,
            relative_tolerance=1.0e-6,
        )
        accepted = compare_fields(
            reference,
            candidate,
            absolute_tolerance=1.0e-3,
            relative_tolerance=0.0,
        )

        self.assertFalse(rejected.passed)
        self.assertEqual(rejected.mismatch_count, 2)
        self.assertEqual(rejected.worst_point.grid_index, (2, 0, 0))
        self.assertAlmostEqual(rejected.worst_point.absolute_difference, 3.0e-8)
        self.assertAlmostEqual(rejected.linf_difference, 5.0e-4)
        self.assertEqual(
            rejected.largest_difference_point.grid_index, (1, 0, 0)
        )
        self.assertAlmostEqual(
            rejected.largest_difference_point.absolute_difference, 5.0e-4
        )
        self.assertTrue(accepted.passed)
        self.assertEqual(accepted.worst_point.grid_index, (1, 0, 0))
        self.assertEqual(
            accepted.largest_difference_point, accepted.worst_point
        )
        self.assertAlmostEqual(accepted.worst_point.absolute_difference, 5.0e-4)

    def test_parser_rejects_nonfinite_values_even_in_unselected_columns(self) -> None:
        text = _field_csv(
            UNIT_AXES,
            ("numerical", "exact"),
            lambda ix, iy, iz, _x, _y, _z: (
                1.0,
                "nan" if (ix, iy, iz) == (1, 1, 1) else 1.0,
            ),
        )
        with self.assertRaisesRegex(FieldComparisonError, "must be finite"):
            parse_regular_grid_csv(text, shape=3, components=("numerical",))

        infinity = text.replace("nan", "inf")
        with self.assertRaisesRegex(FieldComparisonError, "must be finite"):
            parse_regular_grid_csv(infinity, shape=3)

    def test_parser_rejects_wrong_count_order_header_and_irregular_grid(self) -> None:
        valid = _field_csv(
            UNIT_AXES,
            ("u",),
            lambda _ix, _iy, _iz, x, y, z: (x + y + z,),
        )
        with self.subTest("row count"):
            truncated = "\n".join(valid.splitlines()[:-1]) + "\n"
            with self.assertRaisesRegex(FieldComparisonError, "26 data rows"):
                parse_regular_grid_csv(truncated, shape=3)
        with self.subTest("canonical order"):
            rows = list(csv.reader(io.StringIO(valid)))
            rows[2][0] = "0.75"
            stream = io.StringIO(newline="")
            csv.writer(stream, lineterminator="\n").writerows(rows)
            with self.assertRaisesRegex(FieldComparisonError, "canonical"):
                parse_regular_grid_csv(stream.getvalue(), shape=3)
        with self.subTest("duplicate header"):
            duplicate = valid.replace("x,y,z,u", "x,y,z,u,u", 1)
            with self.assertRaisesRegex(FieldComparisonError, "unique"):
                parse_regular_grid_csv(duplicate, shape=3)
        with self.subTest("irregular axis"):
            axes = ((0.0, 0.2, 0.7, 1.0), (0.0, 1.0), (0.0, 1.0))
            irregular = _field_csv(
                axes,
                ("u",),
                lambda _ix, _iy, _iz, x, y, z: (x + y + z,),
            )
            with self.assertRaisesRegex(FieldComparisonError, "not regular"):
                parse_regular_grid_csv(irregular, shape=(4, 2, 2))

    def test_comparison_rejects_incompatible_grids_and_invalid_tolerances(self) -> None:
        reference = parse_regular_grid_csv(
            _field_csv(
                UNIT_AXES,
                ("u",),
                lambda _ix, _iy, _iz, x, y, z: (x + y + z,),
            ),
            shape=3,
        )
        other_axes = ((0.0, 1.0, 2.0),) * 3
        candidate = parse_regular_grid_csv(
            _field_csv(
                other_axes,
                ("u",),
                lambda _ix, _iy, _iz, x, y, z: (x + y + z,),
            ),
            shape=3,
        )
        with self.assertRaisesRegex(FieldComparisonError, "coordinate grids differ"):
            compare_fields(
                reference,
                candidate,
                absolute_tolerance=0.0,
                relative_tolerance=0.0,
            )
        with self.assertRaisesRegex(FieldComparisonError, "finite nonnegative"):
            compare_fields(
                reference,
                reference,
                absolute_tolerance=math.nan,
                relative_tolerance=0.0,
            )


if __name__ == "__main__":
    unittest.main()
