from __future__ import annotations

import unittest

from ads_benchmark.analysis.base import AnalysisError
from ads_benchmark.analysis.statistics import (
    scaling_metrics,
    summarize_samples,
    weak_scaling_efficiency,
)


class ScalingStatisticsTests(unittest.TestCase):
    def test_median_range_and_mad_preserve_odd_sample_set(self) -> None:
        statistics = summarize_samples((8.0, 1.0, 4.0, 2.0, 16.0, 3.0, 5.0))
        self.assertEqual(statistics.count, 7)
        self.assertEqual(statistics.minimum, 1.0)
        self.assertEqual(statistics.maximum, 16.0)
        self.assertEqual(statistics.median, 4.0)
        self.assertEqual(statistics.median_absolute_deviation, 2.0)
        self.assertEqual(statistics.relative_median_absolute_deviation, 0.5)

    def test_speedup_and_efficiency_use_relative_resource_count(self) -> None:
        metrics = scaling_metrics(
            baseline_seconds=12.0,
            measured_seconds=4.0,
            baseline_resources=2,
            resources=8,
        )
        self.assertEqual(metrics.speedup, 3.0)
        self.assertEqual(metrics.efficiency, 0.75)

    def test_weak_efficiency_is_only_the_baseline_time_ratio(self) -> None:
        self.assertEqual(
            weak_scaling_efficiency(
                baseline_seconds=12.0, measured_seconds=15.0
            ),
            0.8,
        )

    def test_invalid_samples_and_scaling_inputs_are_rejected(self) -> None:
        for samples in (
            (),
            (-1.0,),
            (float("nan"),),
            (float("inf"),),
            (True,),
            ("1.0",),
        ):
            with self.subTest(samples=samples), self.assertRaises(AnalysisError):
                summarize_samples(samples)
        with self.assertRaises(AnalysisError):
            scaling_metrics(
                baseline_seconds=1.0,
                measured_seconds=0.0,
                baseline_resources=1,
                resources=2,
            )
        with self.assertRaises(AnalysisError):
            scaling_metrics(
                baseline_seconds=1.0,
                measured_seconds=0.5,
                baseline_resources=4,
                resources=2,
            )
        for invalid_time in (True, "1.0", float("nan"), float("inf")):
            with self.subTest(invalid_time=invalid_time), self.assertRaises(
                AnalysisError
            ):
                scaling_metrics(
                    baseline_seconds=invalid_time,  # type: ignore[arg-type]
                    measured_seconds=1.0,
                    baseline_resources=1,
                    resources=2,
                )
            with self.subTest(weak_invalid_time=invalid_time), self.assertRaises(
                AnalysisError
            ):
                weak_scaling_efficiency(
                    baseline_seconds=invalid_time,  # type: ignore[arg-type]
                    measured_seconds=1.0,
                )


if __name__ == "__main__":
    unittest.main(verbosity=2)
