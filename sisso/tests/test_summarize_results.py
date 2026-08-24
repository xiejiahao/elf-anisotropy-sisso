from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from sisso.summarize_results import formula_family, summarize_dimension


class SummarizeResultsTest(unittest.TestCase):
    def setUp(self) -> None:
        self.ancestry = {
            "f001": {
                "feature_id": "f001",
                "candidate_key": "key1",
                "local_expression": "mult(h,r)",
                "within_reducer": "normalized_range",
                "across_reducer": "max",
                "outer_transform": "sq",
            },
            "f002": {
                "feature_id": "f002",
                "candidate_key": "key2",
                "local_expression": "div(r,d)",
                "within_reducer": "cv",
                "across_reducer": "mean",
                "outer_transform": "identity",
            },
        }

    def test_family_label_uses_ancestry(self) -> None:
        self.assertEqual(
            formula_family(self.ancestry["f001"]),
            "mult(h,r)|within=normalized_range|across=max|outer=sq",
        )

    def test_summary_does_not_assume_a_feature_id(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "D1_outer_predictions.csv"
            rows = [
                {
                    "outer_fold": 0,
                    "sample_id": "m00",
                    "material": "a",
                    "observed_FWHM_meV": 10,
                    "predicted_FWHM_meV": 11,
                    "error_meV": 1,
                    "selected_n_sis": 5,
                    "model_expression": "c0 + a0 * f002",
                },
                {
                    "outer_fold": 1,
                    "sample_id": "m01",
                    "material": "b",
                    "observed_FWHM_meV": 20,
                    "predicted_FWHM_meV": 19,
                    "error_meV": -1,
                    "selected_n_sis": 5,
                    "model_expression": "c0 + a0 * f002",
                },
            ]
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            summary, selected = summarize_dimension(1, path, self.ancestry)
            self.assertEqual(summary["most_frequent_formula_family_set_folds"], 2)
            self.assertEqual(selected[0]["feature_ids"], "f002")


if __name__ == "__main__":
    unittest.main()

