from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from sisso.prepare_candidates import validate_and_prepare


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


class PrepareCandidatesTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)
        self.matrix = self.base / "candidate_matrix.csv"
        self.ancestry = self.base / "candidate_ancestry.csv"
        self.targets = self.base / "targets.csv"
        self.sample_map = self.base / "sample_map.csv"
        self.config = self.base / "search.json"
        self.output = self.base / "output"
        self.config.write_text(
            json.dumps(
                {
                    "analysis_name": "test",
                    "candidate_space": {
                        "expected_materials": 3,
                        "expected_candidates": 2,
                        "q_local": 1,
                    },
                    "sissopp": {"commit": "abc"},
                }
            )
        )
        write_csv(
            self.matrix,
            ["material", "f001", "f002"],
            [
                {"material": "a", "f001": 1.0, "f002": 3.0},
                {"material": "b", "f001": 2.0, "f002": 2.0},
                {"material": "c", "f001": 3.0, "f002": 1.0},
            ],
        )
        ancestry_rows = []
        for feature_id in ["f001", "f002"]:
            ancestry_rows.append(
                {
                    "feature_id": feature_id,
                    "candidate_key": f"key:{feature_id}",
                    "required_modules": "A",
                    "object_type": "sn_i_bond",
                    "local_expression": "mult(h,r)",
                    "within_reducer": "normalized_range",
                    "across_reducer": "max",
                    "outer_transform": "sq",
                    "unit": "1",
                    "grid_stable": "True",
                    "is_family_representative": "True",
                }
            )
        write_csv(self.ancestry, list(ancestry_rows[0]), ancestry_rows)
        write_csv(
            self.targets,
            ["material", "FWHM_meV"],
            [
                {"material": "a", "FWHM_meV": 10},
                {"material": "b", "FWHM_meV": 20},
                {"material": "c", "FWHM_meV": 30},
            ],
        )
        write_csv(
            self.sample_map,
            ["sample_id", "material"],
            [
                {"sample_id": "m00", "material": "a"},
                {"sample_id": "m01", "material": "b"},
                {"sample_id": "m02", "material": "c"},
            ],
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_preparation_preserves_all_candidates_without_target_screening(self) -> None:
        summary = validate_and_prepare(
            self.matrix,
            self.ancestry,
            self.targets,
            self.sample_map,
            self.config,
            self.output,
        )
        self.assertEqual(summary["n_candidates"], 2)
        self.assertFalse(summary["target_used_for_candidate_filtering"])
        with (self.output / "data.csv").open(newline="") as handle:
            header = next(csv.reader(handle))
        self.assertEqual(header, ["sample", "FWHM (meV)", "f001", "f002"])

    def test_rejects_non_dimensionless_candidate(self) -> None:
        with self.ancestry.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        rows[0]["unit"] = "A"
        write_csv(self.ancestry, list(rows[0]), rows)
        with self.assertRaisesRegex(ValueError, "not dimensionless"):
            validate_and_prepare(
                self.matrix,
                self.ancestry,
                self.targets,
                self.sample_map,
                self.config,
                self.output,
            )


if __name__ == "__main__":
    unittest.main()
