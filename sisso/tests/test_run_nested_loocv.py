from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from sisso.run_nested_loocv import make_sissopp_config, parse_dimensions
from sisso.workflow import read_prediction, regression_metrics


class NestedWorkflowTest(unittest.TestCase):
    def test_inner_config_excludes_outer_and_inner_materials(self) -> None:
        settings = {
            "property_key": "FWHM",
            "opset": ["add", "sub"],
            "max_rung": 0,
            "random_seed": 17,
        }
        config = make_sissopp_config(Path("data.csv"), 1, 20, [3, 8], settings)
        self.assertEqual(config["leave_out_inds"], [3, 8])
        self.assertEqual(config["max_rung"], 0)
        self.assertEqual(config["desc_dim"], 1)
        self.assertEqual(config["n_sis_select"], 20)

    def test_outer_config_excludes_only_outer_material(self) -> None:
        settings = {
            "property_key": "FWHM",
            "opset": [],
            "max_rung": 0,
            "random_seed": 17,
        }
        config = make_sissopp_config(Path("data.csv"), 2, 5, [4], settings)
        self.assertEqual(config["leave_out_inds"], [4])
        self.assertEqual(config["desc_dim"], 2)

    def test_duplicate_held_out_index_is_rejected(self) -> None:
        settings = {
            "property_key": "FWHM",
            "opset": [],
            "max_rung": 0,
            "random_seed": 17,
        }
        with self.assertRaisesRegex(ValueError, "unique"):
            make_sissopp_config(Path("data.csv"), 1, 5, [2, 2], settings)

    def test_model_prediction_parser(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "model.dat"
            path.write_text(
                "# c0 + a0 * f001\n"
                "# other metadata\n"
                "m00, 10.0, 11.5, 0.3\n"
                "m01, 20.0, 19.0, 0.8\n"
            )
            self.assertEqual(read_prediction(path, "m01"), (20.0, 19.0, "c0 + a0 * f001"))

    def test_metrics(self) -> None:
        result = regression_metrics([10.0, 20.0, 30.0], [11.0, 18.0, 31.0])
        self.assertAlmostEqual(result["RMSE_meV"], (6.0 / 3.0) ** 0.5)
        self.assertAlmostEqual(result["MAE_meV"], 4.0 / 3.0)
        self.assertAlmostEqual(result["Q2_LOO"], 0.97)

    def test_dimension_parser(self) -> None:
        self.assertEqual(parse_dimensions("1,2"), [1, 2])


if __name__ == "__main__":
    unittest.main()

