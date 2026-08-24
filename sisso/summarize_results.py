#!/usr/bin/env python3
"""Summarize fold-wise SISSO++ selections without assuming a preferred formula."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

try:
    from .workflow import (
        dump_json,
        expression_feature_ids,
        load_json,
        read_csv,
        regression_metrics,
        write_csv,
    )
except ImportError:  # Direct execution: python sisso/summarize_results.py
    from workflow import (  # type: ignore
        dump_json,
        expression_feature_ids,
        load_json,
        read_csv,
        regression_metrics,
        write_csv,
    )


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_DATA = SCRIPT_DIR.parent / "data" / "sisso"
DEFAULT_INPUTS = SCRIPT_DIR / "work" / "inputs"
DEFAULT_RESULTS = SCRIPT_DIR / "work" / "results"
DEFAULT_OUTPUT = SCRIPT_DIR / "work" / "summary"


def formula_family(row: dict[str, str]) -> str:
    return "|".join(
        [
            row["local_expression"],
            f"within={row['within_reducer']}",
            f"across={row['across_reducer']}",
            f"outer={row['outer_transform']}",
        ]
    )


def summarize_dimension(
    dimension: int,
    prediction_path: Path,
    ancestry_by_id: dict[str, dict[str, str]],
) -> tuple[dict, list[dict]]:
    predictions = read_csv(prediction_path)
    computed = regression_metrics(
        [float(row["observed_FWHM_meV"]) for row in predictions],
        [float(row["predicted_FWHM_meV"]) for row in predictions],
    )
    family_rows: list[dict] = []
    family_sets: Counter[str] = Counter()
    for row in predictions:
        feature_ids = expression_feature_ids(row["model_expression"], ancestry_by_id)
        if len(feature_ids) != dimension:
            raise ValueError(
                f"fold {row['outer_fold']} in {prediction_path} contains "
                f"{len(feature_ids)} recognized features; expected {dimension}"
            )
        families = [formula_family(ancestry_by_id[feature_id]) for feature_id in feature_ids]
        family_set = " + ".join(sorted(families))
        family_sets[family_set] += 1
        family_rows.append(
            {
                "dimension": dimension,
                "outer_fold": row["outer_fold"],
                "held_out_material": row["material"],
                "feature_ids": " + ".join(feature_ids),
                "candidate_keys": " + ".join(
                    ancestry_by_id[feature_id]["candidate_key"] for feature_id in feature_ids
                ),
                "formula_families": family_set,
            }
        )
    summary = {
        "dimension": dimension,
        "n_outer_folds": len(predictions),
        **computed,
        "unique_formula_family_sets": len(family_sets),
        "most_frequent_formula_family_set": family_sets.most_common(1)[0][0],
        "most_frequent_formula_family_set_folds": family_sets.most_common(1)[0][1],
    }
    return summary, family_rows


def verify_reference(generated: dict, reference_path: Path, tolerance: float) -> None:
    reference = load_json(reference_path)
    for key in ["RMSE_meV", "MAE_meV", "Q2_LOO"]:
        difference = abs(float(generated[key]) - float(reference[key]))
        if difference > tolerance:
            raise ValueError(
                f"{key} differs from {reference_path}: {generated[key]} versus {reference[key]}"
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ancestry", type=Path, default=DEFAULT_DATA / "candidate_ancestry.csv")
    parser.add_argument("--inputs-dir", type=Path, default=DEFAULT_INPUTS)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--reference-dir", type=Path, default=DEFAULT_DATA / "reference_results"
    )
    parser.add_argument("--verify-reference", action="store_true")
    parser.add_argument("--tolerance", type=float, default=1e-8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ancestry = read_csv(args.ancestry)
    ancestry_by_id = {row["feature_id"]: row for row in ancestry}
    if len(ancestry_by_id) != len(ancestry):
        raise ValueError("candidate ancestry contains duplicate feature IDs")
    preparation = load_json(args.inputs_dir / "preparation_summary.json")

    summaries = []
    all_family_rows: list[dict] = []
    for dimension in [1, 2]:
        summary, family_rows = summarize_dimension(
            dimension,
            args.results_dir / f"D{dimension}_outer_predictions.csv",
            ancestry_by_id,
        )
        if args.verify_reference:
            verify_reference(
                summary,
                args.reference_dir / f"D{dimension}_metrics.json",
                args.tolerance,
            )
        summaries.append(summary)
        all_family_rows.extend(family_rows)

    output = {
        "analysis_name": preparation["analysis_name"],
        "n_materials": preparation["n_materials"],
        "n_candidates": preparation["n_candidates"],
        "sissopp_commit": preparation["sissopp_commit"],
        "dimensions": summaries,
    }
    dump_json(args.output_dir / "summary.json", output)
    write_csv(args.output_dir / "selected_formula_families.csv", all_family_rows)
    print(f"Summary written to {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()

