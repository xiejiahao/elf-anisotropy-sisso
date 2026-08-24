#!/usr/bin/env python3
"""Validate the published candidate bank and create a SISSO++ input table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

try:
    from .workflow import (
        dump_json,
        finite_float,
        load_json,
        parse_bool,
        read_csv,
        read_header,
        require_columns,
        require_unique,
        write_csv,
    )
except ImportError:  # Direct execution: python sisso/prepare_candidates.py
    from workflow import (  # type: ignore
        dump_json,
        finite_float,
        load_json,
        parse_bool,
        read_csv,
        read_header,
        require_columns,
        require_unique,
        write_csv,
    )


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = REPO_ROOT / "data" / "sisso"
DEFAULT_CONFIG = Path(__file__).resolve().parent / "configs" / "search.json"
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "work" / "inputs"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate the 1,539-candidate public search space and write SISSO++ data.csv."
    )
    parser.add_argument("--matrix", type=Path, default=DEFAULT_DATA / "candidate_matrix.csv")
    parser.add_argument("--ancestry", type=Path, default=DEFAULT_DATA / "candidate_ancestry.csv")
    parser.add_argument("--targets", type=Path, default=DEFAULT_DATA / "targets.csv")
    parser.add_argument("--sample-map", type=Path, default=DEFAULT_DATA / "sample_map.csv")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def validate_and_prepare(
    matrix_path: Path,
    ancestry_path: Path,
    targets_path: Path,
    sample_map_path: Path,
    config_path: Path,
    output_dir: Path,
) -> dict[str, object]:
    config = load_json(config_path)
    candidate_config = config["candidate_space"]
    expected_materials = int(candidate_config["expected_materials"])
    expected_candidates = int(candidate_config["expected_candidates"])

    matrix_header = read_header(matrix_path)
    require_columns(matrix_path, matrix_header, ["material"])
    feature_ids = [column for column in matrix_header if column not in {"material", "sample_id"}]
    require_unique(feature_ids, "candidate_matrix feature columns")
    if len(feature_ids) != expected_candidates:
        raise ValueError(
            f"{matrix_path} contains {len(feature_ids)} features; expected {expected_candidates}"
        )

    matrix_rows = read_csv(matrix_path)
    if len(matrix_rows) != expected_materials:
        raise ValueError(
            f"{matrix_path} contains {len(matrix_rows)} materials; expected {expected_materials}"
        )
    materials = [row["material"] for row in matrix_rows]
    require_unique(materials, "candidate_matrix material column")
    for row in matrix_rows:
        for feature_id in feature_ids:
            finite_float(row[feature_id], label=f"{row['material']}:{feature_id}")

    ancestry = read_csv(ancestry_path)
    ancestry_columns = read_header(ancestry_path)
    required_ancestry = [
        "feature_id",
        "candidate_key",
        "required_modules",
        "object_type",
        "local_expression",
        "within_reducer",
        "across_reducer",
        "outer_transform",
        "unit",
        "grid_stable",
        "is_family_representative",
    ]
    require_columns(ancestry_path, ancestry_columns, required_ancestry)
    ancestry_ids = [row["feature_id"] for row in ancestry]
    require_unique(ancestry_ids, "candidate_ancestry feature_id column")
    if set(ancestry_ids) != set(feature_ids):
        missing_metadata = sorted(set(feature_ids) - set(ancestry_ids))
        missing_values = sorted(set(ancestry_ids) - set(feature_ids))
        raise ValueError(
            "candidate matrix and ancestry feature IDs differ; "
            f"metadata missing for {len(missing_metadata)}, values missing for {len(missing_values)}"
        )
    for row in ancestry:
        if row["unit"] != "1":
            raise ValueError(f"candidate {row['feature_id']} is not dimensionless")
        if row["required_modules"] != "A":
            raise ValueError(f"candidate {row['feature_id']} is not an A-module feature")
        if not parse_bool(row["grid_stable"], field="grid_stable"):
            raise ValueError(f"candidate {row['feature_id']} is not grid-stable")
        if not parse_bool(row["is_family_representative"], field="is_family_representative"):
            raise ValueError(f"candidate {row['feature_id']} is not a family representative")

    targets = read_csv(targets_path)
    require_columns(targets_path, read_header(targets_path), ["material", "FWHM_meV"])
    target_by_material = {row["material"]: row for row in targets}
    if len(target_by_material) != len(targets):
        raise ValueError(f"{targets_path} contains duplicate materials")
    if set(target_by_material) != set(materials):
        raise ValueError("targets and candidate matrix contain different materials")
    for material, row in target_by_material.items():
        finite_float(row["FWHM_meV"], label=f"{material}:FWHM_meV")

    sample_rows = read_csv(sample_map_path)
    require_columns(sample_map_path, read_header(sample_map_path), ["sample_id", "material"])
    sample_by_material = {row["material"]: row["sample_id"] for row in sample_rows}
    if len(sample_by_material) != len(sample_rows):
        raise ValueError(f"{sample_map_path} contains duplicate materials")
    if set(sample_by_material) != set(materials):
        raise ValueError("sample map and candidate matrix contain different materials")
    require_unique(list(sample_by_material.values()), "sample_map sample_id column")
    if "sample_id" in matrix_header:
        for row in matrix_rows:
            if row["sample_id"] != sample_by_material[row["material"]]:
                raise ValueError(
                    f"candidate matrix sample_id does not match sample_map for {row['material']}"
                )
    target_header = read_header(targets_path)
    if "sample_id" in target_header:
        for row in targets:
            if row["sample_id"] != sample_by_material[row["material"]]:
                raise ValueError(f"target sample_id does not match sample_map for {row['material']}")

    output_dir.mkdir(parents=True, exist_ok=True)
    data_path = output_dir / "data.csv"
    with data_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample", "FWHM (meV)", *feature_ids])
        for row in matrix_rows:
            material = row["material"]
            writer.writerow(
                [sample_by_material[material], target_by_material[material]["FWHM_meV"]]
                + [row[feature_id] for feature_id in feature_ids]
            )

    ordered_samples = [
        {
            "sample_id": sample_by_material[material],
            "material": material,
            "FWHM_meV": target_by_material[material]["FWHM_meV"],
        }
        for material in materials
    ]
    write_csv(output_dir / "sample_map.csv", ordered_samples)
    ancestry_by_id = {row["feature_id"]: row for row in ancestry}
    write_csv(output_dir / "feature_map.csv", [ancestry_by_id[feature_id] for feature_id in feature_ids])

    summary: dict[str, object] = {
        "analysis_name": config["analysis_name"],
        "n_materials": len(materials),
        "n_candidates": len(feature_ids),
        "dimensionless": True,
        "modules": ["A"],
        "q_local": candidate_config["q_local"],
        "sissopp_commit": config["sissopp"]["commit"],
        "target_used_for_candidate_filtering": False,
    }
    dump_json(output_dir / "preparation_summary.json", summary)
    return summary


def main() -> None:
    args = parse_args()
    summary = validate_and_prepare(
        args.matrix,
        args.ancestry,
        args.targets,
        args.sample_map,
        args.config,
        args.output_dir,
    )
    print(f"Validated {summary['n_materials']} materials and {summary['n_candidates']} candidates.")
    print(f"Wrote SISSO++ inputs to {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
