from __future__ import annotations

import csv
import json
from pathlib import Path

from .models import AnalysisResult


def format_text(result: AnalysisResult) -> str:
    return "\n".join(
        [
            f"Input: {result.input_path}",
            f"Center species: {result.center_species}",
            f"{result.center_species} sites: {result.site_count}",
            "Coordination check: passed (six periodic I neighbors per metal site)",
            (
                f"Most anisotropic site: {result.center_species} site "
                f"{result.most_anisotropic_site} "
                f"(structure index {result.most_anisotropic_center_index})"
            ),
            "",
            f"A_ELF   = {result.A_ELF:.8f}",
            f"A_ELF^2 = {result.A_ELF_squared:.10f}",
            "",
            f"Reference interpretation: {result.reference_statement}.",
            f"Physical interpretation: {result.physical_interpretation}",
        ]
    )


def format_json(result: AnalysisResult) -> str:
    return json.dumps(result.summary_dict(), indent=2, sort_keys=True)


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows are available for {path.name}.")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(result: AnalysisResult, output_dir: str | Path) -> Path:
    destination = Path(output_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    with (destination / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(result.summary_dict(), handle, indent=2, sort_keys=True)
        handle.write("\n")
    _write_csv(destination / "per_site.csv", result.site_dicts())
    _write_csv(destination / "per_bond.csv", result.bond_dicts())
    return destination
