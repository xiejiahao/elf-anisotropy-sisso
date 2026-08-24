from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class Bond:
    center_index: int
    iodine_index: int
    image: tuple[int, int, int]
    delta_frac: tuple[float, float, float]
    distance_A: float


@dataclass(frozen=True)
class Peak:
    first_min_A: float
    topological_max_A: float
    second_min_A: float
    prominence: float
    r_peak_A: float
    h_peak: float


@dataclass(frozen=True)
class BondResult:
    site_number: int
    center_index: int
    direction: int
    iodine_index: int
    image_a: int
    image_b: int
    image_c: int
    distance_A: float
    r_peak_A: float
    h_peak: float
    M: float
    first_min_A: float
    topological_max_A: float
    second_min_A: float
    prominence: float


@dataclass(frozen=True)
class SiteResult:
    site_number: int
    center_index: int
    n_iodine_neighbors: int
    sixth_neighbor_A: float
    seventh_neighbor_A: float
    shell_gap_A: float
    M_mean: float
    M_min: float
    M_max: float
    A_j: float
    A_j_squared: float
    is_most_anisotropic: bool


@dataclass(frozen=True)
class AnalysisResult:
    input_path: Path
    center_species: str
    site_count: int
    coordination_passed: bool
    most_anisotropic_site: int
    most_anisotropic_center_index: int
    A_ELF: float
    A_ELF_squared: float
    reference_category: str
    reference_statement: str
    physical_interpretation: str
    sites: tuple[SiteResult, ...]
    bonds: tuple[BondResult, ...]

    def summary_dict(self) -> dict[str, Any]:
        return {
            "schema_version": "1.0",
            "input_path": str(self.input_path),
            "center_species": self.center_species,
            "site_count": self.site_count,
            "coordination_passed": self.coordination_passed,
            "most_anisotropic_site": self.most_anisotropic_site,
            "most_anisotropic_center_index": self.most_anisotropic_center_index,
            "A_ELF": self.A_ELF,
            "A_ELF_squared": self.A_ELF_squared,
            "reference_category": self.reference_category,
            "reference_statement": self.reference_statement,
            "physical_interpretation": self.physical_interpretation,
        }

    def site_dicts(self) -> list[dict[str, Any]]:
        return [asdict(site) for site in self.sites]

    def bond_dicts(self) -> list[dict[str, Any]]:
        return [asdict(bond) for bond in self.bonds]
