from __future__ import annotations

from pathlib import Path

import numpy as np
from pymatgen.io.vasp.outputs import Chgcar

from .coordination import center_species, periodic_coordination
from .errors import DescriptorError
from .models import AnalysisResult, BondResult, SiteResult
from .profiles import DEFAULT_SAMPLES, sample_bond, topological_peak


def reference_category(A_ELF_squared: float) -> str:
    if A_ELF_squared < 0.001:
        return "weak"
    if A_ELF_squared < 0.03:
        return "moderate"
    return "strong"


def interpretation(species: str, category: str) -> tuple[str, str]:
    orbital = "5s2" if species == "Sn" else "6s2"
    reference = (
        f"{category} static ELF anisotropy relative to the 13-compound "
        "Sn DJ reference set"
    )
    physical = (
        f"The static structure shows {category} directional expression of the "
        f"{species}(II) {orbital} lone pair on this reference scale. A_ELF measures "
        "static directionality and does not describe the full finite-temperature "
        "dynamics of the lone pair."
    )
    return reference, physical


def analyze_elfcar(path: str | Path, *, samples: int = DEFAULT_SAMPLES) -> AnalysisResult:
    input_path = Path(path).expanduser().resolve()
    if not input_path.is_file():
        raise DescriptorError(f"ELFCAR file not found: {input_path}")
    try:
        chgcar = Chgcar.from_file(str(input_path))
    except Exception as exc:
        raise DescriptorError(f"Could not read ELFCAR: {exc}") from exc
    if "total" not in chgcar.data:
        raise DescriptorError("The ELFCAR does not contain a total ELF grid.")

    structure = chgcar.structure
    grid = np.asarray(chgcar.data["total"], dtype=float)
    if grid.ndim != 3 or not np.all(np.isfinite(grid)):
        raise DescriptorError("The total ELF grid is not a finite three-dimensional array.")
    species = center_species(structure)
    coordination = periodic_coordination(structure, species)
    center_indices = sorted(coordination)
    fft_spacing_A = float(
        np.mean(np.asarray(structure.lattice.abc, dtype=float) / np.asarray(grid.shape))
    )

    provisional_sites = []
    bond_results: list[BondResult] = []
    for site_number, center_index in enumerate(center_indices, start=1):
        bonds, seventh_neighbor_A = coordination[center_index]
        M_values = []
        for direction, bond in enumerate(bonds, start=1):
            distances_A, elf = sample_bond(
                grid,
                structure[center_index].frac_coords,
                bond.delta_frac,
                structure.lattice.matrix,
                samples=samples,
            )
            try:
                peak = topological_peak(distances_A, elf, fft_spacing_A)
            except DescriptorError as exc:
                raise DescriptorError(
                    f"Peak identification failed for {species} site {site_number}, "
                    f"direction {direction}: {exc}"
                ) from exc
            M = peak.r_peak_A * peak.h_peak
            M_values.append(M)
            bond_results.append(
                BondResult(
                    site_number=site_number,
                    center_index=center_index,
                    direction=direction,
                    iodine_index=bond.iodine_index,
                    image_a=bond.image[0],
                    image_b=bond.image[1],
                    image_c=bond.image[2],
                    distance_A=bond.distance_A,
                    r_peak_A=peak.r_peak_A,
                    h_peak=peak.h_peak,
                    M=M,
                    first_min_A=peak.first_min_A,
                    topological_max_A=peak.topological_max_A,
                    second_min_A=peak.second_min_A,
                    prominence=peak.prominence,
                )
            )

        M_array = np.asarray(M_values, dtype=float)
        M_mean = float(np.mean(M_array))
        if not np.isfinite(M_mean) or M_mean <= 0:
            raise DescriptorError(
                f"{species} site {site_number} has a non-positive mean r*h value."
            )
        A_j = float(np.ptp(M_array) / M_mean)
        sixth_neighbor_A = max(bond.distance_A for bond in bonds)
        provisional_sites.append(
            {
                "site_number": site_number,
                "center_index": center_index,
                "n_iodine_neighbors": len(bonds),
                "sixth_neighbor_A": sixth_neighbor_A,
                "seventh_neighbor_A": seventh_neighbor_A,
                "shell_gap_A": seventh_neighbor_A - sixth_neighbor_A,
                "M_mean": M_mean,
                "M_min": float(np.min(M_array)),
                "M_max": float(np.max(M_array)),
                "A_j": A_j,
                "A_j_squared": A_j * A_j,
            }
        )

    selected = max(provisional_sites, key=lambda row: row["A_j"])
    sites = tuple(
        SiteResult(
            **row,
            is_most_anisotropic=row["site_number"] == selected["site_number"],
        )
        for row in provisional_sites
    )
    A_ELF = float(selected["A_j"])
    A_ELF_squared = A_ELF * A_ELF
    category = reference_category(A_ELF_squared)
    reference, physical = interpretation(species, category)
    return AnalysisResult(
        input_path=input_path,
        center_species=species,
        site_count=len(sites),
        coordination_passed=True,
        most_anisotropic_site=int(selected["site_number"]),
        most_anisotropic_center_index=int(selected["center_index"]),
        A_ELF=A_ELF,
        A_ELF_squared=A_ELF_squared,
        reference_category=category,
        reference_statement=reference,
        physical_interpretation=physical,
        sites=sites,
        bonds=tuple(bond_results),
    )
