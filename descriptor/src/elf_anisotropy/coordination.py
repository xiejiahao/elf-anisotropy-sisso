from __future__ import annotations

import math

import numpy as np

from .errors import DescriptorError
from .models import Bond


TARGET_COORDINATION = 6


def element_symbol(site) -> str:
    return str(site.specie.symbol)


def center_species(structure) -> str:
    species = {element_symbol(site) for site in structure}
    centers = [symbol for symbol in ("Sn", "Pb") if symbol in species]
    if not centers:
        raise DescriptorError("No Sn or Pb center was found in the ELFCAR structure.")
    if len(centers) > 1:
        raise DescriptorError("Both Sn and Pb are present; a single center species is required.")
    if "I" not in species:
        raise DescriptorError("No iodine sites were found in the ELFCAR structure.")
    return centers[0]


def periodic_coordination(
    structure,
    center_element: str,
    *,
    initial_cutoff_A: float = 3.5,
    max_cutoff_A: float = 8.0,
) -> dict[int, tuple[list[Bond], float]]:
    """Select the six nearest periodic iodine images around each metal center."""
    center_indices = [
        index for index, site in enumerate(structure)
        if element_symbol(site) == center_element
    ]
    if not center_indices:
        raise DescriptorError(f"No {center_element} sites were found.")

    result: dict[int, tuple[list[Bond], float]] = {}
    for center_index in center_indices:
        center = structure[center_index]
        radius = initial_cutoff_A
        candidates = []
        while True:
            candidates = [
                neighbor for neighbor in structure.get_neighbors(center, radius)
                if element_symbol(neighbor) == "I"
            ]
            if len(candidates) >= TARGET_COORDINATION + 1 or radius >= max_cutoff_A:
                break
            radius = min(max_cutoff_A, radius * 1.35)

        if len(candidates) < TARGET_COORDINATION:
            raise DescriptorError(
                f"{center_element} site {center_index} has fewer than six periodic "
                f"iodine directions within {max_cutoff_A:.2f} A."
            )

        unique = {}
        for neighbor in candidates:
            image = tuple(int(value) for value in np.asarray(neighbor.image, dtype=int))
            unique[(int(neighbor.index), image)] = neighbor
        ordered = sorted(
            unique.values(),
            key=lambda neighbor: (
                float(neighbor.nn_distance),
                int(neighbor.index),
                tuple(int(value) for value in neighbor.image),
            ),
        )
        if len(ordered) < TARGET_COORDINATION:
            raise DescriptorError(
                f"{center_element} site {center_index} has fewer than six distinct "
                "periodic iodine directions."
            )

        selected = ordered[:TARGET_COORDINATION]
        seventh_neighbor_A = (
            float(ordered[TARGET_COORDINATION].nn_distance)
            if len(ordered) > TARGET_COORDINATION
            else math.nan
        )

        bonds: list[Bond] = []
        for neighbor in selected:
            image = tuple(int(value) for value in np.asarray(neighbor.image, dtype=int))
            delta_frac = (
                np.asarray(structure[int(neighbor.index)].frac_coords, dtype=float)
                + np.asarray(image, dtype=float)
                - np.asarray(center.frac_coords, dtype=float)
            )
            distance_A = float(np.linalg.norm(delta_frac @ structure.lattice.matrix))
            if not math.isclose(
                distance_A,
                float(neighbor.nn_distance),
                rel_tol=0.0,
                abs_tol=1e-7,
            ):
                raise DescriptorError(
                    f"Periodic image mismatch for {center_element} site {center_index}."
                )
            bonds.append(
                Bond(
                    center_index=center_index,
                    iodine_index=int(neighbor.index),
                    image=image,
                    delta_frac=tuple(float(value) for value in delta_frac),
                    distance_A=distance_A,
                )
            )
        result[center_index] = (bonds, seventh_neighbor_A)
    return result
