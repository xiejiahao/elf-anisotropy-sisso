from __future__ import annotations

import math

import numpy as np
import pytest

from elf_anisotropy.analysis import interpretation, reference_category
from elf_anisotropy.coordination import periodic_coordination
from elf_anisotropy.errors import DescriptorError
from elf_anisotropy.profiles import topological_peak, tricubic_periodic


class _Specie:
    def __init__(self, symbol: str) -> None:
        self.symbol = symbol


class _Site:
    def __init__(self, symbol: str, frac_coords) -> None:
        self.specie = _Specie(symbol)
        self.frac_coords = np.asarray(frac_coords, dtype=float)


class _Neighbor(_Site):
    def __init__(self, index: int, frac_coords) -> None:
        super().__init__("I", frac_coords)
        self.index = index
        self.image = np.zeros(3, dtype=int)
        self.nn_distance = float(np.linalg.norm(self.frac_coords))


class _Lattice:
    matrix = np.eye(3)


class _SixCoordinateStructure:
    def __init__(self) -> None:
        vectors = (
            (0.2, 0.0, 0.0),
            (-0.2, 0.0, 0.0),
            (0.0, 0.2, 0.0),
            (0.0, -0.2, 0.0),
            (0.0, 0.0, 0.2),
            (0.0, 0.0, -0.2),
        )
        self.sites = [_Site("Sn", (0.0, 0.0, 0.0))]
        self.sites.extend(_Site("I", vector) for vector in vectors)
        self.neighbors = [
            _Neighbor(index, vector) for index, vector in enumerate(vectors, start=1)
        ]
        self.lattice = _Lattice()

    def __iter__(self):
        return iter(self.sites)

    def __getitem__(self, index: int):
        return self.sites[index]

    def get_neighbors(self, center, radius: float):
        return self.neighbors


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (0.0, "weak"),
        (0.000999, "weak"),
        (0.001, "moderate"),
        (0.029999, "moderate"),
        (0.03, "strong"),
    ],
)
def test_reference_category_boundaries(value: float, expected: str) -> None:
    assert reference_category(value) == expected


def test_interpretation_uses_species_orbital_and_reference_scope() -> None:
    reference, physical = interpretation("Pb", "moderate")
    assert "13-compound Sn DJ reference set" in reference
    assert "Pb(II) 6s2" in physical
    assert "finite-temperature" in physical


def test_tricubic_interpolation_is_periodic_for_constant_grid() -> None:
    grid = np.full((5, 6, 7), 0.42)
    assert tricubic_periodic(grid, (0.2, 0.3, 0.4)) == pytest.approx(0.42)
    assert tricubic_periodic(grid, (1.2, -0.7, 2.4)) == pytest.approx(0.42)


def test_exactly_six_periodic_iodine_directions_are_accepted() -> None:
    coordination = periodic_coordination(_SixCoordinateStructure(), "Sn")
    bonds, seventh_neighbor_A = coordination[0]
    assert len(bonds) == 6
    assert math.isnan(seventh_neighbor_A)


def test_topological_peak_finds_first_minimum_maximum_minimum_basin() -> None:
    x = np.linspace(0.0, 3.0, 2000)
    y = (
        0.4
        + 0.5 * np.exp(-((x - 0.02) / 0.08) ** 2)
        + 0.35 * np.exp(-((x - 0.72) / 0.16) ** 2)
        + 0.20 * np.exp(-((x - 2.25) / 0.25) ** 2)
    )
    peak = topological_peak(x, y, fft_spacing_A=0.08)
    assert peak.first_min_A < peak.r_peak_A < peak.second_min_A
    assert peak.r_peak_A == pytest.approx(0.72, abs=0.01)
    assert peak.h_peak > 0.7


def test_topological_peak_rejects_monotonic_profile() -> None:
    x = np.linspace(0.0, 3.0, 2000)
    with pytest.raises(DescriptorError, match="post-nuclear minimum"):
        topological_peak(x, np.linspace(1.0, 0.0, 2000), fft_spacing_A=0.08)
