from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from .errors import DescriptorError
from .models import Peak


DEFAULT_SAMPLES = 2000
REFINEMENT_HALF_WIDTH_A = 0.13


def _cubic_1d(p0: float, p1: float, p2: float, p3: float, t: float) -> float:
    a0 = -0.5 * p0 + 1.5 * p1 - 1.5 * p2 + 0.5 * p3
    a1 = p0 - 2.5 * p1 + 2.0 * p2 - 0.5 * p3
    a2 = -0.5 * p0 + 0.5 * p2
    return float(((a0 * t + a1) * t + a2) * t + p1)


def tricubic_periodic(grid: np.ndarray, frac: Sequence[float]) -> float:
    """Evaluate a periodic volumetric grid using Catmull-Rom tricubic interpolation."""
    shape = np.asarray(grid.shape, dtype=int)
    scaled = np.asarray(frac, dtype=float) * shape
    base = np.floor(scaled).astype(int)
    offset = scaled - np.floor(scaled)
    z_values = []
    for dz in (-1, 0, 1, 2):
        y_values = []
        for dy in (-1, 0, 1, 2):
            x_values = [
                grid[
                    (base[0] + dx) % shape[0],
                    (base[1] + dy) % shape[1],
                    (base[2] + dz) % shape[2],
                ]
                for dx in (-1, 0, 1, 2)
            ]
            y_values.append(_cubic_1d(*x_values, float(offset[0])))
        z_values.append(_cubic_1d(*y_values, float(offset[1])))
    return _cubic_1d(*z_values, float(offset[2]))


def sample_bond(
    grid: np.ndarray,
    center_frac: Sequence[float],
    delta_frac: Sequence[float],
    lattice_matrix: np.ndarray,
    *,
    samples: int = DEFAULT_SAMPLES,
) -> tuple[np.ndarray, np.ndarray]:
    if samples < 20:
        raise DescriptorError("At least 20 samples per Sn/Pb-I direction are required.")
    center = np.asarray(center_frac, dtype=float)
    delta = np.asarray(delta_frac, dtype=float)
    length_A = float(np.linalg.norm(delta @ lattice_matrix))
    fraction = np.linspace(0.0, 1.0, samples)
    points = (center[None, :] + fraction[:, None] * delta[None, :]) % 1.0
    elf = np.asarray([tricubic_periodic(grid, point) for point in points], dtype=float)
    return fraction * length_A, elf


def gaussian_smooth(values: np.ndarray, dx_A: float, fwhm_A: float) -> np.ndarray:
    sigma_points = max(1.0, fwhm_A / (2.354820045 * dx_A))
    radius = int(np.ceil(4.0 * sigma_points))
    offsets = np.arange(-radius, radius + 1, dtype=float)
    kernel = np.exp(-0.5 * (offsets / sigma_points) ** 2)
    kernel /= kernel.sum()
    padded = np.pad(values, radius, mode="reflect")
    return np.convolve(padded, kernel, mode="valid")


def extrema_indices(values: np.ndarray, *, maxima: bool) -> np.ndarray:
    signed = values if maxima else -values
    return np.flatnonzero(
        (signed[1:-1] > signed[:-2]) & (signed[1:-1] >= signed[2:])
    ) + 1


def refine_quadratic(
    distances_A: np.ndarray,
    elf: np.ndarray,
    center_index: int,
    *,
    half_width_A: float = REFINEMENT_HALF_WIDTH_A,
) -> float:
    x0 = float(distances_A[center_index])
    mask = np.abs(distances_A - x0) <= half_width_A
    if int(mask.sum()) < 4:
        raise DescriptorError("Too few raw-profile points for quadratic peak refinement.")
    local_x = distances_A[mask] - x0
    coefficients = np.polyfit(local_x, elf[mask], 2)
    roots = np.roots(np.polyder(coefficients))
    roots = roots[np.isreal(roots)].real
    roots = roots[np.abs(roots) <= half_width_A]
    maxima = [
        root for root in roots
        if np.polyval(np.polyder(coefficients, 2), root) < 0
    ]
    if not maxima:
        raise DescriptorError("Quadratic refinement did not identify a local maximum.")
    vertex = max(maxima, key=lambda root: np.polyval(coefficients, root))
    return float(x0 + vertex)


def topological_peak(
    distances_A: np.ndarray,
    elf: np.ndarray,
    fft_spacing_A: float,
) -> Peak:
    """Find the first post-nuclear minimum-maximum-minimum ELF basin."""
    if distances_A.ndim != 1 or elf.ndim != 1 or distances_A.shape != elf.shape:
        raise DescriptorError("Each ELF profile must contain matching one-dimensional arrays.")
    if len(distances_A) < 20 or not np.all(np.isfinite(elf)):
        raise DescriptorError("The ELF profile is incomplete or contains non-finite values.")
    smooth = gaussian_smooth(
        elf,
        float(np.median(np.diff(distances_A))),
        fft_spacing_A,
    )
    minima = extrema_indices(smooth, maxima=False)
    maxima = extrema_indices(smooth, maxima=True)
    if not len(minima):
        raise DescriptorError("No post-nuclear minimum was found in an ELF profile.")
    first_min = int(minima[0])
    following_maxima = maxima[maxima > first_min]
    if not len(following_maxima):
        raise DescriptorError("No maximum follows the first post-nuclear minimum.")
    topological_max = int(following_maxima[0])
    following_minima = minima[minima > topological_max]
    if not len(following_minima):
        raise DescriptorError("No minimum follows the first metal-centered ELF maximum.")
    second_min = int(following_minima[0])

    basin = np.arange(first_min, second_min + 1)
    raw_max = int(basin[np.argmax(elf[basin])])
    r_peak_A = refine_quadratic(distances_A, elf, raw_max)
    h_peak = float(np.interp(r_peak_A, distances_A, elf))
    prominence = min(
        float(smooth[topological_max] - smooth[first_min]),
        float(smooth[topological_max] - smooth[second_min]),
    )
    return Peak(
        first_min_A=float(distances_A[first_min]),
        topological_max_A=float(distances_A[topological_max]),
        second_min_A=float(distances_A[second_min]),
        prominence=prominence,
        r_peak_A=r_peak_A,
        h_peak=h_peak,
    )
