# ELF anisotropy descriptor

This package calculates a static electron-localization-function (ELF)
anisotropy descriptor for six-coordinate Sn(II) and Pb(II) iodide structures.
It reads a VASP `ELFCAR`, identifies the six periodic Sn/Pb-I coordination
directions at every metal site, and reports the structure-level descriptor
`A_ELF` and its squared form `A_ELF^2`.

## Installation

Python 3.11 or newer is required.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install ./descriptor
```

For a development installation with tests:

```bash
python -m pip install -e './descriptor[test]'
pytest descriptor/tests
```

## Command-line use

Analyze one file:

```bash
elf-anisotropy analyze /path/to/ELFCAR
```

Write machine-readable results:

```bash
elf-anisotropy analyze /path/to/ELFCAR --output result
```

The output directory contains:

- `summary.json`: structure-level descriptor and interpretation;
- `per_site.csv`: `A_j` for every inequivalent metal site;
- `per_bond.csv`: periodic neighbor image, bond length, peak coordinates, and
  `M_ij` for every metal-I direction.

JSON can also be printed directly:

```bash
elf-anisotropy analyze /path/to/ELFCAR --format json
```

## Descriptor definition

For direction `i` around metal site `j`, the first metal-centered ELF maximum
provides its distance from the metal, `r_ij`, and its height, `h_ij`:

```text
M_ij = r_ij h_ij
A_j = [max_i(M_ij) - min_i(M_ij)] / mean_i(M_ij)
A_ELF = max_j(A_j)
```

`A_ELF = 0` is the isotropic limit. Increasing `A_ELF` indicates stronger
directional variation of the metal-centered ELF. The command also reports
`A_ELF^2`, the form used for comparison with the experimental material series.

The first peak is defined without an absolute radial search interval. A
smoothed profile identifies the first post-nuclear minimum-maximum-minimum
basin. Peak position is then refined on the unsmoothed profile by a quadratic
fit within 0.13 Angstrom of the basin maximum, and peak height is interpolated
from the unsmoothed profile.

## Reference-scale interpretation

The text output assigns an approximate category using the descriptor
distribution of the 13-compound Sn Dion-Jacobson reference set:

| `A_ELF^2` | Reference category |
|---:|:---|
| `< 0.001` | weak |
| `0.001` to `< 0.03` | moderate |
| `>= 0.03` | strong |

These categories describe static ELF anisotropy relative to that reference
set; they are not universal phase boundaries. For Sn and Pb inputs, the text
output relates the result to the directional expression of the Sn(II) `5s2`
or Pb(II) `6s2` lone pair, respectively. The descriptor is calculated from a
static structure and does not represent the full finite-temperature dynamics
of the lone pair.

## Input checks

The calculation stops without returning a descriptor when:

- neither Sn nor Pb is present;
- iodine is absent;
- six distinct periodic iodine directions cannot be identified at every metal
  site; or
- the first post-nuclear ELF peak cannot be identified on any direction.

The periodic iodine image selected for coordination is also used for ELF-line
sampling, including in non-orthogonal cells. Sixth- and seventh-neighbor
distances and their separation are reported as diagnostics when a seventh
periodic iodine is available; no empirical distance-gap threshold is used to
accept or reject the octahedral six-coordinate environment.
