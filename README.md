# ELF Anisotropy Descriptor and SISSO++ Workflow

This repository contains the code used to identify and evaluate an electron-localization-function (ELF) anisotropy descriptor for two-dimensional Sn(II) and Pb(II) iodide perovskites.

The repository has two independent entry points:

- `descriptor/`: an installable command-line program that calculates the final ELF anisotropy descriptor directly from a VASP `ELFCAR` file.
- `sisso/`: the physically constrained SISSO++ candidate-search and nested leave-one-material-out validation workflow.

## Descriptor definition

For Sn-I or Pb-I direction `i` at metal site `j`,

```text
M_ij = r_ij h_ij
A_j = [max_i(M_ij) - min_i(M_ij)] / mean_i(M_ij)
A_ELF = max_j(A_j)
```

where `r_ij` and `h_ij` are the position and height of the first metal-centered ELF maximum. The squared quantity `A_ELF^2` is the final form identified by the symbolic-regression analysis.

`A_ELF = 0` is the isotropic limit. Increasing values indicate stronger directional asymmetry of the metal-centered ELF in the analyzed static structure.

## Quick start

Install the descriptor program:

```bash
python -m venv .venv
source .venv/bin/activate
pip install ./descriptor
```

Analyze one `ELFCAR` file:

```bash
elf-anisotropy analyze /path/to/ELFCAR
```

Write machine-readable and detailed outputs:

```bash
elf-anisotropy analyze /path/to/ELFCAR --output results/example
```

See [`descriptor/README.md`](descriptor/README.md) for the algorithm, complete command reference, output schema, and interpretation guidance.

## Reproducing the descriptor search

The SISSO++ workflow expects its input tables under a local `data/` directory. These tables are not included in the current public snapshot, and the workflow does not require the original `ELFCAR` files.

```bash
python3 sisso/prepare_candidates.py
python3 sisso/run_nested_loocv.py --binary /path/to/sisso++ --dimensions 1,2
python3 sisso/summarize_results.py --verify-reference
```

See [`sisso/README.md`](sisso/README.md) for the pinned SISSO++ version, build instructions, nested-validation design, and expected results.

## Scientific scope

The descriptor measures static directional ELF anisotropy. Its strong relation to photoluminescence linewidth across the studied material series shows that it captures a dominant cross-material tendency, despite finite-temperature lone-pair dynamics not being explicitly included.

## License and citation

The code in this repository is released under the MIT License. SISSO++ is an external project and retains its own license. Citation metadata are provided in [`CITATION.cff`](CITATION.cff).
