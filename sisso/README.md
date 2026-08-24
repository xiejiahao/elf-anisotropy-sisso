# SISSO++ descriptor search

This directory reproduces the sparse descriptor search from the reported
13-material data set. It operates only on tables supplied locally under
`data/sisso/`; the tables are not included in the current public snapshot, and
the workflow does not read ELFCAR files or recalculate the one-dimensional ELF
profiles.

The search space contains 1,539 dimensionless candidates derived from the
bond-level quantities `d`, `r`, and `h`. One dimension (`D=1`) is the primary
search, while `D=2` tests whether a second descriptor improves held-out
predictions. Model selection uses nested leave-one-material-out validation:
each outer material is excluded from every inner search and from the final fit
used to predict that material.

## Data files

The commands below use these local repository-relative files by default:

```text
data/sisso/candidate_matrix.csv
data/sisso/candidate_ancestry.csv
data/sisso/targets.csv
data/sisso/sample_map.csv
data/sisso/reference_results/
```

`candidate_matrix.csv` contains one row per material and one column per
candidate. `candidate_ancestry.csv` records the local expression, within-site
reducer, across-site reducer, outer transform, and unit for each feature. The
preparation step verifies that both tables contain the same 1,539 feature IDs,
that all values are finite, and that all candidates are dimensionless A-module
family representatives. It does not use FWHM to filter candidates.

## 1. Build the pinned SISSO++ executable

The calculations used official SISSO++ commit
`43b99110118a51b9f4983b02c8d781ae6f25456c`. SISSO++ is not included in this
repository.

```bash
bash sisso/install_sissopp.sh /path/to/sissopp-src /path/to/sissopp-build
export SISSOPP_BINARY=/path/to/sissopp-build/bin/sisso++
```

The build requires a C++17 compiler, CMake 3.20.2 or newer, MPI, OpenMP,
BLAS/LAPACK, Boost, and the linear-programming dependencies required by
SISSO++. If these dependencies are managed by a cluster or package manager,
load them before running the script. See the upstream
[SISSO++ documentation](https://sissopp_developers.gitlab.io/sissopp/) for
platform-specific toolchain files.

## 2. Prepare and validate the candidate bank

```bash
python3 sisso/prepare_candidates.py
```

This writes the SISSO++ `data.csv`, an ordered sample map, a feature map, and a
preparation summary under `sisso/work/inputs/`. Input and output paths can be
changed with command-line options; run `--help` for details.

## 3. Run nested validation

```bash
python3 sisso/run_nested_loocv.py \
  --binary "$SISSOPP_BINARY" \
  --dimensions 1,2 \
  --workers 8
```

For each outer fold, the runner evaluates `n_sis = 5, 10, 20, 50, 100` using
inner leave-one-material-out predictions on the remaining 12 materials. It
then refits with the selected `n_sis` and predicts the untouched outer
material. `--resume` reuses completed SISSO++ model files.

The composed candidates are already present as input columns, so SISSO++ uses
`max_rung = 0`; it selects sparse one- or two-dimensional models without
generating an additional algebraic rung.

## 4. Summarize the selected formulas

```bash
python3 sisso/summarize_results.py --verify-reference
```

The summary reports D=1 and D=2 outer-fold metrics and maps every selected
feature back to its formula ancestry. `--verify-reference` compares the
recomputed metrics with the archived reference results in `data/sisso/`.

## Tests

The tests exercise table validation, fold exclusion, configuration generation,
model parsing, and metrics without compiling SISSO++:

```bash
python3 -m unittest discover -s sisso/tests -v
```

## SISSO++ references

- Ouyang, R.; Curtarolo, S.; Ahmetcik, E.; Scheffler, M.; Ghiringhelli, L. M. SISSO: A Compressed-Sensing Method for Identifying the Best Low-Dimensional Descriptor in an Immensity of Offered Candidates. *Physical Review Materials* **2018**, *2*, 083802. https://doi.org/10.1103/PhysRevMaterials.2.083802
- Purcell, T. A. R.; Scheffler, M.; Carbogno, C.; Ghiringhelli, L. M. SISSO++: A C++ Implementation of the Sure-Independence Screening and Sparsifying Operator Approach. *Journal of Open Source Software* **2022**, *7*, 3960. https://doi.org/10.21105/joss.03960
