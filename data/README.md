# Source data

This directory contains the numerical data used for the ELF-anisotropy descriptor and the associated SISSO++ analysis. The release starts from one-dimensional ELF profiles; volumetric ELFCAR files and crystal structures are not included.

The profiles were extracted from structures whose internal coordinates were relaxed while retaining the experimental lattice parameters. For every crystallographically represented Sn site, the six nearest periodic I neighbors were identified from the full periodic structure. Each Sn–I segment was sampled at 2,000 equally spaced positions using periodic tricubic interpolation of the ELF grid.

## ELF profiles and descriptor values

### `elf_profiles/<material>.csv.gz`

One gzip-compressed CSV file is provided for each of the 13 Sn iodide materials. Every file has five columns:

| Column | Definition |
| --- | --- |
| `material` | Material identifier used throughout this repository. |
| `site` | One-based Sn-site number within the structure. |
| `direction` | One-based index of one of the six periodic Sn–I coordination directions. Directions are ordered by Sn–I distance, I-site index, and periodic image. |
| `distance_A` | Distance from Sn along the Sn–I segment, in Å. |
| `ELF` | Interpolated electron localization function. |

The 13 files contain 312 profiles and 624,000 data rows in total.

### `bond_features.csv`

One row is provided for each of the 312 Sn–I directions.

| Column | Definition |
| --- | --- |
| `material`, `site`, `direction` | Keys linking the row to an ELF profile. |
| `bond_length_A` | Periodic Sn–I distance, in Å. |
| `fft_spacing_A` | Mean real-space spacing of the parent ELF grid, in Å; this sets the smoothing width used only to determine extremum ordering. |
| `r_peak_A` | Refined position of the first post-nuclear ELF maximum, in Å from Sn. |
| `h_peak` | ELF value at the refined maximum. |
| `M` | Directional peak quantity, `M = r_peak_A × h_peak`. |

The peak basin is identified as the first minimum–maximum–minimum sequence after the Sn nucleus on a profile smoothed only for extremum ordering. The maximum is then located on the unsmoothed interpolated profile and refined by a local quadratic fit. No fixed radial search interval is used to select the basin.

### `site_descriptor_values.csv`

For Sn site `j`, the six directional values are summarized as

```text
A_j = (max_i M_ij - min_i M_ij) / mean_i M_ij
```

| Column | Definition |
| --- | --- |
| `material`, `site` | Material and one-based Sn-site number. |
| `M_min`, `M_max`, `M_mean` | Minimum, maximum, and mean of the six directional `M` values. |
| `A_j` | Dimensionless ELF anisotropy of Sn site `j`. |
| `is_material_max` | `True` when the site attains the material-level maximum. |

### `material_descriptor_values.csv`

The material descriptor is

```text
A_ELF = max_j A_j
A_ELF_squared = A_ELF^2
```

| Column | Definition |
| --- | --- |
| `sample_id` | Stable SISSO++ sample identifier. |
| `material` | Material identifier. |
| `n_Sn_sites` | Number of represented Sn sites. |
| `max_site` | Semicolon-separated site number or numbers attaining `A_ELF`. |
| `A_ELF` | Maximum site-level ELF anisotropy. |
| `A_ELF_squared` | Squared material descriptor used in the reported relation. |
| `FWHM_meV` | Experimental PL full width at half maximum, in meV. |

## SISSO++ input data

The `sisso/` directory contains the complete 1,539-feature, dimensionless A-only candidate space used for the reported search.

### `sisso/candidate_matrix.csv`

Rows are materials and columns are candidate features. `sample_id` and `material` identify each row; the remaining columns are feature IDs that map to `candidate_ancestry.csv`. The target is deliberately stored separately.

### `sisso/candidate_ancestry.csv`

This table defines every candidate feature and its construction history.

| Column | Definition |
| --- | --- |
| `candidate_key` | Canonical expression key. |
| `required_modules` | Primitive-data module required by the expression; all released candidates use module A. |
| `object_type` | Local object over which primitives are defined. |
| `local_expression` | Local algebraic expression applied before aggregation. |
| `local_operation` | Operation used to construct the local expression. |
| `source_primitives` | Primitive variables entering the local expression (`d`, `r`, and/or `h`). |
| `within_reducer` | Aggregation over the six directions at one Sn site. |
| `across_reducer` | Aggregation over inequivalent Sn sites in one material. |
| `unit` | Dimensional unit; `1` denotes a dimensionless feature. |
| `outer_transform` | Final transformation applied after aggregation. |
| `numeric_family_id` | Identifier grouping numerically equivalent candidates. |
| `is_family_representative` | Whether the candidate represents its numeric family. |
| `feature_id` | Column identifier used in `candidate_matrix.csv`. |
| `grid_dependent` | Whether the feature depends on ELF-grid-derived quantities. |
| `max_grid_spread_pct_of_full_range` | Grid-sensitivity estimate expressed as a percentage of the feature's full material range. |
| `grid_stable` | Whether the feature passed the numerical-stability filter. |

### `sisso/targets.csv` and `sisso/sample_map.csv`

`targets.csv` provides `sample_id`, `material`, and experimental `FWHM_meV`. `sample_map.csv` provides the one-to-one mapping between the stable sample IDs and material identifiers.

## Reference results

`sisso/reference_results/` contains the reported nested leave-one-out results:

- `D1_metrics.json` and `D2_metrics.json`: aggregate outer-loop metrics for one- and two-dimensional models.
- `D1_outer_predictions.csv` and `D2_outer_predictions.csv`: held-out predictions, selected SIS subspace size, and fitted expression for each outer fold.
- `D1_selected_formula_families.csv`: ancestry of the one-dimensional descriptor selected in each outer fold.

These files are reference outputs for checking a fresh run of the public workflow; they are not additional inputs to model selection.
