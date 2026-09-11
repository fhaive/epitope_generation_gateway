# Supplementary Figure S6 source data

This directory contains the numerical source data underlying Supplementary Figure S6.

Supplementary Figure S6 compares recovery of positive experimental epitope evidence by EGG/Borda and individual ranking features relative to IC50-only ranking. The plotted endpoint uses positive experimental evidence with relaxed min-8 peptide matching.

## Files needed to reproduce the plotted figure

The minimum numerical source-data files for the plotted Supplementary Figure S6 panels are:

- `plot_numbers/strict_positive_experimental__relaxed_min8_peptide__topn.numbers.tsv`  
  Source data for panel A.

- `plot_numbers/strict_positive_experimental__relaxed_min8_peptide__percent.numbers.tsv`  
  Source data for panel B.

These files contain the exact values plotted in the figure, including:

- `match_scope`
- `top_type`
- `cutoff`
- `ranking_method`
- `ranking_label`
- `n_patients`
- `median_delta_hit_rate_vs_ic50_pct_points`
- `patients_better_than_ic50`
- `patients_same_as_ic50`
- `patients_worse_than_ic50`

## Figure panels

- **Panel A:** positive experimental evidence, relaxed min-8 peptide match, top-N candidate cutoffs.  
  Uses `plot_numbers/strict_positive_experimental__relaxed_min8_peptide__topn.numbers.tsv`.

- **Panel B:** positive experimental evidence, relaxed min-8 peptide match, top-percent candidate cutoffs.  
  Uses `plot_numbers/strict_positive_experimental__relaxed_min8_peptide__percent.numbers.tsv`.

## Additional source-data and provenance files

- `plot_numbers/numbers_used_for_selected_strict_positive_plots.tsv`  
  Combined table containing the plotted numerical values for selected strict-positive evidence plots.

- `plot_numbers/strict_positive_database_match_counts_used_for_reply.tsv`  
  Database-level match counts used to summarise external evidence support.

- `tables/summary_ranking_delta_vs_ic50.tsv`  
  Full ranking-delta summary table from which the selected plotted values were derived.

- `tables/per_patient_ranking_delta_vs_ic50.tsv`  
  Per-patient ranking-delta results.

- `tables/strict_positive_exact_vs_relaxed_by_database.tsv`  
  Summary of exact versus relaxed min-8 matching by evidence database.

- `current_analysis_scripts/`  
  Scripts used to generate evidence flags, ranking comparisons, summary tables, and selected plots.

## Reproducing the plots

The numerical values required to reproduce the published Supplementary Figure S6 are provided directly in the two panel-specific `.numbers.tsv` files listed above.

The analysis scripts are included for provenance. Full regeneration from the original external evidence databases requires the normalised reference database inputs used in the analysis. These large upstream database files are not included here; instead, this directory provides the processed numerical source data underlying the plotted figure.

## Licence

The numerical source-data files in this directory are released under the Creative Commons CC0 1.0 Universal Public Domain Dedication.

This CC0 dedication applies to the processed numerical source-data tables in this directory. Source code in the main repository remains under the repository-level GPL-3.0 licence unless otherwise stated.
