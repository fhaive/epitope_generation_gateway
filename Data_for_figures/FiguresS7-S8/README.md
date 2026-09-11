# Supplementary Figures S7-S8 source data

This package contains the numerical source data underlying Supplementary Figures S7 and S8.

## Supplementary Figure S7

Supplementary Figure S7 shows database-stratified counts of EGG candidates supported by external evidence under peptide-coverage thresholds.

Minimum source-data file:

- `tables/Supplementary_Figure_S7_database_counts_no_any_plotted_values.tsv`
- `tables/supported_candidate_counts_by_database_and_coverage_threshold.tsv`

The plotted figure files are provided in `figures/` as:

- `Supplementary_Figure_S7_database_counts_no_any.pdf`
- `Supplementary_Figure_S7_database_counts_no_any.png`
- `Supplementary_Figure_S7_database_counts_no_any.svg`

## Supplementary Figure S8

Supplementary Figure S8 shows the effect of peptide-length-normalised coverage thresholds on ranking-based recovery of external evidence.

Minimum source-data files:

- `tables/summary_ranking_delta_vs_ic50_coverage_thresholds.tsv`
- `tables/per_patient_ranking_delta_vs_ic50_coverage_thresholds.tsv`
- `tables/ranking_methods_used.tsv`

Panel-specific plotted-value files:

- `tables/strict_positive_experimental__relaxed_min8_peptide__topn_original_logic.numbers.tsv`
- `tables/strict_positive_experimental__relaxed_min8_peptide__percent_original_logic.numbers.tsv`
- `tables/strict_positive_experimental__relaxed_min8_peptide_cov50__topn_original_logic.numbers.tsv`
- `tables/strict_positive_experimental__relaxed_min8_peptide_cov50__percent_original_logic.numbers.tsv`
- `tables/strict_positive_experimental__relaxed_min8_peptide_cov70__topn_original_logic.numbers.tsv`
- `tables/strict_positive_experimental__relaxed_min8_peptide_cov70__percent_original_logic.numbers.tsv`
- `tables/strict_positive_experimental__relaxed_min8_peptide_cov80__topn_original_logic.numbers.tsv`
- `tables/strict_positive_experimental__relaxed_min8_peptide_cov80__percent_original_logic.numbers.tsv`

The plotted figure files are provided in `figures/` as:

- `coverage_thresholds_all_methods_previous_style_original_logic.pdf`
- `coverage_thresholds_all_methods_previous_style_original_logic.png`
- `coverage_thresholds_all_methods_previous_style_original_logic.svg`

## Notes

The `.numbers.tsv` files contain the exact values used for the plotted panels. The summary and per-patient tables provide the underlying ranking comparison data. The processed candidate-level support table is also included for provenance.

## Licence

The numerical source-data files in this package are released under the Creative Commons CC0 1.0 Universal Public Domain Dedication.

This CC0 dedication applies to the processed numerical source-data tables. Source code in the main repository remains under the repository-level GPL-3.0 licence unless otherwise stated.
