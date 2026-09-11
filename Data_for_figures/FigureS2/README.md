# Supplementary Figure S2 source data

This directory contains the numerical source data underlying Supplementary Figure S2 from the manuscript:

**Epitope Generation Gateway (EGG): A biological context-dependent approach to cancer vaccine target prioritisation**

Supplementary Figure S2 evaluates whether peptide length is associated with final EGG prioritisation. The figure contains two plotted components:

- **Panel A:** peptide-length distributions for all candidates versus top-10 candidates per patient, stratified by MHC class.
- **Panel B:** peptide length versus final prioritisation-rank percentile, stratified by MHC class.

## Files

### Source-data files

- `Supplementary_Figure_S2_source_data_per_candidate.tsv`  
  Per-candidate source-data table used to generate the plotted values in Supplementary Figure S2. This file contains peptide lengths, MHC class annotations, Borda ranks, top-10 status, and rank-percentile values.

- `Supplementary_Figure_S2_summary_statistics.tsv`  
  Summary statistics used for the figure annotations, including Mann–Whitney and Spearman test results.

- `Supplementary_Figure_S2_input_columns_used.tsv`  
  Provenance table recording which columns from the original final epitope tables were used for peptide sequence, MHC class, and rank.

### Optional scripts

- `Supplementary_Figure_S2_generate_source_data.py`  
  Python script documenting how the per-candidate source-data table and summary-statistics table were generated from the final EGG epitope tables.

- `Supplementary_Figure_S2_plot_panel_A_peptide_length_distribution.py`  
  Python script used to generate the Panel A peptide-length distribution plot.

- `Supplementary_Figure_S2_plot_panel_B_length_vs_rank_percentile.py`  
  Python script used to generate the Panel B peptide-length versus rank-percentile plot.

## Figure panels

### Panel A

Panel A compares peptide-length distributions between:

- all candidates; and
- top-10 candidates per patient.

The relevant columns in `Supplementary_Figure_S2_source_data_per_candidate.tsv` are:

- `peptide_length`
- `MHC_Class_clean`
- `is_top10_per_sample`

The Mann–Whitney p-values shown in the figure are provided in:

- `Supplementary_Figure_S2_summary_statistics.tsv`

under the columns:

- `mannwhitney_top10_vs_non_top10_p`

### Panel B

Panel B shows peptide length versus final prioritisation-rank percentile within each MHC class.

The relevant columns in `Supplementary_Figure_S2_source_data_per_candidate.tsv` are:

- `peptide_length`
- `MHC_Class_clean`
- `rank_percentile_within_sample_class`

The Spearman correlation coefficients, p-values, and sample sizes shown in the figure are provided in:

- `Supplementary_Figure_S2_summary_statistics.tsv`

under the columns:

- `spearman_n_length_vs_class_rank_percentile`
- `spearman_rho_length_vs_class_rank_percentile`
- `spearman_p_length_vs_class_rank_percentile`

## Reproducing the plotted components

The provided source-data TSV files are sufficient to reproduce the numerical content of Supplementary Figure S2.

To regenerate the plotted components from the provided source-data files, run:

```bash
python Supplementary_Figure_S2_plot_panel_A_peptide_length_distribution.py
python Supplementary_Figure_S2_plot_panel_B_length_vs_rank_percentile.py

```

These scripts read the source-data and summary-statistics files from the current directory and write output plots to:

plots/

The final multi-panel Supplementary Figure S2 was assembled from the two plotted components.

Reproducing the source-data tables

The optional source-data generation script expects final EGG epitope tables to be placed in:

input_final_epitope_tables/

with filenames matching:

*_epitopes_final.csv

The source-data files provided in this directory are already derived from these upstream outputs, so the upstream final epitope tables are not required to inspect or reuse the numerical values underlying Supplementary Figure S2.

Software requirements

The scripts require Python 3 with the following packages:

numpy
pandas
matplotlib
scipy


## Licence

The numerical source-data files in this directory are released under the Creative Commons CC0 1.0 Universal Public Domain Dedication.

This CC0 dedication applies to the processed numerical source-data tables in this directory. Source code in the main repository remains under the repository-level GPL-3.0 licence unless otherwise stated.

