# Supplementary Figure S5 source data

This directory contains the numerical source data underlying Supplementary Figure S5.

Supplementary Figure S5 evaluates the stability of EGG/Borda prioritisation under alternative weighting schemes, leave-one-network-metric-out analyses, and changes in total network-feature weight.

## Files

- `Supplementary_Figure_S5_rank_stability_summary.tsv`  
  Scenario-level summary statistics used for panels A, B, C, and E.

- `Supplementary_Figure_S5_rank_stability_by_sample.tsv`  
  Per-sample rank-stability values used for the panel D boxplots.

- `Supplementary_Figure_S5_leave_one_network_metric_summary.tsv`  
  Focused summary table for leave-one-network-metric-out analyses.

- `Supplementary_Figure_S5_weighting_scenarios.tsv`  
  Weighting scheme used for each sensitivity scenario.

- `Supplementary_Figure_S5_epitope_scores_and_sensitivity_ranks.tsv.gz`  
  Per-candidate ranks and scores across sensitivity scenarios.

- `Supplementary_Figure_S5_columns_used.tsv`  
  Input-column provenance table.

- `Supplementary_Figure_S5_generate_source_data.py`  
  Optional script used to generate the source-data tables from final EGG epitope tables.

- `Supplementary_Figure_S5_plot_clean_panels.py`  
  Optional script used to generate the cleaned plot components.

## Figure panels

- **Panel A:** leave-one-network-metric-out full-rank similarity.  
  Uses `median_spearman_rho_vs_default`, `q1_spearman_rho_vs_default`, and `q3_spearman_rho_vs_default`.

- **Panel B:** leave-one-network-metric-out top-10 and top-20 overlap.  
  Uses `median_top10_overlap_fraction` and `median_top20_overlap_fraction`.

- **Panel C:** full-rank similarity across feature/block sensitivity scenarios.  
  Uses `median_spearman_rho_vs_default`.

- **Panel D:** top-10 overlap distributions across alternative weighting scenarios.  
  Uses `top10_overlap_with_default_fraction` from the per-sample table.

- **Panel E:** effect of changing total network weight.  
  Uses network-weight scenarios from `Supplementary_Figure_S5_rank_stability_summary.tsv`.

## Optional provenance files

The feature-correlation tables are included as supporting provenance for the sensitivity analysis:

- `Supplementary_Figure_S5_feature_pairwise_correlations_by_sample.tsv`
- `Supplementary_Figure_S5_feature_pairwise_correlations_summary.tsv`
- `Supplementary_Figure_S5_feature_pairwise_correlations_overall.tsv`

## Reproducing the plots

The provided TSV files are sufficient to reproduce the numerical content of Supplementary Figure S5.

To regenerate the cleaned plot components from the provided source-data files, run:

```bash
python Supplementary_Figure_S5_plot_clean_panels.py
```

The source-data generation script expects final EGG epitope tables in:

input_final_epitope_tables/

with filenames matching:

*_epitopes_final.csv
Software requirements
numpy
pandas
matplotlib
scipy