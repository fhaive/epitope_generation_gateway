# Supplementary Figure S4 source data

This directory contains the numerical source data underlying Supplementary Figure S4.

Supplementary Figure S4 evaluates the effect of removing PPI/topology-derived metrics from the final EGG ranking.

## Files

- `Supplementary_Figure_S4_source_data_by_sample.tsv`  
  Per-sample values used to generate the two histogram panels.

- `Supplementary_Figure_S4_cohort_summary.tsv`  
  Cohort-level summary statistics for the ablation analysis.

- `Supplementary_Figure_S4_weighting_scheme.tsv`  
  Original and PPI/topology-ablated weighting schemes.

- `Supplementary_Figure_S4_detailed_candidate_ranks.tsv.gz`  
  Optional per-candidate rank-level comparison underlying the per-sample summaries.

- `Supplementary_Figure_S4_generate_source_data_and_plots.py`  
  Optional script documenting how the source-data tables and plots were generated.

## Figure panels

Panel A uses the column:

```text
spearman_rho_original_vs_ppi_topology_ablated_rank
```

from:

Supplementary_Figure_S4_source_data_by_sample.tsv

Panel B uses the column:

top10_overlap_n

from the same file.

Reproducing the plots

The provided source-data table is sufficient to reproduce the numerical content of Supplementary Figure S4.

The optional script expects final EGG epitope tables as input and regenerates the ablation tables and plot components.
