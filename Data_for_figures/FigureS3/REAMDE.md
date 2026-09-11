# Supplementary Figure S3 source data

This directory contains the numerical source data underlying Supplementary Figure S3 from the manuscript:

**Epitope Generation Gateway (EGG): A biological context-dependent approach to cancer vaccine target prioritisation**

Supplementary Figure S3 evaluates whether tumour DNA variant allele fraction (VAF) differs between all somatic epitope candidates and top-ranked somatic candidates, and whether VAF is associated with final prioritisation rank.

## Files

### Source-data files

- `Supplementary_Figure_S3_source_data_somatic_candidates_with_vaf.tsv`  
  Per-candidate source-data table used to generate the plotted values in Supplementary Figure S3. This file contains somatic epitope candidates matched to tumour DNA VAF values.

- `Supplementary_Figure_S3_summary_statistics.tsv`  
  Summary statistics used for the figure annotations, including Mann–Whitney and Spearman test results.

- `Supplementary_Figure_S3_vaf_match_summary.tsv`  
  Summary of VAF matching completeness and the VAF threshold used in the analysis.

### Optional provenance script

- `Supplementary_Figure_S3_generate_source_data_and_plots.py`  
  Python script documenting how the VAF source-data tables and plot components were generated from EGG final epitope tables and pVACseq somatic epitope outputs.

## Figure panels

### Panel A

Panel A compares tumour DNA VAF between:

- all final somatic epitope candidates; and
- the somatic subset of top-10 candidates per patient.

The relevant columns in `Supplementary_Figure_S3_source_data_somatic_candidates_with_vaf.tsv` are:

- `tumor_dna_vaf`
- `is_top10_overall_per_sample`

The Mann–Whitney p-value shown in the figure is provided in:

- `Supplementary_Figure_S3_summary_statistics.tsv`

under the row:

- `analysis = all_somatic_vs_top10_overall`

### Panel B

Panel B shows tumour DNA VAF versus final prioritisation-rank percentile, stratified by MHC class.

The relevant columns in `Supplementary_Figure_S3_source_data_somatic_candidates_with_vaf.tsv` are:

- `tumor_dna_vaf`
- `rank_percentile_within_sample`
- `MHC_Class_clean`

The Spearman correlation coefficients, p-values, and sample sizes shown in the figure are provided in:

- `Supplementary_Figure_S3_summary_statistics.tsv`

under the rows:

- `analysis = MHC I_vaf_vs_overall_rank_percentile`
- `analysis = MHC II_vaf_vs_overall_rank_percentile`

## Notes on the VAF threshold

The analysis uses a tumour DNA VAF threshold of 0.25 as a simplified clonality proxy. This threshold is recorded in:

- `Supplementary_Figure_S3_vaf_match_summary.tsv`

and is also represented in the source-data table by the column:

- `clonality_proxy`

## Reproducing the analysis

The provided source-data TSV files are sufficient to reproduce the numerical content of Supplementary Figure S3.

The optional Python script can be used to regenerate the source-data tables and plot components from upstream EGG outputs. It expects the following input directory structure:

```text
TCGA_melanoma/
  epitopes_prioritisation/
    final_epitopes/
      *_epitopes_final.csv

  2A_somatic_mutation_epitopes/
    pvacSeq/
      Sample_*/
        MHC_Class_*/
          *.filtered.tsv
```

## Licence

The numerical source-data files in this directory are released under the Creative Commons CC0 1.0 Universal Public Domain Dedication.

This CC0 dedication applies to the processed numerical source-data tables in this directory. Source code in the main repository remains under the repository-level GPL-3.0 licence unless otherwise stated.
