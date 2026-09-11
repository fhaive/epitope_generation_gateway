EGG reviewer round 2 coverage-sensitivity analysis package

Purpose
-------
This package contains the scripts, summary tables, diagnostics, and figures used to assess whether relaxed min-8 external-evidence matching was influenced by peptide length.

Analysis summary
----------------
The analysis compared the original relaxed min-8 containment endpoint with stricter peptide-length-normalised coverage thresholds. For each relaxed containment match between an EGG candidate peptide and an experimentally supported external reference peptide, the matched fraction of the longer peptide was calculated. Candidate-level support was then recalculated using:
  1. Min-8 containment
  2. Min-8 + >=50% longer-peptide coverage
  3. Min-8 + >=70% longer-peptide coverage
  4. Min-8 + >=80% longer-peptide coverage

Key validation
--------------
The rebuilt min-8 endpoint from the full match-level table exactly reproduced the original candidate-level evidence flags:
  n_min8_rebuilt_vs_original_mismatches = 0

Key output
----------
Main follow-up figure:
  figures/Supplementary_Figure_S7_coverage_sensitivity_combined.pdf
  figures/Supplementary_Figure_S7_coverage_sensitivity_combined.png
  figures/Supplementary_Figure_S7_coverage_sensitivity_combined.svg

Main summary table:
  tables/summary_hit_rate_difference_vs_ic50_coverage_threshold.tsv

Candidate support counts:
  tables/supported_candidate_counts_by_coverage_threshold.tsv

Scripts
-------
scripts/00_patch_builder_for_all_matches.py
  Creates a patched copy of the original evidence-flag builder script without modifying the original script.

scripts/01_build_exact_relaxed_evidence_flags_complete_rescue_vaf_ALL_MATCHES.py
  Rebuilds exact/relaxed evidence flags and additionally writes a full match-level table with peptide coverage values.

scripts/07_percent_coverage_sensitivity_relaxed_matches.py
  Performs the coverage-threshold sensitivity analysis and generates summary tables and figures.

Original scripts are included for provenance:
  scripts/original_01_build_exact_relaxed_evidence_flags_complete_rescue_vaf.py
  scripts/original_02_compare_rankings_vs_ic50_exact_relaxed_complete_rescue_vaf.py
