# update_s5_filenames.ps1
# Run from the Supplementary_Figure_S5 folder.

$ErrorActionPreference = "Stop"

# Scripts to modify. It is okay if only the old or renamed versions exist.
$scriptFiles = @(
    "analyze_borda_weight_sensitivity_complete_rescue_v2.py",
    "Supplementary_Figure_S5_generate_source_data.py",
    "plot_reviewer3_sensitivity_clean_v3.py",
    "Supplementary_Figure_S5_plot_clean_panels.py"
) | Where-Object { Test-Path -LiteralPath $_ }

if ($scriptFiles.Count -eq 0) {
    throw "No S5 Python scripts found in this folder."
}

# Ordered literal replacements inside scripts.
$replacements = @(
    ,@('BASE = Path("TCGA_melanoma")',
       'BASE = Path(__file__).resolve().parent')

    ,@('FINAL_DIR = BASE / "epitopes_prioritisation_complete_rescue" / "final_epitopes"',
       'FINAL_DIR = BASE / "input_final_epitope_tables"')

    ,@('OUT_DIR = BASE / "rescue_final_analysis" / "reviewer_updates" / "reviewer_3_borda_sensitivity_v2"',
       'OUT_DIR = BASE')

    ,@('TABLE_DIR = OUT_DIR / "tables"',
       'TABLE_DIR = OUT_DIR')

    ,@('R3 = BASE / "rescue_final_analysis" / "reviewer_updates" / "reviewer_3_borda_sensitivity_v2"',
       'R3 = BASE')

    ,@('SUMMARY = R3 / "tables" / "reviewer3_rank_stability_summary.tsv"',
       'SUMMARY = R3 / "Supplementary_Figure_S5_rank_stability_summary.tsv"')

    ,@('leave.to_csv(R3 / "tables" / "reviewer3_leave_one_network_metric_summary.tsv", sep="\t", index=False)',
       'leave.to_csv(R3 / "Supplementary_Figure_S5_leave_one_network_metric_summary.tsv", sep="\t", index=False)')

    ,@('reviewer3_weighting_scenarios.tsv',
       'Supplementary_Figure_S5_weighting_scenarios.tsv')

    ,@('reviewer3_columns_used.tsv',
       'Supplementary_Figure_S5_columns_used.tsv')

    ,@('reviewer3_epitope_scores_and_sensitivity_ranks.tsv.gz',
       'Supplementary_Figure_S5_epitope_scores_and_sensitivity_ranks.tsv.gz')

    ,@('reviewer3_rank_stability_by_sample.tsv',
       'Supplementary_Figure_S5_rank_stability_by_sample.tsv')

    ,@('reviewer3_rank_stability_summary.tsv',
       'Supplementary_Figure_S5_rank_stability_summary.tsv')

    ,@('reviewer3_feature_pairwise_correlations_by_sample.tsv',
       'Supplementary_Figure_S5_feature_pairwise_correlations_by_sample.tsv')

    ,@('reviewer3_feature_pairwise_correlations_summary.tsv',
       'Supplementary_Figure_S5_feature_pairwise_correlations_summary.tsv')

    ,@('reviewer3_feature_pairwise_correlations_overall.tsv',
       'Supplementary_Figure_S5_feature_pairwise_correlations_overall.tsv')

    ,@('reviewer3_leave_one_network_metric_summary.tsv',
       'Supplementary_Figure_S5_leave_one_network_metric_summary.tsv')
)

foreach ($file in $scriptFiles) {
    $backup = "$file.bak"
    if (-not (Test-Path -LiteralPath $backup)) {
        Copy-Item -LiteralPath $file -Destination $backup
    }

    $content = Get-Content -LiteralPath $file -Raw

    foreach ($pair in $replacements) {
        $content = $content.Replace($pair[0], $pair[1])
    }

    Set-Content -LiteralPath $file -Value $content -NoNewline
    Write-Host "Updated script: $file"
}

# Rename data files if they still have the old names.
$renames = [ordered]@{
    "reviewer3_weighting_scenarios.tsv" = "Supplementary_Figure_S5_weighting_scenarios.tsv"
    "reviewer3_columns_used.tsv" = "Supplementary_Figure_S5_columns_used.tsv"
    "reviewer3_epitope_scores_and_sensitivity_ranks.tsv.gz" = "Supplementary_Figure_S5_epitope_scores_and_sensitivity_ranks.tsv.gz"
    "reviewer3_rank_stability_by_sample.tsv" = "Supplementary_Figure_S5_rank_stability_by_sample.tsv"
    "reviewer3_rank_stability_summary.tsv" = "Supplementary_Figure_S5_rank_stability_summary.tsv"
    "reviewer3_feature_pairwise_correlations_by_sample.tsv" = "Supplementary_Figure_S5_feature_pairwise_correlations_by_sample.tsv"
    "reviewer3_feature_pairwise_correlations_summary.tsv" = "Supplementary_Figure_S5_feature_pairwise_correlations_summary.tsv"
    "reviewer3_feature_pairwise_correlations_overall.tsv" = "Supplementary_Figure_S5_feature_pairwise_correlations_overall.tsv"
    "reviewer3_leave_one_network_metric_summary.tsv" = "Supplementary_Figure_S5_leave_one_network_metric_summary.tsv"
}

foreach ($old in $renames.Keys) {
    $new = $renames[$old]

    if ((Test-Path -LiteralPath $old) -and -not (Test-Path -LiteralPath $new)) {
        Rename-Item -LiteralPath $old -NewName $new
        Write-Host "Renamed file: $old -> $new"
    }
}

# Rename scripts if they still have old names.
if ((Test-Path -LiteralPath "analyze_borda_weight_sensitivity_complete_rescue_v2.py") -and
    -not (Test-Path -LiteralPath "Supplementary_Figure_S5_generate_source_data.py")) {
    Rename-Item -LiteralPath "analyze_borda_weight_sensitivity_complete_rescue_v2.py" -NewName "Supplementary_Figure_S5_generate_source_data.py"
}

if ((Test-Path -LiteralPath "plot_reviewer3_sensitivity_clean_v3.py") -and
    -not (Test-Path -LiteralPath "Supplementary_Figure_S5_plot_clean_panels.py")) {
    Rename-Item -LiteralPath "plot_reviewer3_sensitivity_clean_v3.py" -NewName "Supplementary_Figure_S5_plot_clean_panels.py"
}

Write-Host ""
Write-Host "Done. Backups were created as *.bak files."
Write-Host "Now check for leftover internal names:"
Write-Host 'Select-String -Path *.py -Pattern "reviewer3_|complete_rescue|reviewer_3_borda|TCGA_melanoma|/ `"tables`""'