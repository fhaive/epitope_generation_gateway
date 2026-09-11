#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

BASE = Path.cwd()

IN = BASE / "reviewer2_coverage_sensitivity" / "database_stratified_plots"
OUT = BASE / "reviewer2_coverage_sensitivity" / "final_s7_counts_and_all_methods"
OUT.mkdir(parents=True, exist_ok=True)

SUMMARY_IN = IN / "database_stratified_coverage_all_methods_summary.tsv"
COUNTS_IN = IN / "database_stratified_supported_candidate_counts.tsv"

CRITERIA = ["min8", "min8_cov50", "min8_cov70", "min8_cov80"]

CRITERION_LABELS = {
    "min8": "Min-8 containment",
    "min8_cov50": "Min-8 + ≥50%",
    "min8_cov70": "Min-8 + ≥70%",
    "min8_cov80": "Min-8 + ≥80%",
}

DATABASE_ORDER = ["Any", "CEDAR", "IEDB", "TSNAdb_v2"]

DATABASE_LABELS = {
    "Any": "Any database",
    "CEDAR": "CEDAR",
    "IEDB": "IEDB",
    "TSNAdb_v2": "TSNAdb v2.0",
}

METHOD_ORDER = [
    "Borda",
    "DepMap survivability",
    "Network betweenness",
    "Network degree",
    "Network impact",
    "Network strength",
    "Network WCI",
]


def plot_counts_panel(ax, counts):
    counts = counts.copy()
    counts["criterion"] = pd.Categorical(counts["criterion"], categories=CRITERIA, ordered=True)
    counts["database"] = pd.Categorical(counts["database"], categories=DATABASE_ORDER, ordered=True)
    counts = counts.sort_values(["criterion", "database"])

    x = np.arange(len(CRITERIA))
    width = 0.20

    for i, database in enumerate(DATABASE_ORDER):
        sub = counts[counts["database"].eq(database)].copy()
        values = [
            int(sub.loc[sub["criterion"].eq(c), "supported_candidates"].iloc[0])
            for c in CRITERIA
        ]

        offset = (i - 1.5) * width
        bars = ax.bar(
            x + offset,
            values,
            width,
            label=DATABASE_LABELS[database],
        )

        for bar, value in zip(bars, values):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height(),
                str(value),
                ha="center",
                va="bottom",
                fontsize=7,
            )

    ax.set_xticks(x)
    ax.set_xticklabels([CRITERION_LABELS[c] for c in CRITERIA], rotation=15, ha="right")
    ax.set_ylabel("Supported EGG candidates")
    ax.set_xlabel("Matching criterion")
    ax.set_title("A. Candidate support by database and coverage threshold")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, fontsize=8, ncol=4)


def plot_method_panel(ax, summary, criterion, cutoff_type, panel_label):
    d = summary[
        summary["database"].eq("Any")
        & summary["criterion"].eq(criterion)
        & summary["cutoff_type"].eq(cutoff_type)
    ].copy()

    for method in METHOD_ORDER:
        sub = d[d["method"].eq(method)].sort_values("cutoff")
        if sub.empty:
            continue

        ax.plot(
            sub["cutoff"],
            sub["median_difference_percentage_points"],
            marker="o",
            linewidth=1.5,
            markersize=4,
            label=method,
        )

    ax.axhline(0, linestyle="--", linewidth=1)
    ax.grid(True, alpha=0.25)

    ax.set_title(f"{panel_label}. {CRITERION_LABELS[criterion]}; {cutoff_type}")
    ax.set_xlabel(
        "Top-N candidates within each patient"
        if cutoff_type == "Top-N"
        else "Top-percent cutoff within each patient"
    )
    ax.set_ylabel("Median delta hit rate vs IC50\npercentage points")


def main():
    if not SUMMARY_IN.exists():
        raise FileNotFoundError(f"Missing summary file: {SUMMARY_IN}")
    if not COUNTS_IN.exists():
        raise FileNotFoundError(f"Missing counts file: {COUNTS_IN}")

    summary = pd.read_csv(SUMMARY_IN, sep="\t")
    counts = pd.read_csv(COUNTS_IN, sep="\t")

    # Save the source tables used for this final S7.
    summary.to_csv(OUT / "source_database_stratified_coverage_all_methods_summary.tsv", sep="\t", index=False)
    counts.to_csv(OUT / "source_database_stratified_supported_candidate_counts.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(15, 22))
    gs = GridSpec(
        nrows=5,
        ncols=2,
        height_ratios=[1.25, 1, 1, 1, 1],
        hspace=0.55,
        wspace=0.28,
        figure=fig,
    )

    ax_counts = fig.add_subplot(gs[0, :])
    plot_counts_panel(ax_counts, counts)

    panel_letters = [
        ("B", "C"),
        ("D", "E"),
        ("F", "G"),
        ("H", "I"),
    ]

    for row_idx, criterion in enumerate(CRITERIA, start=1):
        left_letter, right_letter = panel_letters[row_idx - 1]

        ax_left = fig.add_subplot(gs[row_idx, 0])
        ax_right = fig.add_subplot(gs[row_idx, 1])

        plot_method_panel(
            ax_left,
            summary,
            criterion=criterion,
            cutoff_type="Top-N",
            panel_label=left_letter,
        )

        plot_method_panel(
            ax_right,
            summary,
            criterion=criterion,
            cutoff_type="Top-percentile",
            panel_label=right_letter,
        )

    # One shared legend for ranking methods, using first method panel.
    handles, labels = fig.axes[1].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.985),
        ncol=4,
        frameon=False,
        fontsize=9,
    )

    fig.suptitle(
        "Sensitivity of external evidence support and ranking recovery to peptide-length-normalised coverage thresholds",
        y=0.998,
        fontsize=14,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.965])

    fig.savefig(OUT / "Supplementary_Figure_S7_counts_and_all_methods_coverage_effect.pdf")
    fig.savefig(OUT / "Supplementary_Figure_S7_counts_and_all_methods_coverage_effect.png", dpi=300)
    fig.savefig(OUT / "Supplementary_Figure_S7_counts_and_all_methods_coverage_effect.svg")
    plt.close(fig)

    # Also make a more compact all-methods figure without the count panel.
    fig, axes = plt.subplots(
        nrows=4,
        ncols=2,
        figsize=(14, 16),
        sharey=False,
    )

    for row_idx, criterion in enumerate(CRITERIA):
        left_letter, right_letter = panel_letters[row_idx]

        plot_method_panel(
            axes[row_idx, 0],
            summary,
            criterion=criterion,
            cutoff_type="Top-N",
            panel_label=left_letter,
        )

        plot_method_panel(
            axes[row_idx, 1],
            summary,
            criterion=criterion,
            cutoff_type="Top-percentile",
            panel_label=right_letter,
        )

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.985),
        ncol=4,
        frameon=False,
        fontsize=9,
    )

    fig.suptitle(
        "Ranking-method recovery versus IC50-only across coverage thresholds; all databases pooled",
        y=0.998,
        fontsize=14,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.965])

    fig.savefig(OUT / "all_methods_coverage_thresholds_previous_style.pdf")
    fig.savefig(OUT / "all_methods_coverage_thresholds_previous_style.png", dpi=300)
    fig.savefig(OUT / "all_methods_coverage_thresholds_previous_style.svg")
    plt.close(fig)

    print("Done.")
    print(f"Output directory: {OUT.resolve()}")
    print()
    print("Main figure:")
    print(OUT / "Supplementary_Figure_S7_counts_and_all_methods_coverage_effect.png")
    print()
    print("Compact all-methods figure:")
    print(OUT / "all_methods_coverage_thresholds_previous_style.png")


if __name__ == "__main__":
    main()
