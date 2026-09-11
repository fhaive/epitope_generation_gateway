#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path.cwd()

IN = BASE / "reviewer2_coverage_sensitivity" / "database_stratified_plots"
OUT = BASE / "reviewer2_coverage_sensitivity" / "requested_final_plots"
OUT.mkdir(parents=True, exist_ok=True)

SUMMARY_IN = IN / "database_stratified_coverage_all_methods_summary.tsv"
COUNTS_IN = IN / "database_stratified_supported_candidate_counts.tsv"

CRITERIA = ["min8", "min8_cov50", "min8_cov70", "min8_cov80"]

CRITERION_LABELS = {
    "min8": "Min-8 containment",
    "min8_cov50": "Min-8 + ≥50% coverage",
    "min8_cov70": "Min-8 + ≥70% coverage",
    "min8_cov80": "Min-8 + ≥80% coverage",
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

DATABASE_ORDER = ["Any", "CEDAR", "IEDB", "TSNAdb_v2"]

DATABASE_LABELS = {
    "Any": "Any database",
    "CEDAR": "CEDAR",
    "IEDB": "IEDB",
    "TSNAdb_v2": "TSNAdb v2.0",
}


def plot_all_methods_any_database_by_coverage(summary):
    """
    Previous-style figure:
    - all databases pooled only
    - rows = matching threshold
    - columns = top-N and top-percentile
    - lines = ranking methods
    """

    d = summary[summary["database"].eq("Any")].copy()

    fig, axes = plt.subplots(
        nrows=len(CRITERIA),
        ncols=2,
        figsize=(14, 4.0 * len(CRITERIA)),
        sharey=False,
    )

    for row_idx, criterion in enumerate(CRITERIA):
        for col_idx, cutoff_type in enumerate(["Top-N", "Top-percentile"]):
            ax = axes[row_idx, col_idx]

            sub = d[
                d["criterion"].eq(criterion)
                & d["cutoff_type"].eq(cutoff_type)
            ].copy()

            for method in METHOD_ORDER:
                m = sub[sub["method"].eq(method)].sort_values("cutoff")
                if m.empty:
                    continue

                ax.plot(
                    m["cutoff"],
                    m["median_difference_percentage_points"],
                    marker="o",
                    label=method,
                )

            ax.axhline(0, linestyle="--", linewidth=1)
            ax.grid(True, alpha=0.25)

            ax.set_title(f"{CRITERION_LABELS[criterion]}; {cutoff_type}")
            ax.set_xlabel(
                "Top-N candidates within each patient"
                if cutoff_type == "Top-N"
                else "Top-percent cutoff within each patient"
            )
            ax.set_ylabel("Median delta hit rate vs IC50 ranking\npercentage points")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False)

    fig.suptitle(
        "Positive experimental evidence recovery across coverage thresholds; all databases pooled",
        y=0.995,
        fontsize=14,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.965])

    for ext in ["pdf", "png", "svg"]:
        fig.savefig(
            OUT / f"coverage_thresholds_all_methods_any_database_previous_style.{ext}",
            dpi=300 if ext == "png" else None,
        )

    plt.close(fig)


def plot_supported_candidate_counts(counts):
    """
    Simple grouped bar plot:
    - x-axis = matching threshold
    - bars = database
    - y-axis = number of supported EGG candidates
    """

    d = counts.copy()
    d["criterion"] = pd.Categorical(d["criterion"], categories=CRITERIA, ordered=True)
    d["database"] = pd.Categorical(d["database"], categories=DATABASE_ORDER, ordered=True)
    d = d.sort_values(["criterion", "database"])

    x = np.arange(len(CRITERIA))
    width = 0.20

    fig, ax = plt.subplots(figsize=(9, 5))

    for i, database in enumerate(DATABASE_ORDER):
        sub = d[d["database"].eq(database)].copy()
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
                fontsize=8,
                rotation=0,
            )

    ax.set_xticks(x)
    ax.set_xticklabels([CRITERION_LABELS[c] for c in CRITERIA], rotation=20, ha="right")
    ax.set_ylabel("Supported EGG candidates")
    ax.set_xlabel("External-evidence matching criterion")
    ax.set_title("Number of supported EGG candidates retained under coverage thresholds")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)

    fig.tight_layout()

    for ext in ["pdf", "png", "svg"]:
        fig.savefig(
            OUT / f"supported_candidate_counts_by_database_and_coverage_threshold.{ext}",
            dpi=300 if ext == "png" else None,
        )

    plt.close(fig)


def main():
    if not SUMMARY_IN.exists():
        raise FileNotFoundError(f"Missing summary file: {SUMMARY_IN}")
    if not COUNTS_IN.exists():
        raise FileNotFoundError(f"Missing counts file: {COUNTS_IN}")

    summary = pd.read_csv(SUMMARY_IN, sep="\t")
    counts = pd.read_csv(COUNTS_IN, sep="\t")

    plot_all_methods_any_database_by_coverage(summary)
    plot_supported_candidate_counts(counts)

    print("Done.")
    print(f"Output directory: {OUT.resolve()}")
    print()
    print("Created:")
    print(OUT / "coverage_thresholds_all_methods_any_database_previous_style.png")
    print(OUT / "supported_candidate_counts_by_database_and_coverage_threshold.png")


if __name__ == "__main__":
    main()
