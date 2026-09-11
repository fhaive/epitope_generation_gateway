#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(".").resolve()

IN_COUNTS = (
    ROOT
    / "reviewer2_coverage_sensitivity"
    / "original_logic_coverage_threshold_plots"
    / "supported_candidate_counts_by_database_and_coverage_threshold.tsv"
)

OUTDIR = (
    ROOT
    / "reviewer2_coverage_sensitivity"
    / "original_logic_coverage_threshold_plots"
    / "figures"
)
OUTDIR.mkdir(parents=True, exist_ok=True)

DATABASE_ORDER = ["CEDAR", "IEDB", "TSNAdb_v2"]

DATABASE_LABELS = {
    "CEDAR": "CEDAR",
    "IEDB": "IEDB",
    "TSNAdb_v2": "TSNAdb v2.0",
}

THRESHOLD_ORDER = [
    "relaxed_min8_peptide",
    "relaxed_min8_peptide_cov50",
    "relaxed_min8_peptide_cov70",
    "relaxed_min8_peptide_cov80",
]

THRESHOLD_LABELS = {
    "relaxed_min8_peptide": "Min-8 containment",
    "relaxed_min8_peptide_cov50": "Min-8 + ≥50%",
    "relaxed_min8_peptide_cov70": "Min-8 + ≥70%",
    "relaxed_min8_peptide_cov80": "Min-8 + ≥80%",
}


def main():
    if not IN_COUNTS.exists():
        raise FileNotFoundError(f"Missing counts table: {IN_COUNTS}")

    counts = pd.read_csv(IN_COUNTS, sep="\t")

    # Remove Any database from the displayed figure.
    counts = counts[counts["database"].isin(DATABASE_ORDER)].copy()

    counts["database"] = pd.Categorical(
        counts["database"],
        categories=DATABASE_ORDER,
        ordered=True,
    )

    counts["match_scope"] = pd.Categorical(
        counts["match_scope"],
        categories=THRESHOLD_ORDER,
        ordered=True,
    )

    counts = counts.sort_values(["match_scope", "database"])

    # Save the exact plotted values.
    plotted = counts[
        [
            "database",
            "database_label",
            "match_scope",
            "criterion_label",
            "supported_candidates",
            "total_candidates",
            "supported_fraction",
        ]
    ].copy()

    plotted.to_csv(
        OUTDIR / "Supplementary_Figure_S7_database_counts_no_any_plotted_values.tsv",
        sep="\t",
        index=False,
    )

    x = np.arange(len(THRESHOLD_ORDER))
    width = 0.24

    fig, ax = plt.subplots(figsize=(8.8, 5.2))

    for i, database in enumerate(DATABASE_ORDER):
        sub = counts[counts["database"].eq(database)].copy()

        values = [
            int(sub.loc[sub["match_scope"].eq(scope), "supported_candidates"].iloc[0])
            for scope in THRESHOLD_ORDER
        ]

        offset = (i - 1) * width

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
            )

    ax.set_xticks(x)
    ax.set_xticklabels(
        [THRESHOLD_LABELS[s] for s in THRESHOLD_ORDER],
        rotation=15,
        ha="right",
    )

    ax.set_ylabel("Supported EGG candidates")
    ax.set_xlabel("External-evidence matching criterion")
    ax.set_title("Supported EGG candidates retained under peptide-coverage thresholds")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(axis="y", alpha=0.25)

    fig.tight_layout()

    outstem = OUTDIR / "Supplementary_Figure_S7_database_counts_no_any"

    fig.savefig(outstem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outstem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(outstem.with_suffix(".svg"), bbox_inches="tight")

    plt.close(fig)

    print("Done.")
    print("Wrote:")
    print(outstem.with_suffix(".png"))
    print(outstem.with_suffix(".pdf"))
    print(outstem.with_suffix(".svg"))
    print()
    print("Plotted values:")
    print(plotted.to_string(index=False))


if __name__ == "__main__":
    main()
