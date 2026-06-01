#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(".").resolve()

SUMMARY = ROOT / "ranking_vs_ic50_exact_relaxed_outputs" / "summary_ranking_delta_vs_ic50.tsv"
INVENTORY = ROOT / "candidate_evidence_inventory.tsv"
DB_SUMMARY = ROOT / "strict_positive_match_database_summary" / "strict_positive_match_counts_by_database.tsv"

OUTDIR = ROOT / "reply_selected_strict_positive_plots_clean"
OUTDIR.mkdir(exist_ok=True)

ENDPOINT = "strict_positive_experimental"

PLOTS = [
    ("relaxed_min8_peptide", "topn", "strict_positive_experimental__relaxed_min8_peptide__topn"),
    ("exact_peptide", "topn", "strict_positive_experimental__exact_peptide__topn"),
    ("relaxed_min8_peptide", "percent", "strict_positive_experimental__relaxed_min8_peptide__percent"),
    ("exact_peptide", "percent", "strict_positive_experimental__exact_peptide__percent"),
]

METHOD_ORDER = [
    "borda",
    "depmap_survivability",
    "net_betweenness",
    "net_degree",
    "net_impact",
    "net_strength",
    "net_wci",
]

METHOD_LABELS = {
    "borda": "Borda",
    "depmap_survivability": "DepMap survivability",
    "net_betweenness": "Network betweenness",
    "net_degree": "Network degree",
    "net_impact": "Network impact",
    "net_strength": "Network strength",
    "net_wci": "Network WCI",
}

MARKERS = {
    "borda": "o",
    "depmap_survivability": "s",
    "net_betweenness": "^",
    "net_degree": "D",
    "net_impact": "v",
    "net_strength": "P",
    "net_wci": "X",
}

LINESTYLES = {
    "borda": "-",
    "depmap_survivability": "--",
    "net_betweenness": "-.",
    "net_degree": ":",
    "net_impact": "-",
    "net_strength": "--",
    "net_wci": "-.",
}


def savefig(path):
    plt.tight_layout()
    plt.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()


def clean_inventory_value(inv, metric):
    hit = inv.loc[inv["metric"] == metric, "value"]
    if hit.empty:
        return None
    return int(float(hit.iloc[0]))


def main():
    df = pd.read_csv(SUMMARY, sep="\t")

    # Remove Tumor DNA VAF from plots and tables.
    df = df[df["ranking_method"] != "tumor_dna_vaf"].copy()

    # Keep only strict positive experimental endpoint.
    df = df[df["endpoint"] == ENDPOINT].copy()

    # Convert delta from fraction to percentage points.
    df["median_delta_hit_rate_vs_ic50_pct_points"] = (
        100.0 * pd.to_numeric(df["median_delta_hit_rate_vs_ic50"], errors="coerce")
    )

    # Save all exact plotted values.
    plot_source = df[
        [
            "endpoint",
            "match_scope",
            "top_type",
            "cutoff",
            "ranking_method",
            "ranking_label",
            "ranking_column",
            "ranking_direction",
            "n_patients",
            "median_delta_hit_rate_vs_ic50",
            "median_delta_hit_rate_vs_ic50_pct_points",
            "patients_better_than_ic50",
            "patients_same_as_ic50",
            "patients_worse_than_ic50",
        ]
    ].copy()

    plot_source.to_csv(
        OUTDIR / "numbers_used_for_selected_strict_positive_plots.tsv",
        sep="\t",
        index=False,
    )

    # Optional concise summary of database evidence counts.
    if DB_SUMMARY.exists():
        db = pd.read_csv(DB_SUMMARY, sep="\t")
        db.to_csv(OUTDIR / "strict_positive_database_match_counts_used_for_reply.tsv", sep="\t", index=False)

    # Candidate evidence inventory copied for methods transparency.
    if INVENTORY.exists():
        inv = pd.read_csv(INVENTORY, sep="\t")
        inv.to_csv(OUTDIR / "candidate_evidence_inventory_used_for_reply.tsv", sep="\t", index=False)

    for match_scope, top_type, outstem in PLOTS:
        sub = df[
            (df["match_scope"] == match_scope)
            & (df["top_type"] == top_type)
        ].copy()

        if sub.empty:
            print(f"[WARN] no data for {outstem}")
            continue

        methods = [m for m in METHOD_ORDER if m in set(sub["ranking_method"])]
        cutoffs = sorted(sub["cutoff"].unique())

        plt.figure(figsize=(9.5, 5.8))

        for method in methods:
            mdf = sub[sub["ranking_method"] == method].sort_values("cutoff")
            x = pd.to_numeric(mdf["cutoff"], errors="coerce")
            y = pd.to_numeric(mdf["median_delta_hit_rate_vs_ic50_pct_points"], errors="coerce")

            plt.plot(
                x,
                y,
                label=METHOD_LABELS.get(method, method),
                marker=MARKERS.get(method, "o"),
                linestyle=LINESTYLES.get(method, "-"),
                linewidth=2.0,
                markersize=6.0,
                alpha=0.9,
            )

        plt.axhline(0, linestyle="--", linewidth=1.2)

        if top_type == "topn":
            xlabel = "Top-N candidates within each patient"
            title_cutoff = "top-N cutoffs"
        else:
            xlabel = "Top-percent cutoff within each patient"
            title_cutoff = "top-percent cutoffs"

        if match_scope == "relaxed_min8_peptide":
            scope_title = "relaxed min-8 peptide match"
        else:
            scope_title = "exact peptide match"

        plt.xlabel(xlabel)
        plt.ylabel("Median delta hit rate vs IC50 ranking\npercentage points")
        plt.title(f"Strict positive experimental evidence\n{scope_title}; {title_cutoff}")
        plt.grid(axis="y", alpha=0.25)
        plt.legend(loc="best", fontsize=8, frameon=True)

        savefig(OUTDIR / outstem)

        # Write a small per-plot table with exactly what went into the image.
        sub[
            [
                "match_scope",
                "top_type",
                "cutoff",
                "ranking_method",
                "ranking_label",
                "n_patients",
                "median_delta_hit_rate_vs_ic50_pct_points",
                "patients_better_than_ic50",
                "patients_same_as_ic50",
                "patients_worse_than_ic50",
            ]
        ].sort_values(["cutoff", "ranking_method"]).to_csv(
            OUTDIR / f"{outstem}.numbers.tsv",
            sep="\t",
            index=False,
        )

        print(f"Wrote {OUTDIR / (outstem + '.png')}")
        print(f"Wrote {OUTDIR / (outstem + '.numbers.tsv')}")

    # Package.
    import tarfile
    tar_path = ROOT / "reply_selected_strict_positive_plots_clean.tar.gz"
    with tarfile.open(tar_path, "w:gz") as tar:
        tar.add(OUTDIR, arcname=OUTDIR.name)

    print()
    print("Output folder:")
    print(OUTDIR)
    print()
    print("Compressed package:")
    print(tar_path)


if __name__ == "__main__":
    main()
