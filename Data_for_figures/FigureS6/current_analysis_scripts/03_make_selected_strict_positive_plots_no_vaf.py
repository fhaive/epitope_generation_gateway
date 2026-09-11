#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(".").resolve()
IN_SUMMARY = ROOT / "ranking_vs_ic50_exact_relaxed_outputs" / "summary_ranking_delta_vs_ic50.tsv"

OUTDIR = ROOT / "selected_strict_positive_plots_no_vaf"
OUTDIR.mkdir(exist_ok=True)

ENDPOINT = "strict_positive_experimental"

PLOTS_TO_MAKE = [
    ("exact_peptide", "percent", "strict_positive_experimental__exact_peptide__percent"),
    ("relaxed_min8_peptide", "percent", "strict_positive_experimental__relaxed_min8_peptide__percent"),
    ("exact_peptide", "topn", "strict_positive_experimental__exact_peptide__topn"),
    ("relaxed_min8_peptide", "topn", "strict_positive_experimental__relaxed_min8_peptide__topn"),
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

MARKERS = ["o", "s", "^", "D", "v", "P", "X"]
LINESTYLES = ["-", "--", "-.", ":", "-", "--", "-."]


def savefig(path):
    plt.tight_layout()
    plt.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()


def main():
    df = pd.read_csv(IN_SUMMARY, sep="\t")

    # Remove Tumor DNA VAF completely.
    df = df[df["ranking_method"] != "tumor_dna_vaf"].copy()

    # Keep only strict positive experimental endpoint.
    df = df[df["endpoint"] == ENDPOINT].copy()

    # Save the filtered plotting source table for transparency.
    df.to_csv(OUTDIR / "selected_strict_positive_plot_source.tsv", sep="\t", index=False)

    for match_scope, top_type, outstem in PLOTS_TO_MAKE:
        sub = df[
            (df["match_scope"] == match_scope)
            & (df["top_type"] == top_type)
        ].copy()

        if sub.empty:
            print(f"[WARN] No data for {outstem}")
            continue

        methods = [m for m in METHOD_ORDER if m in set(sub["ranking_method"])]
        if not methods:
            print(f"[WARN] No ranking methods found for {outstem}")
            continue

        cutoffs = sorted(sub["cutoff"].unique())
        x_values = np.array(cutoffs, dtype=float)

        # Small horizontal offsets make overlapping points/lines visible.
        if len(methods) > 1:
            step = np.median(np.diff(sorted(cutoffs))) if len(cutoffs) > 1 else 1.0
            offsets = np.linspace(-0.18 * step, 0.18 * step, len(methods))
        else:
            offsets = [0.0]

        plt.figure(figsize=(10.5, 6.2))

        for i, method in enumerate(methods):
            mdf = sub[sub["ranking_method"] == method].sort_values("cutoff").copy()

            # Convert from proportion to percentage points.
            y = 100.0 * pd.to_numeric(mdf["median_delta_hit_rate_vs_ic50"], errors="coerce")

            x = pd.to_numeric(mdf["cutoff"], errors="coerce").to_numpy(dtype=float) + offsets[i]

            label = METHOD_LABELS.get(method, mdf["ranking_label"].iloc[0])

            plt.plot(
                x,
                y,
                marker=MARKERS[i % len(MARKERS)],
                linestyle=LINESTYLES[i % len(LINESTYLES)],
                linewidth=2.0,
                markersize=6.5,
                alpha=0.90,
                label=label,
            )

            # Add final-point label to help where lines overlap.
            if len(x) > 0 and not pd.isna(y.iloc[-1]):
                plt.text(
                    x[-1] + 0.15,
                    y.iloc[-1],
                    label,
                    fontsize=8,
                    va="center",
                )

        plt.axhline(0, linestyle="--", linewidth=1.2)

        if top_type == "percent":
            xlabel = "Top-percent cutoff within each patient"
            title_suffix = "top-percent cutoffs"
        else:
            xlabel = "Top-N candidates within each patient"
            title_suffix = "top-N cutoffs"

        if match_scope == "exact_peptide":
            scope_label = "exact peptide match"
        elif match_scope == "relaxed_min8_peptide":
            scope_label = "relaxed min-8 peptide match"
        else:
            scope_label = match_scope

        plt.xlabel(xlabel)
        plt.ylabel("Median delta hit rate vs IC50 ranking\npercentage points")
        plt.title(
            "Strict positive experimental evidence\n"
            f"{scope_label}; {title_suffix}"
        )
        plt.legend(loc="best", fontsize=8, frameon=True)
        plt.grid(axis="y", alpha=0.25)

        savefig(OUTDIR / outstem)

        print(f"Wrote: {OUTDIR / (outstem + '.png')}")
        print(f"Wrote: {OUTDIR / (outstem + '.pdf')}")

    # Package only these selected plots.
    import tarfile
    tar_path = ROOT / "selected_strict_positive_plots_no_vaf.tar.gz"
    with tarfile.open(tar_path, "w:gz") as tar:
        tar.add(OUTDIR, arcname=OUTDIR.name)

    print()
    print("Selected plot folder:", OUTDIR)
    print("Compressed package:", tar_path)


if __name__ == "__main__":
    main()
