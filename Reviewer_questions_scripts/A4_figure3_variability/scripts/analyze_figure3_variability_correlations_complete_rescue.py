#!/usr/bin/env python3

from pathlib import Path
import math
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


BASE = Path("TCGA_melanoma")

OUT_DIR = BASE / "rescue_final_analysis" / "reviewer_updates" / "reviewer_A4_fig3_variability"
PLOT_DIR = OUT_DIR / "plots"
TABLE_DIR = OUT_DIR / "tables"

PLOT_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)

FIG2_COUNTS_CANDIDATES = [
    BASE / "rescue_final_analysis/plots/figure2_old_logic_complete_rescue/tcga_complete_rescue_mutation_fusion_counts_merged.tsv",
    BASE / "rescue_final_analysis/plots/figure2_complete_rescue_export/tables/Figure2_complete_rescue_event_counts.tsv",
]

FIG3_TABLE_DIR = BASE / "rescue_final_analysis/plots/figure3_complete_rescue/tables"

FIG3_FILES = {
    "somatic": FIG3_TABLE_DIR / "figure3_neoantigens_somatic_complete_rescue.csv",
    "fusion": FIG3_TABLE_DIR / "figure3_neoantigens_fusion_complete_rescue.csv",
    "splicing": FIG3_TABLE_DIR / "figure3_neoantigens_splicing_complete_rescue.csv",
}


def find_existing(paths):
    for p in paths:
        if p.exists():
            return p
    raise FileNotFoundError("None of these files exist:\n" + "\n".join(str(p) for p in paths))


def norm_sample(x):
    x = str(x)
    x = x.replace("Sample_", "")
    return x


def read_fig2_counts():
    p = find_existing(FIG2_COUNTS_CANDIDATES)

    if p.suffix == ".tsv":
        df = pd.read_csv(p, sep="\t")
    else:
        df = pd.read_csv(p)

    df = df.copy()
    df["sample"] = df["sample"].map(norm_sample)

    required = ["sample", "snv_count", "indel_count", "fusion_count"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Figure 2 counts file {p} missing columns: {missing}. Columns: {df.columns.tolist()}")

    df["snv_count"] = pd.to_numeric(df["snv_count"], errors="coerce").fillna(0).astype(int)
    df["indel_count"] = pd.to_numeric(df["indel_count"], errors="coerce").fillna(0).astype(int)
    df["fusion_count"] = pd.to_numeric(df["fusion_count"], errors="coerce").fillna(0).astype(int)
    df["snv_indel_count"] = df["snv_count"] + df["indel_count"]

    return df[["sample", "snv_count", "indel_count", "snv_indel_count", "fusion_count"]], p


def read_category_counts(category):
    p = FIG3_FILES[category]
    if not p.exists():
        raise FileNotFoundError(f"Missing Figure 3 table for {category}: {p}")

    df = pd.read_csv(p)
    df = df.copy()
    df["sample"] = df["sample"].map(norm_sample)
    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(0).astype(int)

    wide = (
        df.pivot_table(
            index="sample",
            columns="MHC_Class",
            values="count",
            aggfunc="sum",
            fill_value=0,
        )
        .reset_index()
    )

    for col in ["MHC I", "MHC II"]:
        if col not in wide.columns:
            wide[col] = 0

    wide = wide.rename(columns={
        "MHC I": f"{category}_mhcI_neoantigens",
        "MHC II": f"{category}_mhcII_neoantigens",
    })

    wide[f"{category}_total_neoantigens"] = (
        wide[f"{category}_mhcI_neoantigens"] + wide[f"{category}_mhcII_neoantigens"]
    )

    return wide[
        [
            "sample",
            f"{category}_mhcI_neoantigens",
            f"{category}_mhcII_neoantigens",
            f"{category}_total_neoantigens",
        ]
    ], p


def set_pub_style():
    mpl.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 300,
        "font.size": 12,
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "grid.linestyle": "--",
        "grid.linewidth": 0.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def spearman_safe(df, x, y):
    sub = df[[x, y]].dropna().copy()
    n = len(sub)

    if n < 3:
        return n, np.nan, np.nan

    if sub[x].nunique() < 2 or sub[y].nunique() < 2:
        return n, np.nan, np.nan

    rho, p = spearmanr(sub[x], sub[y])
    return n, rho, p


def add_lm_line(ax, x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    ok = np.isfinite(x) & np.isfinite(y)

    if ok.sum() < 2:
        return

    x = x[ok]
    y = y[ok]

    if len(np.unique(x)) < 2:
        return

    coef = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 100)
    ys = coef[0] * xs + coef[1]
    ax.plot(xs, ys, color="black", linewidth=1.2, linestyle="--")


def scatter_plot(df, x, y, title, xlabel, ylabel, outbase, panel=None):
    n, rho, p = spearman_safe(df, x, y)

    fig = plt.figure(figsize=(5.2, 4.6))
    ax = plt.gca()

    ax.scatter(df[x], df[y], s=42, edgecolor="black", linewidth=0.5, alpha=0.85)
    add_lm_line(ax, df[x], df[y])

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if not math.isnan(rho):
        txt = f"Spearman ρ = {rho:.2f}\np = {p:.2g}\nn = {n}"
    else:
        txt = f"Spearman ρ = NA\np = NA\nn = {n}"

    ax.text(
        0.04,
        0.96,
        txt,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="0.7", alpha=0.9),
    )

    if panel:
        ax.text(
            -0.12,
            1.04,
            panel,
            transform=ax.transAxes,
            fontweight="bold",
            fontsize=14,
            va="bottom",
            ha="right",
        )

    fig.tight_layout()

    png = outbase.with_suffix(".png")
    svg = outbase.with_suffix(".svg")

    fig.savefig(png, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote: {png}")
    print(f"Wrote: {svg}")


def make_correlation_table(df):
    tests = [
        ("somatic_total_neoantigens", "snv_count", "Somatic total neoantigens", "SNV count"),
        ("somatic_total_neoantigens", "indel_count", "Somatic total neoantigens", "Indel count"),
        ("somatic_total_neoantigens", "snv_indel_count", "Somatic total neoantigens", "SNV + indel count"),
        ("somatic_mhcI_neoantigens", "snv_indel_count", "Somatic MHC I neoantigens", "SNV + indel count"),
        ("somatic_mhcII_neoantigens", "snv_indel_count", "Somatic MHC II neoantigens", "SNV + indel count"),
        ("fusion_total_neoantigens", "fusion_count", "Fusion total neoantigens", "Fusion event count"),
        ("fusion_mhcI_neoantigens", "fusion_count", "Fusion MHC I neoantigens", "Fusion event count"),
        ("fusion_mhcII_neoantigens", "fusion_count", "Fusion MHC II neoantigens", "Fusion event count"),
        ("fusion_total_neoantigens", "snv_count", "Fusion total neoantigens", "SNV count"),
        ("fusion_total_neoantigens", "snv_indel_count", "Fusion total neoantigens", "SNV + indel count"),
    ]

    rows = []

    for y, x, y_label, x_label in tests:
        if x not in df.columns or y not in df.columns:
            continue

        n, rho, p = spearman_safe(df, x, y)

        rows.append({
            "response": y,
            "response_label": y_label,
            "predictor": x,
            "predictor_label": x_label,
            "n_samples": n,
            "spearman_rho": rho,
            "p_value": p,
        })

    return pd.DataFrame(rows)


def main():
    set_pub_style()

    fig2, fig2_path = read_fig2_counts()

    somatic, somatic_path = read_category_counts("somatic")
    fusion, fusion_path = read_category_counts("fusion")
    splicing, splicing_path = read_category_counts("splicing")

    df = fig2.merge(somatic, on="sample", how="left")
    df = df.merge(fusion, on="sample", how="left")
    df = df.merge(splicing, on="sample", how="left")

    neo_cols = [c for c in df.columns if c.endswith("_neoantigens")]
    for c in neo_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)

    df = df.sort_values("sample").reset_index(drop=True)

    merged_path = TABLE_DIR / "A4_complete_rescue_mutation_fusion_neoantigen_counts_by_sample.tsv"
    df.to_csv(merged_path, sep="\t", index=False)

    corr = make_correlation_table(df)
    corr_path = TABLE_DIR / "A4_complete_rescue_spearman_correlations.tsv"
    corr.to_csv(corr_path, sep="\t", index=False)

    input_paths = pd.DataFrame([
        {"input": "figure2_mutation_fusion_counts", "path": str(fig2_path)},
        {"input": "figure3_somatic_counts", "path": str(somatic_path)},
        {"input": "figure3_fusion_counts", "path": str(fusion_path)},
        {"input": "figure3_splicing_counts", "path": str(splicing_path)},
    ])
    input_paths.to_csv(TABLE_DIR / "A4_input_files_used.tsv", sep="\t", index=False)

    scatter_plot(
        df,
        x="snv_indel_count",
        y="somatic_total_neoantigens",
        title="Somatic neoantigen yield vs SNV + indel mutations",
        xlabel="SNV + indel mutations",
        ylabel="Somatic neoantigen count",
        outbase=PLOT_DIR / "A4A_somatic_neoantigens_vs_snv_indel_count",
        panel="A",
    )

    scatter_plot(
        df,
        x="snv_count",
        y="somatic_total_neoantigens",
        title="Somatic neoantigen yield vs SNV mutations",
        xlabel="SNV mutations",
        ylabel="Somatic neoantigen count",
        outbase=PLOT_DIR / "A4B_somatic_neoantigens_vs_snv_count",
        panel="B",
    )

    scatter_plot(
        df,
        x="fusion_count",
        y="fusion_total_neoantigens",
        title="Fusion neoantigen yield vs fusion burden",
        xlabel="STAR-Fusion event count",
        ylabel="Fusion neoantigen count",
        outbase=PLOT_DIR / "A4C_fusion_neoantigens_vs_fusion_count",
        panel="C",
    )

    scatter_plot(
        df,
        x="snv_count",
        y="fusion_total_neoantigens",
        title="Fusion neoantigen yield vs SNV burden",
        xlabel="SNV mutations",
        ylabel="Fusion neoantigen count",
        outbase=PLOT_DIR / "A4D_fusion_neoantigens_vs_snv_count",
        panel="D",
    )

    print()
    print("Wrote tables:")
    print(merged_path)
    print(corr_path)
    print(TABLE_DIR / "A4_input_files_used.tsv")
    print()
    print("Correlation table:")
    print(corr.to_string(index=False))


if __name__ == "__main__":
    main()
