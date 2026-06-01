#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr


BASE = Path("TCGA_melanoma")

FINAL_DIR = BASE / "epitopes_prioritisation_complete_rescue" / "final_epitopes"

OUT_DIR = BASE / "rescue_final_analysis" / "reviewer_updates" / "reviewer_A8_peptide_length"
PLOT_DIR = OUT_DIR / "plots"
TABLE_DIR = OUT_DIR / "tables"

PLOT_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)

PEPTIDE_COL_CANDIDATES = [
    "MT.Epitope.Seq",
    "MT Epitope Seq",
    "Epitope.Seq",
    "Epitope Seq",
    "Peptide",
    "peptide",
    "Mutation",
]

MHC_COL_CANDIDATES = [
    "HLA.Class",
    "HLA Class",
    "MHC_Class",
    "MHC Class",
    "mhc_class",
]

RANK_COL_CANDIDATES = [
    "Borda_Rank",
    "Borda Rank",
    "Borda.Rank",
    "Final_Rank",
    "Final Rank",
    "Rank",
]


def find_col(df, candidates, required=True):
    for c in candidates:
        if c in df.columns:
            return c
    if required:
        raise ValueError(f"Could not find any of {candidates}. Available columns: {df.columns.tolist()}")
    return None


def sample_from_filename(path):
    name = path.name
    return name.replace("_epitopes_final.csv", "")


def norm_mhc(x):
    s = str(x).strip()

    if s in {"MHC_Class_I", "MHC I", "MHC Class I", "Class I", "I"}:
        return "MHC I"

    if s in {"MHC_Class_II", "MHC II", "MHC Class II", "Class II", "II"}:
        return "MHC II"

    s_upper = s.upper()

    if "CLASS_II" in s_upper or "CLASS II" in s_upper or s_upper.endswith("II"):
        return "MHC II"

    if "CLASS_I" in s_upper or "CLASS I" in s_upper or s_upper.endswith("I"):
        return "MHC I"

    return s


def peptide_len(x):
    s = str(x).strip()
    if s == "" or s.lower() in {"nan", "none", "na"}:
        return np.nan
    return len(s)


def read_final_tables():
    files = sorted(FINAL_DIR.glob("*_epitopes_final.csv"))

    if not files:
        raise FileNotFoundError(f"No final epitope tables found in {FINAL_DIR}")

    rows = []

    for f in files:
        df = pd.read_csv(f, dtype=str, low_memory=False).fillna("")
        df["sample"] = sample_from_filename(f)
        df["source_file"] = str(f)
        rows.append(df)

    all_df = pd.concat(rows, ignore_index=True)

    peptide_col = find_col(all_df, PEPTIDE_COL_CANDIDATES)
    mhc_col = find_col(all_df, MHC_COL_CANDIDATES)
    rank_col = find_col(all_df, RANK_COL_CANDIDATES)

    all_df["peptide_sequence_for_length"] = all_df[peptide_col].astype(str).str.strip()
    all_df["peptide_length"] = all_df["peptide_sequence_for_length"].map(peptide_len)
    all_df["MHC_Class_clean"] = all_df[mhc_col].map(norm_mhc)
    all_df["rank_numeric"] = pd.to_numeric(all_df[rank_col], errors="coerce")

    all_df = all_df[
        all_df["MHC_Class_clean"].isin(["MHC I", "MHC II"])
        & all_df["peptide_length"].notna()
        & all_df["rank_numeric"].notna()
    ].copy()

    all_df["peptide_length"] = all_df["peptide_length"].astype(int)

    all_df = all_df.sort_values(["sample", "rank_numeric"]).reset_index(drop=True)
    all_df["rank_order_within_sample"] = all_df.groupby("sample").cumcount() + 1
    all_df["is_top10_per_sample"] = all_df["rank_order_within_sample"] <= 10

    all_df = all_df.sort_values(["sample", "MHC_Class_clean", "rank_numeric"]).reset_index(drop=True)
    all_df["rank_order_within_sample_class"] = all_df.groupby(["sample", "MHC_Class_clean"]).cumcount() + 1
    all_df["n_within_sample_class"] = all_df.groupby(["sample", "MHC_Class_clean"])["rank_order_within_sample_class"].transform("max")

    all_df["rank_percentile_within_sample_class"] = np.where(
        all_df["n_within_sample_class"] > 1,
        (all_df["rank_order_within_sample_class"] - 1) / (all_df["n_within_sample_class"] - 1),
        0.0,
    )

    all_df["n_within_sample"] = all_df.groupby("sample")["rank_order_within_sample"].transform("max")
    all_df["rank_percentile_within_sample"] = np.where(
        all_df["n_within_sample"] > 1,
        (all_df["rank_order_within_sample"] - 1) / (all_df["n_within_sample"] - 1),
        0.0,
    )

    meta = {
        "peptide_col": peptide_col,
        "mhc_col": mhc_col,
        "rank_col": rank_col,
    }

    return all_df, meta


def safe_iqr(x):
    x = pd.Series(x).dropna()
    if len(x) == 0:
        return np.nan, np.nan
    return np.percentile(x, 25), np.percentile(x, 75)


def safe_mannwhitney(a, b):
    a = pd.Series(a).dropna()
    b = pd.Series(b).dropna()

    if len(a) < 1 or len(b) < 1:
        return np.nan

    return mannwhitneyu(a, b, alternative="two-sided").pvalue


def safe_spearman(x, y):
    d = pd.DataFrame({"x": x, "y": y}).dropna()

    if len(d) < 3:
        return len(d), np.nan, np.nan

    if d["x"].nunique() < 2 or d["y"].nunique() < 2:
        return len(d), np.nan, np.nan

    rho, p = spearmanr(d["x"], d["y"])
    return len(d), rho, p


def make_summary(df):
    rows = []

    for mhc_class, g in df.groupby("MHC_Class_clean", dropna=False):
        top = g[g["is_top10_per_sample"]].copy()
        non_top = g[~g["is_top10_per_sample"]].copy()

        all_q1, all_q3 = safe_iqr(g["peptide_length"])
        top_q1, top_q3 = safe_iqr(top["peptide_length"])
        non_q1, non_q3 = safe_iqr(non_top["peptide_length"])

        mw_p = safe_mannwhitney(top["peptide_length"], non_top["peptide_length"])

        n_class_pct, rho_class_pct, p_class_pct = safe_spearman(
            g["peptide_length"],
            g["rank_percentile_within_sample_class"],
        )

        n_overall_pct, rho_overall_pct, p_overall_pct = safe_spearman(
            g["peptide_length"],
            g["rank_percentile_within_sample"],
        )

        rows.append({
            "MHC_Class": mhc_class,
            "n_all_candidates": int(len(g)),
            "median_length_all": float(np.median(g["peptide_length"])) if len(g) else np.nan,
            "q1_length_all": all_q1,
            "q3_length_all": all_q3,
            "n_top10_candidates": int(len(top)),
            "median_length_top10": float(np.median(top["peptide_length"])) if len(top) else np.nan,
            "q1_length_top10": top_q1,
            "q3_length_top10": top_q3,
            "n_non_top10_candidates": int(len(non_top)),
            "median_length_non_top10": float(np.median(non_top["peptide_length"])) if len(non_top) else np.nan,
            "q1_length_non_top10": non_q1,
            "q3_length_non_top10": non_q3,
            "mannwhitney_top10_vs_non_top10_p": mw_p,
            "spearman_n_length_vs_class_rank_percentile": n_class_pct,
            "spearman_rho_length_vs_class_rank_percentile": rho_class_pct,
            "spearman_p_length_vs_class_rank_percentile": p_class_pct,
            "spearman_n_length_vs_overall_rank_percentile": n_overall_pct,
            "spearman_rho_length_vs_overall_rank_percentile": rho_overall_pct,
            "spearman_p_length_vs_overall_rank_percentile": p_overall_pct,
            "max_length_all": int(g["peptide_length"].max()) if len(g) else np.nan,
            "max_length_top10": int(top["peptide_length"].max()) if len(top) else np.nan,
        })

    return pd.DataFrame(rows)


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


def draw_box(ax, data, pos, color):
    bp = ax.boxplot(
        data,
        positions=[pos],
        widths=0.55,
        patch_artist=True,
        showfliers=True,
    )

    for patch in bp["boxes"]:
        patch.set_facecolor(color)
        patch.set_alpha(0.85)
        patch.set_edgecolor("black")
        patch.set_linewidth(0.9)

    for k in ["medians", "whiskers", "caps"]:
        for artist in bp[k]:
            artist.set_color("black")
            artist.set_linewidth(1.0)

    for flier in bp["fliers"]:
        flier.set_marker("o")
        flier.set_markersize(3)
        flier.set_markerfacecolor(color)
        flier.set_markeredgecolor("black")
        flier.set_alpha(0.5)


def make_plot(df, summary):
    set_pub_style()

    fig = plt.figure(figsize=(7.4, 4.8))
    ax = plt.gca()

    plot_specs = [
        ("MHC I", "All candidates", 0, "#1f77b4"),
        ("MHC I", "Top 10 per sample", 1, "#ff7f0e"),
        ("MHC II", "All candidates", 3, "#1f77b4"),
        ("MHC II", "Top 10 per sample", 4, "#ff7f0e"),
    ]

    for mhc, group, pos, color in plot_specs:
        sub = df[df["MHC_Class_clean"] == mhc].copy()

        if group == "Top 10 per sample":
            sub = sub[sub["is_top10_per_sample"]]

        draw_box(ax, sub["peptide_length"].values, pos, color)

    ax.set_xticks([0.5, 3.5])
    ax.set_xticklabels(["MHC I", "MHC II"])
    ax.set_ylabel("Peptide length (amino acids)")
    ax.set_title("Peptide length distribution before and after prioritisation")

    ax.legend(
        handles=[
            plt.Rectangle((0, 0), 1, 1, facecolor="#1f77b4", edgecolor="black", alpha=0.85, label="All candidates"),
            plt.Rectangle((0, 0), 1, 1, facecolor="#ff7f0e", edgecolor="black", alpha=0.85, label="Top 10 per sample"),
        ],
        frameon=False,
        loc="upper left",
    )

    ymax = max(35, int(df["peptide_length"].max()) + 3)
    ax.set_ylim(0, ymax)

    # Add p-value annotations for top10 vs non-top10.
    for mhc, xpos in [("MHC I", 0.5), ("MHC II", 3.5)]:
        row = summary[summary["MHC_Class"] == mhc]
        if len(row):
            pval = row["mannwhitney_top10_vs_non_top10_p"].iloc[0]
            if pd.notna(pval):
                label = f"MW p = {pval:.3g}"
            else:
                label = "MW p = NA"

            ax.text(
                xpos,
                ymax * 0.96,
                label,
                ha="center",
                va="top",
                fontsize=10,
            )

    fig.tight_layout()

    png = PLOT_DIR / "A8_peptide_length_all_vs_top10_by_MHC_complete_rescue.png"
    svg = PLOT_DIR / "A8_peptide_length_all_vs_top10_by_MHC_complete_rescue.svg"

    fig.savefig(png, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)

    print("Wrote:")
    print(png)
    print(svg)


def main():
    df, meta = read_final_tables()

    dataset_path = TABLE_DIR / "A8_peptide_length_rank_dataset_complete_rescue.tsv"
    summary_path = TABLE_DIR / "A8_peptide_length_bias_summary_complete_rescue.tsv"
    meta_path = TABLE_DIR / "A8_input_columns_used.tsv"

    df.to_csv(dataset_path, sep="\t", index=False)

    summary = make_summary(df)
    summary.to_csv(summary_path, sep="\t", index=False)

    pd.DataFrame([meta]).to_csv(meta_path, sep="\t", index=False)

    make_plot(df, summary)

    print()
    print("Wrote tables:")
    print(dataset_path)
    print(summary_path)
    print(meta_path)
    print()
    print("Summary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
