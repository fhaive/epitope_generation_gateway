#!/usr/bin/env python3

from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr


BASE = Path("TCGA_melanoma")

FINAL_DIR = BASE / "epitopes_prioritisation_complete_rescue" / "final_epitopes"
PVAC_DIR = BASE / "2A_somatic_mutation_epitopes_complete_rescue" / "pvacSeq"

OUT_DIR = BASE / "rescue_final_analysis" / "reviewer_updates" / "reviewer_A9_vaf_clonality"
PLOT_DIR = OUT_DIR / "plots"
TABLE_DIR = OUT_DIR / "tables"

PLOT_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)

# Simplified clonality proxy:
# assuming tumor purity ~1, VAF >= 0.25 is broadly compatible with clonal/near-clonal events.
CLONAL_VAF_THRESHOLD = 0.25


FINAL_SOURCE_COLS = ["Mutation.Source", "Mutation Source", "Source"]
FINAL_MHC_COLS = ["HLA.Class", "HLA Class", "MHC_Class", "MHC Class", "mhc_class"]
FINAL_RANK_COLS = ["Borda_Rank", "Borda Rank", "Borda.Rank", "Final_Rank", "Final Rank", "Rank"]
FINAL_GENE_COLS = ["Gene.Name", "Gene Name", "Gene", "Gene.Symbol"]
FINAL_PEPTIDE_COLS = ["MT.Epitope.Seq", "MT Epitope Seq", "Epitope Seq", "Epitope.Seq", "Peptide"]
FINAL_HLA_COLS = ["HLA.Allele", "HLA Allele", "HLA"]

PVAC_GENE_COLS = ["Gene Name", "Gene", "Gene.Name", "Gene.Symbol"]
PVAC_PEPTIDE_COLS = ["MT Epitope Seq", "MT.Epitope.Seq", "Epitope Seq", "Epitope.Seq", "Peptide"]
PVAC_HLA_COLS = ["HLA Allele", "HLA.Allele", "HLA"]
PVAC_VAF_COLS = ["Tumor DNA VAF", "Tumor.DNA.VAF", "Tumor_DNA_VAF"]


def find_col(df, candidates, required=True):
    for c in candidates:
        if c in df.columns:
            return c
    if required:
        raise ValueError(
            f"Could not find any of {candidates}. Available columns:\n"
            + "\n".join(df.columns.astype(str))
        )
    return None


def norm_sample(x):
    x = str(x).strip()
    x = x.replace("_epitopes_final.csv", "")
    if x.startswith("Sample_"):
        return x
    m = re.search(r"(TCGA-[A-Z0-9]+-[A-Z0-9]+)", x)
    if m:
        return "Sample_" + m.group(1)
    return x


def sample_from_final_path(path):
    return norm_sample(path.name.replace("_epitopes_final.csv", ""))


def sample_from_pvac_path(path):
    # TCGA_melanoma/.../pvacSeq/Sample_TCGA-.../MHC_Class_I/file.filtered.tsv
    return norm_sample(path.parents[1].name)


def mhc_from_pvac_path(path):
    d = path.parent.name
    if d == "MHC_Class_I":
        return "MHC I"
    if d == "MHC_Class_II":
        return "MHC II"
    return d


def norm_mhc(x):
    s = str(x).strip()

    if s in {"MHC_Class_I", "MHC I", "MHC Class I", "Class I", "I"}:
        return "MHC I"
    if s in {"MHC_Class_II", "MHC II", "MHC Class II", "Class II", "II"}:
        return "MHC II"

    su = s.upper()
    if "CLASS_II" in su or "CLASS II" in su or su.endswith("II"):
        return "MHC II"
    if "CLASS_I" in su or "CLASS I" in su or su.endswith("I"):
        return "MHC I"

    return s


def norm_hla(x):
    return str(x).strip()


def clean_string(x):
    return str(x).strip()


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


def read_pvacseq_vaf_lookup():
    files = sorted(PVAC_DIR.glob("Sample_*/MHC_Class_*/*.filtered.tsv"))
    if not files:
        raise FileNotFoundError(f"No pVACseq filtered TSVs found under {PVAC_DIR}")

    rows = []
    column_records = []

    for f in files:
        df = pd.read_csv(f, sep="\t", dtype=str, low_memory=False).fillna("")

        gene_col = find_col(df, PVAC_GENE_COLS)
        peptide_col = find_col(df, PVAC_PEPTIDE_COLS)
        hla_col = find_col(df, PVAC_HLA_COLS)
        vaf_col = find_col(df, PVAC_VAF_COLS)

        column_records.append({
            "file": str(f),
            "gene_col": gene_col,
            "peptide_col": peptide_col,
            "hla_col": hla_col,
            "vaf_col": vaf_col,
        })

        tmp = pd.DataFrame({
            "sample_key": sample_from_pvac_path(f),
            "mhc_key": mhc_from_pvac_path(f),
            "gene_key": df[gene_col].map(clean_string),
            "peptide_key": df[peptide_col].map(clean_string),
            "hla_key": df[hla_col].map(norm_hla),
            "tumor_dna_vaf": pd.to_numeric(df[vaf_col], errors="coerce"),
            "pvac_file": str(f),
        })

        rows.append(tmp)

    pvac = pd.concat(rows, ignore_index=True)
    pvac = pvac[pvac["tumor_dna_vaf"].notna()].copy()

    if pvac["tumor_dna_vaf"].dropna().max() > 1.5:
        pvac["tumor_dna_vaf"] = pvac["tumor_dna_vaf"] / 100.0

    key_cols = ["sample_key", "mhc_key", "gene_key", "peptide_key", "hla_key"]

    lookup = (
        pvac
        .groupby(key_cols, as_index=False)
        .agg(
            tumor_dna_vaf=("tumor_dna_vaf", "median"),
            tumor_dna_vaf_min=("tumor_dna_vaf", "min"),
            tumor_dna_vaf_max=("tumor_dna_vaf", "max"),
            n_pvac_rows_for_key=("tumor_dna_vaf", "size"),
            pvac_files=("pvac_file", lambda x: ";".join(sorted(set(map(str, x)))))
        )
    )

    pd.DataFrame(column_records).to_csv(TABLE_DIR / "A9_pvacseq_columns_used.tsv", sep="\t", index=False)
    pvac.to_csv(TABLE_DIR / "A9_pvacseq_vaf_rows_before_grouping.tsv", sep="\t", index=False)
    lookup.to_csv(TABLE_DIR / "A9_pvacseq_vaf_lookup_grouped.tsv", sep="\t", index=False)

    return lookup


def read_final_somatic_tables():
    files = sorted(FINAL_DIR.glob("*_epitopes_final.csv"))
    if not files:
        raise FileNotFoundError(f"No final epitope tables found in {FINAL_DIR}")

    rows = []
    column_records = []

    for f in files:
        df = pd.read_csv(f, dtype=str, low_memory=False).fillna("")

        source_col = find_col(df, FINAL_SOURCE_COLS)
        mhc_col = find_col(df, FINAL_MHC_COLS)
        rank_col = find_col(df, FINAL_RANK_COLS)
        gene_col = find_col(df, FINAL_GENE_COLS)
        peptide_col = find_col(df, FINAL_PEPTIDE_COLS)
        hla_col = find_col(df, FINAL_HLA_COLS)

        sample = sample_from_final_path(f)

        column_records.append({
            "file": str(f),
            "source_col": source_col,
            "mhc_col": mhc_col,
            "rank_col": rank_col,
            "gene_col": gene_col,
            "peptide_col": peptide_col,
            "hla_col": hla_col,
        })

        df["sample"] = sample
        df["source_file"] = str(f)
        df["source_clean"] = df[source_col].astype(str).str.strip().str.lower()
        df["MHC_Class_clean"] = df[mhc_col].map(norm_mhc)
        df["rank_numeric"] = pd.to_numeric(df[rank_col], errors="coerce")

        df["sample_key"] = sample
        df["mhc_key"] = df["MHC_Class_clean"]
        df["gene_key"] = df[gene_col].map(clean_string)
        df["peptide_key"] = df[peptide_col].map(clean_string)
        df["hla_key"] = df[hla_col].map(norm_hla)

        rows.append(df)

    final = pd.concat(rows, ignore_index=True)
    final = final[final["rank_numeric"].notna()].copy()

    final = final.sort_values(["sample", "rank_numeric"]).reset_index(drop=True)
    final["rank_order_within_sample"] = final.groupby("sample").cumcount() + 1
    final["is_top10_overall_per_sample"] = final["rank_order_within_sample"] <= 10

    final["n_within_sample"] = final.groupby("sample")["rank_order_within_sample"].transform("max")
    final["rank_percentile_within_sample"] = np.where(
        final["n_within_sample"] > 1,
        (final["rank_order_within_sample"] - 1) / (final["n_within_sample"] - 1),
        0.0,
    )

    somatic = final[
        (final["source_clean"] == "somatic")
        & final["MHC_Class_clean"].isin(["MHC I", "MHC II"])
    ].copy()

    pd.DataFrame(column_records).to_csv(TABLE_DIR / "A9_final_table_columns_used.tsv", sep="\t", index=False)

    return final, somatic


def attach_vaf_to_final_somatic(somatic, lookup):
    key_cols = ["sample_key", "mhc_key", "gene_key", "peptide_key", "hla_key"]

    merged = somatic.merge(
        lookup,
        on=key_cols,
        how="left",
    )

    # Fallback relaxed match if exact HLA key fails: sample + MHC + gene + peptide.
    unmatched = merged["tumor_dna_vaf"].isna()

    if unmatched.any():
        relaxed_lookup = (
            lookup
            .groupby(["sample_key", "mhc_key", "gene_key", "peptide_key"], as_index=False)
            .agg(
                tumor_dna_vaf_relaxed=("tumor_dna_vaf", "median"),
                n_relaxed_matches=("tumor_dna_vaf", "size")
            )
        )

        merged = merged.merge(
            relaxed_lookup,
            on=["sample_key", "mhc_key", "gene_key", "peptide_key"],
            how="left",
        )

        use_relaxed = merged["tumor_dna_vaf"].isna() & merged["tumor_dna_vaf_relaxed"].notna()
        merged.loc[use_relaxed, "tumor_dna_vaf"] = merged.loc[use_relaxed, "tumor_dna_vaf_relaxed"]
        merged.loc[use_relaxed, "vaf_match_mode"] = "relaxed_without_hla"

    if "vaf_match_mode" not in merged.columns:
        merged["vaf_match_mode"] = ""

    merged["vaf_match_mode"] = merged["vaf_match_mode"].fillna("")

    merged.loc[merged["vaf_match_mode"] == "", "vaf_match_mode"] = np.where(
        merged["tumor_dna_vaf"].notna(),
        "exact_sample_mhc_gene_peptide_hla",
        "unmatched_to_pvacseq_vaf",
    )

    matched = merged[merged["tumor_dna_vaf"].notna()].copy()
    matched["clonality_proxy"] = np.where(
        matched["tumor_dna_vaf"] >= CLONAL_VAF_THRESHOLD,
        "VAF >= 0.25",
        "VAF < 0.25",
    )

    return merged, matched


def make_summary(matched):
    all_som = matched.copy()
    top_som = matched[matched["is_top10_overall_per_sample"]].copy()
    non_top_som = matched[~matched["is_top10_overall_per_sample"]].copy()

    all_q1, all_q3 = safe_iqr(all_som["tumor_dna_vaf"])
    top_q1, top_q3 = safe_iqr(top_som["tumor_dna_vaf"])
    non_q1, non_q3 = safe_iqr(non_top_som["tumor_dna_vaf"])

    mw_p = safe_mannwhitney(top_som["tumor_dna_vaf"], non_top_som["tumor_dna_vaf"])

    n_all, rho_all, p_all = safe_spearman(
        all_som["tumor_dna_vaf"],
        all_som["rank_percentile_within_sample"],
    )

    rows = [{
        "analysis": "all_somatic_vs_top10_overall",
        "n_all_somatic_matched_vaf": int(len(all_som)),
        "median_vaf_all_somatic": float(np.median(all_som["tumor_dna_vaf"])) if len(all_som) else np.nan,
        "q1_vaf_all_somatic": all_q1,
        "q3_vaf_all_somatic": all_q3,
        "n_top10_overall_somatic_matched_vaf": int(len(top_som)),
        "median_vaf_top10_overall_somatic": float(np.median(top_som["tumor_dna_vaf"])) if len(top_som) else np.nan,
        "q1_vaf_top10_overall_somatic": top_q1,
        "q3_vaf_top10_overall_somatic": top_q3,
        "n_non_top10_somatic_matched_vaf": int(len(non_top_som)),
        "median_vaf_non_top10_somatic": float(np.median(non_top_som["tumor_dna_vaf"])) if len(non_top_som) else np.nan,
        "q1_vaf_non_top10_somatic": non_q1,
        "q3_vaf_non_top10_somatic": non_q3,
        "mannwhitney_top10_vs_non_top10_p": mw_p,
        "spearman_n_vaf_vs_rank_percentile": n_all,
        "spearman_rho_vaf_vs_rank_percentile": rho_all,
        "spearman_p_vaf_vs_rank_percentile": p_all,
        "n_vaf_ge_0_25_all_somatic": int((all_som["tumor_dna_vaf"] >= CLONAL_VAF_THRESHOLD).sum()),
        "pct_vaf_ge_0_25_all_somatic": 100 * (all_som["tumor_dna_vaf"] >= CLONAL_VAF_THRESHOLD).mean(),
        "n_vaf_ge_0_25_top10_somatic": int((top_som["tumor_dna_vaf"] >= CLONAL_VAF_THRESHOLD).sum()),
        "pct_vaf_ge_0_25_top10_somatic": 100 * (top_som["tumor_dna_vaf"] >= CLONAL_VAF_THRESHOLD).mean() if len(top_som) else np.nan,
    }]

    for mhc_class, g in matched.groupby("MHC_Class_clean", dropna=False):
        n, rho, p = safe_spearman(
            g["tumor_dna_vaf"],
            g["rank_percentile_within_sample"],
        )

        top = g[g["is_top10_overall_per_sample"]].copy()
        non_top = g[~g["is_top10_overall_per_sample"]].copy()

        rows.append({
            "analysis": f"{mhc_class}_vaf_vs_overall_rank_percentile",
            "n_all_somatic_matched_vaf": int(len(g)),
            "median_vaf_all_somatic": float(np.median(g["tumor_dna_vaf"])) if len(g) else np.nan,
            "q1_vaf_all_somatic": safe_iqr(g["tumor_dna_vaf"])[0],
            "q3_vaf_all_somatic": safe_iqr(g["tumor_dna_vaf"])[1],
            "n_top10_overall_somatic_matched_vaf": int(len(top)),
            "median_vaf_top10_overall_somatic": float(np.median(top["tumor_dna_vaf"])) if len(top) else np.nan,
            "q1_vaf_top10_overall_somatic": safe_iqr(top["tumor_dna_vaf"])[0],
            "q3_vaf_top10_overall_somatic": safe_iqr(top["tumor_dna_vaf"])[1],
            "n_non_top10_somatic_matched_vaf": int(len(non_top)),
            "median_vaf_non_top10_somatic": float(np.median(non_top["tumor_dna_vaf"])) if len(non_top) else np.nan,
            "q1_vaf_non_top10_somatic": safe_iqr(non_top["tumor_dna_vaf"])[0],
            "q3_vaf_non_top10_somatic": safe_iqr(non_top["tumor_dna_vaf"])[1],
            "mannwhitney_top10_vs_non_top10_p": safe_mannwhitney(top["tumor_dna_vaf"], non_top["tumor_dna_vaf"]),
            "spearman_n_vaf_vs_rank_percentile": n,
            "spearman_rho_vaf_vs_rank_percentile": rho,
            "spearman_p_vaf_vs_rank_percentile": p,
            "n_vaf_ge_0_25_all_somatic": int((g["tumor_dna_vaf"] >= CLONAL_VAF_THRESHOLD).sum()),
            "pct_vaf_ge_0_25_all_somatic": 100 * (g["tumor_dna_vaf"] >= CLONAL_VAF_THRESHOLD).mean(),
            "n_vaf_ge_0_25_top10_somatic": int((top["tumor_dna_vaf"] >= CLONAL_VAF_THRESHOLD).sum()),
            "pct_vaf_ge_0_25_top10_somatic": 100 * (top["tumor_dna_vaf"] >= CLONAL_VAF_THRESHOLD).mean() if len(top) else np.nan,
        })

    return pd.DataFrame(rows)


def add_trend_line(ax, x, y):
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
    ax.plot(xs, ys, color="black", linestyle="--", linewidth=1.2)


def draw_boxplot(matched, summary):
    set_pub_style()

    all_som = matched["tumor_dna_vaf"].dropna().values
    top_som = matched.loc[matched["is_top10_overall_per_sample"], "tumor_dna_vaf"].dropna().values

    fig = plt.figure(figsize=(5.0, 4.8))
    ax = plt.gca()

    bp = ax.boxplot(
        [all_som, top_som],
        labels=["All somatic", "Top 10 per sample\nsomatic subset"],
        patch_artist=True,
        showfliers=True,
    )

    colors = ["#1f77b4", "#ff7f0e"]

    for patch, color in zip(bp["boxes"], colors):
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
        flier.set_markeredgecolor("black")
        flier.set_alpha(0.55)

    ax.axhline(CLONAL_VAF_THRESHOLD, color="black", linestyle="--", linewidth=1.0)
    ax.text(
        0.02,
        CLONAL_VAF_THRESHOLD + 0.015,
        "VAF = 0.25",
        transform=ax.get_yaxis_transform(),
        ha="left",
        va="bottom",
        fontsize=9,
    )

    pval = summary.loc[
        summary["analysis"] == "all_somatic_vs_top10_overall",
        "mannwhitney_top10_vs_non_top10_p"
    ].iloc[0]

    ax.text(
        0.5,
        0.96,
        f"MW p = {pval:.3g}",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=10,
    )

    ax.set_ylabel("Tumor DNA VAF")
    ax.set_title("Tumor DNA VAF in somatic epitope candidates")
    ax.set_ylim(0, min(1.0, max(0.7, matched["tumor_dna_vaf"].max() + 0.08)))

    fig.tight_layout()

    png = PLOT_DIR / "A9_tumor_dna_vaf_all_somatic_vs_top10_complete_rescue.png"
    svg = PLOT_DIR / "A9_tumor_dna_vaf_all_somatic_vs_top10_complete_rescue.svg"

    fig.savefig(png, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)

    print("Wrote:")
    print(png)
    print(svg)


def draw_scatter(matched, summary):
    set_pub_style()

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.6), sharex=True, sharey=True)

    colors = {
        "MHC I": "#1f77b4",
        "MHC II": "#ff7f0e",
    }

    for ax, mhc_class in zip(axes, ["MHC I", "MHC II"]):
        sub = matched[matched["MHC_Class_clean"] == mhc_class].copy()

        x = sub["tumor_dna_vaf"].astype(float).values
        y = sub["rank_percentile_within_sample"].astype(float).values

        ax.scatter(
            x,
            y,
            s=18,
            alpha=0.35,
            color=colors[mhc_class],
            edgecolor="none",
        )

        add_trend_line(ax, x, y)

        ax.axvline(CLONAL_VAF_THRESHOLD, color="black", linestyle=":", linewidth=1.0)

        ax.set_title(mhc_class)
        ax.set_xlabel("Tumor DNA VAF")
        ax.set_ylabel("Overall rank percentile\n(0 = best rank)")
        ax.set_xlim(0, min(1.0, max(0.7, matched["tumor_dna_vaf"].max() + 0.08)))
        ax.set_ylim(-0.02, 1.02)

        row = summary[summary["analysis"] == f"{mhc_class}_vaf_vs_overall_rank_percentile"]
        if len(row):
            rho = row["spearman_rho_vaf_vs_rank_percentile"].iloc[0]
            pval = row["spearman_p_vaf_vs_rank_percentile"].iloc[0]
            n = int(row["spearman_n_vaf_vs_rank_percentile"].iloc[0])

            ax.text(
                0.03,
                0.97,
                f"Spearman ρ = {rho:.3f}\np = {pval:.3g}\nn = {n}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=10,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="0.7", alpha=0.9),
            )

    fig.suptitle("Tumor DNA VAF versus final prioritisation rank", y=1.02)
    fig.tight_layout()

    png = PLOT_DIR / "A9_tumor_dna_vaf_vs_rank_percentile_by_MHC_complete_rescue.png"
    svg = PLOT_DIR / "A9_tumor_dna_vaf_vs_rank_percentile_by_MHC_complete_rescue.svg"

    fig.savefig(png, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)

    print("Wrote:")
    print(png)
    print(svg)


def main():
    lookup = read_pvacseq_vaf_lookup()
    all_final, somatic = read_final_somatic_tables()
    merged, matched = attach_vaf_to_final_somatic(somatic, lookup)

    match_summary = pd.DataFrame([{
        "n_final_somatic_rows": int(len(somatic)),
        "n_final_somatic_rows_with_vaf": int(matched.shape[0]),
        "n_final_somatic_rows_without_vaf": int(merged["tumor_dna_vaf"].isna().sum()),
        "pct_final_somatic_rows_with_vaf": 100 * matched.shape[0] / len(somatic) if len(somatic) else np.nan,
        "clonal_vaf_threshold": CLONAL_VAF_THRESHOLD,
    }])

    summary = make_summary(matched)

    all_final.to_csv(TABLE_DIR / "A9_all_final_epitopes_with_rank_percentile_complete_rescue.tsv", sep="\t", index=False)
    merged.to_csv(TABLE_DIR / "A9_final_somatic_epitopes_with_pvacseq_vaf_complete_rescue.tsv", sep="\t", index=False)
    matched.to_csv(TABLE_DIR / "A9_final_somatic_epitopes_matched_vaf_only_complete_rescue.tsv", sep="\t", index=False)
    match_summary.to_csv(TABLE_DIR / "A9_vaf_match_summary_complete_rescue.tsv", sep="\t", index=False)
    summary.to_csv(TABLE_DIR / "A9_tumor_dna_vaf_clonality_summary_complete_rescue.tsv", sep="\t", index=False)

    draw_boxplot(matched, summary)
    draw_scatter(matched, summary)

    print()
    print("VAF match summary:")
    print(match_summary.to_string(index=False))
    print()
    print("Clonality/VAF summary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
