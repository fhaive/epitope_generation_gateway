#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


BASE = Path("TCGA_melanoma")

FINAL_DIR_WITH_VAF = BASE / "epitopes_prioritisation_complete_rescue_with_vaf" / "final_epitopes"
FINAL_DIR_DEFAULT = BASE / "epitopes_prioritisation_complete_rescue" / "final_epitopes"

OUT_DIR = BASE / "rescue_final_analysis" / "reviewer_updates" / "reviewer_A5_ppi_topology_ablation"
PLOT_DIR = OUT_DIR / "plots"
TABLE_DIR = OUT_DIR / "tables"

PLOT_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)

TOP_N = 10

IC50_CANDIDATES = [
    "Median.MT.IC50.Score",
    "Median MT IC50 Score",
    "Median_MT_IC50_Score",
]

DEPMAP_CANDIDATES = [
    "Depmap_survivability_score",
    "DepMap_survivability_score",
    "Depmap.Survivability.Score",
]

WCI_CANDIDATES = [
    "Net_WCI",
    "WCI",
]

RANK_CANDIDATES = [
    "Borda_Rank",
    "Borda Rank",
    "Borda.Rank",
    "Final_Rank",
    "Rank",
]

SCORE_CANDIDATES = [
    "Borda_Score",
    "Borda Score",
    "Borda.Score",
]

PPI_TOPOLOGY_COLS = [
    "Net_Betweenness",
    "Net_Degree",
    "Net_Impact",
    "Net_Strength",
]

IDENTITY_CANDIDATES = [
    "Gene.Name",
    "Gene Name",
    "Gene",
    "HLA.Allele",
    "HLA Allele",
    "MT.Epitope.Seq",
    "MT Epitope Seq",
    "Mutation.Source",
    "Mutation Source",
    "HLA.Class",
    "HLA Class",
]


def resolve_final_dir():
    if FINAL_DIR_WITH_VAF.exists():
        return FINAL_DIR_WITH_VAF
    if FINAL_DIR_DEFAULT.exists():
        return FINAL_DIR_DEFAULT
    raise FileNotFoundError(
        f"Could not find final epitope tables in either:\n"
        f"{FINAL_DIR_WITH_VAF}\n{FINAL_DIR_DEFAULT}"
    )


def find_col(df, candidates, required=True):
    for c in candidates:
        if c in df.columns:
            return c
    if required:
        raise ValueError(
            f"Could not find any of columns: {candidates}\n"
            f"Available columns:\n" + "\n".join(df.columns.astype(str))
        )
    return None


def sample_from_file(path):
    return path.name.replace("_epitopes_final.csv", "")


def within_sample_score(x, direction):
    """
    Return score where larger is better.
    direction='low': lower raw value is better.
    direction='high': higher raw value is better.
    """
    x = pd.to_numeric(x, errors="coerce")
    out = pd.Series(np.nan, index=x.index, dtype=float)
    valid = x.notna()

    if valid.sum() == 0:
        return out

    if valid.sum() == 1:
        out.loc[x[valid].index] = 1.0
        return out

    if direction == "low":
        ranks = x.rank(method="average", ascending=False, na_option="keep")
    elif direction == "high":
        ranks = x.rank(method="average", ascending=True, na_option="keep")
    else:
        raise ValueError(direction)

    return (ranks - 1) / (valid.sum() - 1)


def compute_ppi_topology_ablated_rank(df, ic50_col, depmap_col, wci_col):
    """
    PPI/topology ablation ranking.

    Removes:
      Net_Betweenness
      Net_Degree
      Net_Impact / LCI
      Net_Strength

    Retains:
      MHC-IC50 = 0.25
      DepMap = 0.25
      WCI = 0.50

    Rationale:
      WCI is retained as a patient-specific co-expression/network perturbation metric,
      while PPI-derived centrality/topology metrics are ablated.
    """
    d = df.copy()

    d["score_mhc_ic50"] = within_sample_score(d[ic50_col], direction="low")
    d["score_depmap"] = within_sample_score(d[depmap_col], direction="low")
    d["score_wci"] = within_sample_score(d[wci_col], direction="high")

    for c in ["score_mhc_ic50", "score_depmap", "score_wci"]:
        d[c] = d[c].fillna(d[c].median())
        d[c] = d[c].fillna(0.5)

    d["ppi_topology_ablated_score"] = (
        0.25 * d["score_mhc_ic50"]
        + 0.25 * d["score_depmap"]
        + 0.50 * d["score_wci"]
    )

    d["ppi_topology_ablated_rank"] = d["ppi_topology_ablated_score"].rank(
        method="first",
        ascending=False,
    )

    return d


def plot_histogram(series, bins, xlabel, ylabel, title, outbase):
    fig, ax = plt.subplots(figsize=(6.3, 5.0), dpi=300)

    ax.hist(series.dropna(), bins=bins, edgecolor="black", linewidth=0.7)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(axis="y", linestyle="--", alpha=0.35)

    fig.tight_layout()

    for ext in ["png", "svg", "pdf"]:
        fig.savefig(outbase.with_suffix(f".{ext}"), bbox_inches="tight")

    plt.close(fig)


def plot_scatter(all_comp):
    valid = all_comp[["original_Borda_Rank", "ppi_topology_ablated_rank"]].dropna().copy()

    fig, ax = plt.subplots(figsize=(6.2, 6.0), dpi=300)

    ax.scatter(
        valid["original_Borda_Rank"],
        valid["ppi_topology_ablated_rank"],
        alpha=0.38,
        s=14,
        linewidths=0,
    )

    mx = max(valid["original_Borda_Rank"].max(), valid["ppi_topology_ablated_rank"].max())
    ax.plot([1, mx], [1, mx], linestyle="--", color="black", linewidth=1.0)

    ax.set_xlabel("Original Borda rank")
    ax.set_ylabel("Rank after removing PPI/topology metrics")
    ax.set_title("Original versus PPI/topology-ablated ranking")
    ax.grid(True, linestyle="--", alpha=0.35)

    fig.tight_layout()

    outbase = PLOT_DIR / "A5_original_vs_ppi_topology_ablated_rank_scatter"
    for ext in ["png", "svg", "pdf"]:
        fig.savefig(outbase.with_suffix(f".{ext}"), bbox_inches="tight")

    plt.close(fig)


def main():
    final_dir = resolve_final_dir()
    files = sorted(final_dir.glob("*_epitopes_final.csv"))

    if not files:
        raise FileNotFoundError(f"No *_epitopes_final.csv files found in {final_dir}")

    summary_rows = []
    all_rows = []

    for f in files:
        sample = sample_from_file(f)
        df = pd.read_csv(f, low_memory=False)
        df["sample"] = sample
        df["row_index_within_sample"] = np.arange(len(df))

        rank_col = find_col(df, RANK_CANDIDATES)
        score_col = find_col(df, SCORE_CANDIDATES, required=False)
        ic50_col = find_col(df, IC50_CANDIDATES)
        depmap_col = find_col(df, DEPMAP_CANDIDATES)
        wci_col = find_col(df, WCI_CANDIDATES)

        d = compute_ppi_topology_ablated_rank(df, ic50_col, depmap_col, wci_col)

        d["original_Borda_Rank"] = pd.to_numeric(d[rank_col], errors="coerce")
        d["original_Borda_Score"] = pd.to_numeric(d[score_col], errors="coerce") if score_col else np.nan

        d["rank_shift_ppi_topology_ablated_minus_original"] = (
            d["ppi_topology_ablated_rank"] - d["original_Borda_Rank"]
        )
        d["abs_rank_shift"] = d["rank_shift_ppi_topology_ablated_minus_original"].abs()

        old_top = set(d.nsmallest(TOP_N, "original_Borda_Rank").index)
        new_top = set(d.nsmallest(TOP_N, "ppi_topology_ablated_rank").index)

        overlap = len(old_top & new_top)
        union = len(old_top | new_top)
        jaccard = overlap / union if union else np.nan

        valid = d[["original_Borda_Rank", "ppi_topology_ablated_rank"]].dropna()

        if len(valid) >= 3 and valid["original_Borda_Rank"].nunique() > 1 and valid["ppi_topology_ablated_rank"].nunique() > 1:
            sp = spearmanr(valid["original_Borda_Rank"], valid["ppi_topology_ablated_rank"])
            rho = float(sp.statistic)
            pval = float(sp.pvalue)
        else:
            rho, pval = np.nan, np.nan

        summary_rows.append({
            "sample": sample,
            "n_epitopes": len(d),
            "spearman_rho_original_vs_ppi_topology_ablated_rank": rho,
            "spearman_p_original_vs_ppi_topology_ablated_rank": pval,
            "top10_overlap_n": overlap,
            "top10_overlap_fraction": overlap / TOP_N,
            "top10_jaccard": jaccard,
            "median_abs_rank_shift": d["abs_rank_shift"].median(),
            "mean_abs_rank_shift": d["abs_rank_shift"].mean(),
            "max_abs_rank_shift": d["abs_rank_shift"].max(),
            "ic50_col": ic50_col,
            "depmap_col": depmap_col,
            "wci_col": wci_col,
            "rank_col": rank_col,
            "score_col": score_col if score_col else "",
        })

        keep_cols = [
            "sample",
            "row_index_within_sample",
            "original_Borda_Rank",
            "original_Borda_Score",
            "ppi_topology_ablated_rank",
            "ppi_topology_ablated_score",
            "score_mhc_ic50",
            "score_depmap",
            "score_wci",
            "rank_shift_ppi_topology_ablated_minus_original",
            "abs_rank_shift",
        ]

        identity_cols = [c for c in IDENTITY_CANDIDATES if c in d.columns]
        ppi_cols_present = [c for c in PPI_TOPOLOGY_COLS if c in d.columns]

        all_rows.append(d[keep_cols + identity_cols + ppi_cols_present].copy())

    summary = pd.DataFrame(summary_rows).sort_values("sample")
    all_comp = pd.concat(all_rows, ignore_index=True)

    summary.to_csv(TABLE_DIR / "A5_ppi_topology_ablation_summary_by_sample.tsv", sep="\t", index=False)
    all_comp.to_csv(TABLE_DIR / "A5_ppi_topology_ablation_detailed_ranks.tsv.gz", sep="\t", index=False, compression="gzip")

    cohort = pd.DataFrame([{
        "n_samples": summary["sample"].nunique(),
        "median_spearman_rho_original_vs_ppi_topology_ablated_rank": summary["spearman_rho_original_vs_ppi_topology_ablated_rank"].median(),
        "q1_spearman_rho_original_vs_ppi_topology_ablated_rank": summary["spearman_rho_original_vs_ppi_topology_ablated_rank"].quantile(0.25),
        "q3_spearman_rho_original_vs_ppi_topology_ablated_rank": summary["spearman_rho_original_vs_ppi_topology_ablated_rank"].quantile(0.75),
        "median_top10_overlap_n": summary["top10_overlap_n"].median(),
        "q1_top10_overlap_n": summary["top10_overlap_n"].quantile(0.25),
        "q3_top10_overlap_n": summary["top10_overlap_n"].quantile(0.75),
        "median_top10_overlap_fraction": summary["top10_overlap_fraction"].median(),
        "median_abs_rank_shift": summary["median_abs_rank_shift"].median(),
        "analysis_note": "GO/Hallmark/pathway annotations are not used in Borda ranking. PPI/topology ablation removes Net_Betweenness, Net_Degree, Net_Impact/LCI, and Net_Strength, while retaining MHC-IC50, DepMap, and WCI. WCI receives the total weight of the removed PPI/topology metrics.",
    }])

    cohort.to_csv(TABLE_DIR / "A5_ppi_topology_ablation_cohort_summary.tsv", sep="\t", index=False)

    pd.DataFrame([{
        "annotation_layer": "GO/Hallmark/pathway annotations",
        "used_in_Borda_ranking": "No",
        "expected_rank_effect_if_removed": "None",
        "interpretation": "These annotations are retained for biological interpretation only; removing them does not alter Borda scores or ranks.",
    }]).to_csv(TABLE_DIR / "A5_GO_pathway_annotation_rank_impact_note.tsv", sep="\t", index=False)

    pd.DataFrame([
        {"scenario": "Original/default", "MHC_IC50_weight": 0.25, "DepMap_weight": 0.25, "WCI_weight": 0.10, "PPI_topology_total_weight": 0.40, "PPI_topology_metrics": "Net_Betweenness;Net_Degree;Net_Impact/LCI;Net_Strength"},
        {"scenario": "PPI/topology ablated", "MHC_IC50_weight": 0.25, "DepMap_weight": 0.25, "WCI_weight": 0.50, "PPI_topology_total_weight": 0.00, "PPI_topology_metrics": "removed"},
    ]).to_csv(TABLE_DIR / "A5_ppi_topology_ablation_weighting_scheme.tsv", sep="\t", index=False)

    plot_histogram(
        summary["top10_overlap_n"],
        bins=np.arange(-0.5, TOP_N + 1.5, 1),
        xlabel="Top-10 overlap count",
        ylabel="Number of samples",
        title="Top-10 overlap after removing PPI/topology metrics",
        outbase=PLOT_DIR / "A5_top10_overlap_ppi_topology_ablated_histogram",
    )

    plot_histogram(
        summary["spearman_rho_original_vs_ppi_topology_ablated_rank"],
        bins=np.linspace(0, 1, 13),
        xlabel="Spearman ρ: original rank vs PPI/topology-ablated rank",
        ylabel="Number of samples",
        title="Full-rank correlation after removing PPI/topology metrics",
        outbase=PLOT_DIR / "A5_spearman_ppi_topology_ablated_histogram",
    )

    plot_scatter(all_comp)

    print("Wrote A5 outputs to:")
    print(OUT_DIR)
    print()
    print("Cohort summary:")
    print(cohort.to_string(index=False))


if __name__ == "__main__":
    main()
