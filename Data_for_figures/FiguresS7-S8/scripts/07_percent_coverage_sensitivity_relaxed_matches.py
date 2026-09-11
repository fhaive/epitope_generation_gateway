#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path.cwd()
INPUT_DIR = BASE / "reviewer2_coverage_sensitivity" / "patched_builder_outputs"
OUT = BASE / "reviewer2_coverage_sensitivity" / "coverage_sensitivity_outputs"
OUT.mkdir(parents=True, exist_ok=True)

CANDIDATES_IN = INPUT_DIR / "egg_candidates.exact_relaxed_evidence_flags.tsv.gz"
MATCHES_IN = INPUT_DIR / "all_matches.exact_relaxed.with_coverage.tsv.gz"

SUPPORTED_EVIDENCE_STRENGTHS = {
    "presentation",
    "strong_immunogenicity",
    "strong_immunogenicity_and_presentation",
}

TOP_N_CUTOFFS = [10, 20, 50]
TOP_PERCENT_CUTOFFS = [1, 5, 10, 20]

CRITERIA = [
    "min8",
    "min8_cov50",
    "min8_cov70",
    "min8_cov80",
]

CRITERION_LABELS = {
    "min8": "Min-8 containment",
    "min8_cov50": "Min-8 + >=50% coverage",
    "min8_cov70": "Min-8 + >=70% coverage",
    "min8_cov80": "Min-8 + >=80% coverage",
}


def safe_numeric(x):
    return pd.to_numeric(x, errors="coerce")


def hit_rate_top_n(df, rank_col, hit_col, n):
    d = df.dropna(subset=[rank_col]).sort_values(rank_col, ascending=True)
    if len(d) == 0:
        return np.nan
    k = min(n, len(d))
    return d.head(k)[hit_col].mean()


def hit_rate_top_percent(df, rank_col, hit_col, percent):
    d = df.dropna(subset=[rank_col]).sort_values(rank_col, ascending=True)
    if len(d) == 0:
        return np.nan
    k = max(1, int(np.ceil(len(d) * percent / 100.0)))
    return d.head(k)[hit_col].mean()


def evaluate(cand, criteria):
    rows = []

    for sample_id, sdf in cand.groupby("sample_id"):
        sdf = sdf.copy()

        for criterion in criteria:
            hit_col = f"supported_{criterion}"

            for method, rank_col in {
                "Borda": "Borda_Rank",
                "IC50-only": "IC50_Rank",
            }.items():

                for n in TOP_N_CUTOFFS:
                    rows.append({
                        "sample_id": sample_id,
                        "criterion": criterion,
                        "criterion_label": CRITERION_LABELS[criterion],
                        "cutoff_type": "Top-N",
                        "cutoff": n,
                        "method": method,
                        "hit_rate": hit_rate_top_n(sdf, rank_col, hit_col, n),
                    })

                for pct in TOP_PERCENT_CUTOFFS:
                    rows.append({
                        "sample_id": sample_id,
                        "criterion": criterion,
                        "criterion_label": CRITERION_LABELS[criterion],
                        "cutoff_type": "Top-percentile",
                        "cutoff": pct,
                        "method": method,
                        "hit_rate": hit_rate_top_percent(sdf, rank_col, hit_col, pct),
                    })

    return pd.DataFrame(rows)


def difference_vs_ic50(hit_rates):
    ic50 = hit_rates[hit_rates["method"] == "IC50-only"].copy()
    ic50 = ic50.rename(columns={"hit_rate": "ic50_hit_rate"})
    ic50 = ic50.drop(columns=["method"])

    d = hit_rates.merge(
        ic50,
        on=["sample_id", "criterion", "criterion_label", "cutoff_type", "cutoff"],
        how="left",
    )

    d = d[d["method"] != "IC50-only"].copy()
    d["hit_rate_difference_vs_ic50"] = d["hit_rate"] - d["ic50_hit_rate"]
    return d


def summarize(diff):
    return (
        diff
        .groupby(["criterion", "criterion_label", "cutoff_type", "cutoff", "method"], as_index=False)
        .agg(
            median_difference=("hit_rate_difference_vs_ic50", "median"),
            q1_difference=("hit_rate_difference_vs_ic50", lambda x: np.nanpercentile(x, 25)),
            q3_difference=("hit_rate_difference_vs_ic50", lambda x: np.nanpercentile(x, 75)),
            n_samples=("hit_rate_difference_vs_ic50", lambda x: x.notna().sum()),
        )
    )


def plot_line(ax, summary, cutoff_type, title):
    d = summary[
        (summary["cutoff_type"] == cutoff_type)
        & (summary["method"] == "Borda")
    ].copy()

    for criterion in CRITERIA:
        sub = d[d["criterion"] == criterion].sort_values("cutoff")
        ax.plot(
            sub["cutoff"],
            sub["median_difference"],
            marker="o",
            label=CRITERION_LABELS[criterion],
        )

    ax.axhline(0, linestyle="--", linewidth=1)
    ax.set_xlabel("Top-N cutoff" if cutoff_type == "Top-N" else "Top-percentile cutoff")
    ax.set_ylabel("Median hit-rate difference vs IC50-only")
    ax.set_title(title)


def main():
    if not CANDIDATES_IN.exists():
        raise FileNotFoundError(f"Missing candidate file: {CANDIDATES_IN}")
    if not MATCHES_IN.exists():
        raise FileNotFoundError(f"Missing all-matches file: {MATCHES_IN}")

    cand = pd.read_csv(CANDIDATES_IN, sep="\t", compression="gzip")
    matches = pd.read_csv(MATCHES_IN, sep="\t", compression="gzip")

    cand["Borda_Rank"] = safe_numeric(cand["Borda_Rank"])
    cand["Median.MT.IC50.Score"] = safe_numeric(cand["Median.MT.IC50.Score"])

    # IC50-only ranking within each patient/sample. Lower IC50 is better.
    cand["IC50_Rank"] = (
        cand
        .groupby("sample_id")["Median.MT.IC50.Score"]
        .rank(method="first", ascending=True)
    )

    matches["longer_peptide_coverage"] = safe_numeric(matches["longer_peptide_coverage"])
    matches["overlap_len"] = safe_numeric(matches["overlap_len"])

    # Keep only evidence categories used for the experimentally supported endpoint.
    # Binding-only, prediction-only, negative, and uncategorized evidence are excluded.
    matches_supported = matches[
        matches["evidence_strength"].isin(SUPPORTED_EVIDENCE_STRENGTHS)
    ].copy()

    relaxed = matches_supported[
        matches_supported["match_scope"].eq("relaxed_min8_peptide")
        & (matches_supported["overlap_len"] >= 8)
    ].copy()

    supported_sets = {
        "min8": set(relaxed["candidate_id"]),
        "min8_cov50": set(relaxed.loc[relaxed["longer_peptide_coverage"] >= 0.50, "candidate_id"]),
        "min8_cov70": set(relaxed.loc[relaxed["longer_peptide_coverage"] >= 0.70, "candidate_id"]),
        "min8_cov80": set(relaxed.loc[relaxed["longer_peptide_coverage"] >= 0.80, "candidate_id"]),
    }

    for criterion, ids in supported_sets.items():
        cand[f"supported_{criterion}"] = cand["candidate_id"].isin(ids)

    # Diagnostic: compare rebuilt min-8 support with original candidate-level flags.
    original_min8_cols = [
        c for c in cand.columns
        if c.startswith("relaxed_min8_peptide__")
        and c.split("__")[-1] in SUPPORTED_EVIDENCE_STRENGTHS
    ]

    if original_min8_cols:
        cand["supported_min8_original_flags"] = cand[original_min8_cols].astype(bool).any(axis=1)
        mismatches = cand[cand["supported_min8"] != cand["supported_min8_original_flags"]].copy()
    else:
        cand["supported_min8_original_flags"] = np.nan
        mismatches = pd.DataFrame()

    mismatches.to_csv(
        OUT / "diagnostic_min8_rebuilt_vs_original_flag_mismatches.tsv",
        sep="\t",
        index=False
    )

    diagnostic = pd.DataFrame([
        {
            "metric": "n_candidates",
            "value": len(cand),
        },
        {
            "metric": "n_pairwise_matches_all_scopes_all_evidence",
            "value": len(matches),
        },
        {
            "metric": "n_pairwise_relaxed_min8_supported_evidence_matches",
            "value": len(relaxed),
        },
        {
            "metric": "n_original_min8_supported_columns_detected",
            "value": len(original_min8_cols),
        },
        {
            "metric": "n_min8_rebuilt_vs_original_mismatches",
            "value": len(mismatches),
        },
    ])
    diagnostic.to_csv(OUT / "diagnostic_summary.tsv", sep="\t", index=False)

    cand.to_csv(
        OUT / "egg_candidates.coverage_sensitivity_evidence_flags.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    counts = []
    for criterion in CRITERIA:
        counts.append({
            "criterion": criterion,
            "criterion_label": CRITERION_LABELS[criterion],
            "supported_candidates": int(cand[f"supported_{criterion}"].sum()),
            "total_candidates": len(cand),
            "supported_fraction": cand[f"supported_{criterion}"].mean(),
        })

    counts = pd.DataFrame(counts)
    counts.to_csv(
        OUT / "supported_candidate_counts_by_coverage_threshold.tsv",
        sep="\t",
        index=False
    )

    hit_rates = evaluate(cand, CRITERIA)
    hit_rates.to_csv(
        OUT / "hit_rates_by_sample_method_coverage_threshold.tsv",
        sep="\t",
        index=False
    )

    diff = difference_vs_ic50(hit_rates)
    diff.to_csv(
        OUT / "hit_rate_difference_vs_ic50_by_sample_coverage_threshold.tsv",
        sep="\t",
        index=False
    )

    summary = summarize(diff)
    summary.to_csv(
        OUT / "summary_hit_rate_difference_vs_ic50_coverage_threshold.tsv",
        sep="\t",
        index=False
    )

    relaxed.to_csv(
        OUT / "relaxed_min8_supported_matches_used_for_coverage_sensitivity.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    # Individual plots.
    fig, ax = plt.subplots(figsize=(7.5, 5))
    plot_line(ax, summary, "Top-N", "Coverage sensitivity across top-N cutoffs")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "coverage_sensitivity_topN_borda_vs_ic50.pdf")
    fig.savefig(OUT / "coverage_sensitivity_topN_borda_vs_ic50.png", dpi=300)
    fig.savefig(OUT / "coverage_sensitivity_topN_borda_vs_ic50.svg")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 5))
    plot_line(ax, summary, "Top-percentile", "Coverage sensitivity across top-percentile cutoffs")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "coverage_sensitivity_topPercent_borda_vs_ic50.pdf")
    fig.savefig(OUT / "coverage_sensitivity_topPercent_borda_vs_ic50.png", dpi=300)
    fig.savefig(OUT / "coverage_sensitivity_topPercent_borda_vs_ic50.svg")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.bar(counts["criterion_label"], counts["supported_candidates"])
    ax.set_ylabel("Supported EGG candidates")
    ax.set_xlabel("Matching criterion")
    ax.set_title("External-evidence support under coverage thresholds")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(OUT / "supported_candidate_counts_by_coverage_threshold.pdf")
    fig.savefig(OUT / "supported_candidate_counts_by_coverage_threshold.png", dpi=300)
    fig.savefig(OUT / "supported_candidate_counts_by_coverage_threshold.svg")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    vals = relaxed["longer_peptide_coverage"].dropna()
    ax.hist(vals, bins=np.linspace(0, 1, 21), edgecolor="black")
    ax.set_xlabel("Matched fraction of longer peptide")
    ax.set_ylabel("Pairwise relaxed min-8 matches")
    ax.set_title("Coverage distribution among relaxed min-8 matches")
    fig.tight_layout()
    fig.savefig(OUT / "longer_peptide_coverage_distribution_relaxed_min8.pdf")
    fig.savefig(OUT / "longer_peptide_coverage_distribution_relaxed_min8.png", dpi=300)
    fig.savefig(OUT / "longer_peptide_coverage_distribution_relaxed_min8.svg")
    plt.close(fig)

    # Combined 4-panel follow-up figure.
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    plot_line(
        axes[0, 0],
        summary,
        "Top-N",
        "A. Top-N recovery versus IC50-only"
    )

    plot_line(
        axes[0, 1],
        summary,
        "Top-percentile",
        "B. Top-percentile recovery versus IC50-only"
    )

    axes[1, 0].bar(counts["criterion_label"], counts["supported_candidates"])
    axes[1, 0].set_ylabel("Supported EGG candidates")
    axes[1, 0].set_xlabel("Matching criterion")
    axes[1, 0].set_title("C. Supported candidates by criterion")
    axes[1, 0].tick_params(axis="x", rotation=30)

    axes[1, 1].hist(vals, bins=np.linspace(0, 1, 21), edgecolor="black")
    axes[1, 1].set_xlabel("Matched fraction of longer peptide")
    axes[1, 1].set_ylabel("Pairwise relaxed min-8 matches")
    axes[1, 1].set_title("D. Longer-peptide coverage distribution")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False)

    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(OUT / "Supplementary_Figure_S7_coverage_sensitivity_combined.pdf")
    fig.savefig(OUT / "Supplementary_Figure_S7_coverage_sensitivity_combined.png", dpi=300)
    fig.savefig(OUT / "Supplementary_Figure_S7_coverage_sensitivity_combined.svg")
    plt.close(fig)

    print("Done.")
    print(f"Output directory: {OUT.resolve()}")
    print()
    print("Diagnostic summary:")
    print(diagnostic.to_string(index=False))
    print()
    print("Supported candidate counts:")
    print(counts.to_string(index=False))


if __name__ == "__main__":
    main()
