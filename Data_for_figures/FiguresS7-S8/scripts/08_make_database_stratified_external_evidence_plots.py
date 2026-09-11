#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path.cwd()

INPUT_DIR = BASE / "reviewer2_coverage_sensitivity" / "patched_builder_outputs"
OUT = BASE / "reviewer2_coverage_sensitivity" / "database_stratified_plots"
OUT.mkdir(parents=True, exist_ok=True)

CANDIDATES_IN = INPUT_DIR / "egg_candidates.exact_relaxed_evidence_flags.tsv.gz"
MATCHES_IN = INPUT_DIR / "all_matches.exact_relaxed.with_coverage.tsv.gz"

SUPPORTED_EVIDENCE_STRENGTHS = {
    "presentation",
    "strong_immunogenicity",
    "strong_immunogenicity_and_presentation",
}

# These match the previous-looking figure you showed.
TOP_N_CUTOFFS_PREVIOUS_STYLE = [10, 20, 30, 40, 50]
TOP_PERCENT_CUTOFFS_PREVIOUS_STYLE = [10, 20, 30, 40, 50, 60]

# These match the coverage-sensitivity figure style.
TOP_N_CUTOFFS_COVERAGE = [10, 20, 50]
TOP_PERCENT_CUTOFFS_COVERAGE = [1, 5, 10, 20]

DATABASE_ORDER = ["Any", "CEDAR", "IEDB", "TSNAdb_v2"]

DATABASE_LABELS = {
    "Any": "Any database",
    "CEDAR": "CEDAR",
    "IEDB": "IEDB",
    "TSNAdb_v2": "TSNAdb v2.0",
}

CRITERIA = ["min8", "min8_cov50", "min8_cov70", "min8_cov80"]

CRITERION_LABELS = {
    "min8": "Min-8 containment",
    "min8_cov50": "Min-8 + ≥50% coverage",
    "min8_cov70": "Min-8 + ≥70% coverage",
    "min8_cov80": "Min-8 + ≥80% coverage",
}

# Ranking methods compared with IC50-only.
# Lower rank/value is better when ascending=True.
# For network metrics, higher centrality is usually better, hence ascending=False.
# If your original plotting script used a different DepMap direction, change it here.
METHODS = {
    "Borda": {
        "column": "Borda_Rank",
        "ascending": True,
        "rank_column": "rank__Borda",
    },
    "DepMap survivability": {
        "column": "Depmap_survivability_score",
        "ascending": True,
        "rank_column": "rank__DepMap_survivability",
    },
    "Network betweenness": {
        "column": "Net_Betweenness",
        "ascending": False,
        "rank_column": "rank__Network_betweenness",
    },
    "Network degree": {
        "column": "Net_Degree",
        "ascending": False,
        "rank_column": "rank__Network_degree",
    },
    "Network impact": {
        "column": "Net_Impact",
        "ascending": False,
        "rank_column": "rank__Network_impact",
    },
    "Network strength": {
        "column": "Net_Strength",
        "ascending": False,
        "rank_column": "rank__Network_strength",
    },
    "Network WCI": {
        "column": "Net_WCI",
        "ascending": False,
        "rank_column": "rank__Network_WCI",
    },
}

IC50_RANK_COL = "rank__IC50_only"


def safe_numeric(x):
    return pd.to_numeric(x, errors="coerce")


def make_ranks(cand):
    cand = cand.copy()

    cand["Median.MT.IC50.Score"] = safe_numeric(cand["Median.MT.IC50.Score"])
    cand[IC50_RANK_COL] = (
        cand.groupby("sample_id")["Median.MT.IC50.Score"]
        .rank(method="first", ascending=True, na_option="bottom")
    )

    for method, cfg in METHODS.items():
        col = cfg["column"]
        rank_col = cfg["rank_column"]

        if col not in cand.columns:
            print(f"WARNING: missing column for {method}: {col}. Skipping this method.")
            cand[rank_col] = np.nan
            continue

        cand[col] = safe_numeric(cand[col])
        cand[rank_col] = (
            cand.groupby("sample_id")[col]
            .rank(method="first", ascending=cfg["ascending"], na_option="bottom")
        )

    return cand


def prepare_matches(matches):
    matches = matches.copy()

    matches["overlap_len"] = safe_numeric(matches["overlap_len"])
    matches["longer_peptide_coverage"] = safe_numeric(matches["longer_peptide_coverage"])

    matches = matches[
        matches["evidence_strength"].isin(SUPPORTED_EVIDENCE_STRENGTHS)
        & matches["match_scope"].eq("relaxed_min8_peptide")
        & (matches["overlap_len"] >= 8)
    ].copy()

    # Remove exact repeated evidence rows, but keep distinct candidate-reference peptide matches.
    dedup_cols = [
        "candidate_id",
        "database",
        "evidence_strength",
        "evidence_class",
        "reference_peptide",
        "reference_hla_normalized",
        "reference_gene",
    ]
    dedup_cols = [c for c in dedup_cols if c in matches.columns]
    matches = matches.drop_duplicates(subset=dedup_cols)

    return matches


def add_support_flags(cand, matches):
    cand = cand.copy()

    for database in DATABASE_ORDER:
        if database == "Any":
            db_matches = matches.copy()
        else:
            db_matches = matches[matches["database"].eq(database)].copy()

        supported_sets = {
            "min8": set(db_matches["candidate_id"]),
            "min8_cov50": set(db_matches.loc[db_matches["longer_peptide_coverage"] >= 0.50, "candidate_id"]),
            "min8_cov70": set(db_matches.loc[db_matches["longer_peptide_coverage"] >= 0.70, "candidate_id"]),
            "min8_cov80": set(db_matches.loc[db_matches["longer_peptide_coverage"] >= 0.80, "candidate_id"]),
        }

        for criterion, ids in supported_sets.items():
            flag_col = f"supported__{database}__{criterion}"
            cand[flag_col] = cand["candidate_id"].isin(ids)

    return cand


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


def evaluate_methods(cand, criteria, databases, top_n_cutoffs, top_percent_cutoffs):
    rows = []

    rank_map = {method: cfg["rank_column"] for method, cfg in METHODS.items()}
    rank_map["IC50-only"] = IC50_RANK_COL

    for sample_id, sdf in cand.groupby("sample_id"):
        sdf = sdf.copy()

        for database in databases:
            for criterion in criteria:
                hit_col = f"supported__{database}__{criterion}"

                for method, rank_col in rank_map.items():
                    for n in top_n_cutoffs:
                        rows.append({
                            "sample_id": sample_id,
                            "database": database,
                            "database_label": DATABASE_LABELS[database],
                            "criterion": criterion,
                            "criterion_label": CRITERION_LABELS[criterion],
                            "cutoff_type": "Top-N",
                            "cutoff": n,
                            "method": method,
                            "hit_rate": hit_rate_top_n(sdf, rank_col, hit_col, n),
                        })

                    for pct in top_percent_cutoffs:
                        rows.append({
                            "sample_id": sample_id,
                            "database": database,
                            "database_label": DATABASE_LABELS[database],
                            "criterion": criterion,
                            "criterion_label": CRITERION_LABELS[criterion],
                            "cutoff_type": "Top-percentile",
                            "cutoff": pct,
                            "method": method,
                            "hit_rate": hit_rate_top_percent(sdf, rank_col, hit_col, pct),
                        })

    return pd.DataFrame(rows)


def difference_vs_ic50(hit_rates):
    ic50 = hit_rates[hit_rates["method"].eq("IC50-only")].copy()
    ic50 = ic50.rename(columns={"hit_rate": "ic50_hit_rate"})
    ic50 = ic50.drop(columns=["method"])

    d = hit_rates.merge(
        ic50,
        on=[
            "sample_id",
            "database",
            "database_label",
            "criterion",
            "criterion_label",
            "cutoff_type",
            "cutoff",
        ],
        how="left",
    )

    d = d[~d["method"].eq("IC50-only")].copy()
    d["hit_rate_difference_vs_ic50"] = d["hit_rate"] - d["ic50_hit_rate"]
    d["hit_rate_difference_vs_ic50_percentage_points"] = 100 * d["hit_rate_difference_vs_ic50"]

    return d


def summarize_diff(diff):
    return (
        diff.groupby(
            [
                "database",
                "database_label",
                "criterion",
                "criterion_label",
                "cutoff_type",
                "cutoff",
                "method",
            ],
            as_index=False,
        )
        .agg(
            median_difference=("hit_rate_difference_vs_ic50", "median"),
            median_difference_percentage_points=("hit_rate_difference_vs_ic50_percentage_points", "median"),
            q1_difference_percentage_points=("hit_rate_difference_vs_ic50_percentage_points", lambda x: np.nanpercentile(x, 25)),
            q3_difference_percentage_points=("hit_rate_difference_vs_ic50_percentage_points", lambda x: np.nanpercentile(x, 75)),
            n_samples=("hit_rate_difference_vs_ic50", lambda x: x.notna().sum()),
        )
    )


def supported_counts(cand):
    rows = []

    for database in DATABASE_ORDER:
        for criterion in CRITERIA:
            flag_col = f"supported__{database}__{criterion}"
            rows.append({
                "database": database,
                "database_label": DATABASE_LABELS[database],
                "criterion": criterion,
                "criterion_label": CRITERION_LABELS[criterion],
                "supported_candidates": int(cand[flag_col].sum()),
                "total_candidates": len(cand),
                "supported_fraction": cand[flag_col].mean(),
            })

    return pd.DataFrame(rows)


def plot_previous_style_database_stratified(summary):
    """
    Similar to the previous figure:
    rows = database
    columns = Top-N and Top-percentile
    lines = ranking methods
    endpoint = relaxed min-8 only
    """
    d = summary[summary["criterion"].eq("min8")].copy()

    fig, axes = plt.subplots(
        nrows=len(DATABASE_ORDER),
        ncols=2,
        figsize=(14, 4.2 * len(DATABASE_ORDER)),
        sharey=False,
    )

    for row_idx, database in enumerate(DATABASE_ORDER):
        for col_idx, cutoff_type in enumerate(["Top-N", "Top-percentile"]):
            ax = axes[row_idx, col_idx]
            sub = d[
                d["database"].eq(database)
                & d["cutoff_type"].eq(cutoff_type)
            ].copy()

            for method in METHODS:
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
            ax.set_title(f"{DATABASE_LABELS[database]}: {cutoff_type}")
            ax.set_xlabel(
                "Top-N candidates within each patient"
                if cutoff_type == "Top-N"
                else "Top-percent cutoff within each patient"
            )
            ax.set_ylabel("Median delta hit rate vs IC50 ranking\npercentage points")
            ax.grid(True, alpha=0.25)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False)

    fig.suptitle(
        "Positive experimental evidence, min-8 peptide match, stratified by database",
        y=0.995,
        fontsize=14,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.965])

    fig.savefig(OUT / "database_stratified_min8_all_methods_previous_style.pdf")
    fig.savefig(OUT / "database_stratified_min8_all_methods_previous_style.png", dpi=300)
    fig.savefig(OUT / "database_stratified_min8_all_methods_previous_style.svg")
    plt.close(fig)


def plot_coverage_database_stratified(summary, counts):
    """
    Like Supplementary Figure S7 but stratified by database.
    rows = database
    columns = Top-N and Top-percentile
    lines = coverage thresholds
    method = Borda only
    """
    d = summary[summary["method"].eq("Borda")].copy()

    fig, axes = plt.subplots(
        nrows=len(DATABASE_ORDER),
        ncols=2,
        figsize=(14, 4.2 * len(DATABASE_ORDER)),
        sharey=False,
    )

    for row_idx, database in enumerate(DATABASE_ORDER):
        for col_idx, cutoff_type in enumerate(["Top-N", "Top-percentile"]):
            ax = axes[row_idx, col_idx]
            sub = d[
                d["database"].eq(database)
                & d["cutoff_type"].eq(cutoff_type)
            ].copy()

            for criterion in CRITERIA:
                c = sub[sub["criterion"].eq(criterion)].sort_values("cutoff")
                if c.empty:
                    continue

                ax.plot(
                    c["cutoff"],
                    c["median_difference_percentage_points"],
                    marker="o",
                    label=CRITERION_LABELS[criterion],
                )

            ax.axhline(0, linestyle="--", linewidth=1)
            ax.set_title(f"{DATABASE_LABELS[database]}: {cutoff_type}")
            ax.set_xlabel(
                "Top-N candidates within each patient"
                if cutoff_type == "Top-N"
                else "Top-percent cutoff within each patient"
            )
            ax.set_ylabel("Median delta hit rate vs IC50 ranking\npercentage points")
            ax.grid(True, alpha=0.25)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False)

    fig.suptitle(
        "Coverage-threshold sensitivity of external evidence recovery, stratified by database",
        y=0.995,
        fontsize=14,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.965])

    fig.savefig(OUT / "database_stratified_coverage_sensitivity_borda_vs_ic50.pdf")
    fig.savefig(OUT / "database_stratified_coverage_sensitivity_borda_vs_ic50.png", dpi=300)
    fig.savefig(OUT / "database_stratified_coverage_sensitivity_borda_vs_ic50.svg")
    plt.close(fig)

    # Bar plot of supported candidates by database and criterion.
    fig, axes = plt.subplots(
        nrows=1,
        ncols=len(DATABASE_ORDER),
        figsize=(4.2 * len(DATABASE_ORDER), 4.5),
        sharey=True,
    )

    for idx, database in enumerate(DATABASE_ORDER):
        ax = axes[idx]
        sub = counts[counts["database"].eq(database)].copy()
        sub["criterion"] = pd.Categorical(sub["criterion"], categories=CRITERIA, ordered=True)
        sub = sub.sort_values("criterion")

        ax.bar(sub["criterion_label"], sub["supported_candidates"])
        ax.set_title(DATABASE_LABELS[database])
        ax.set_xlabel("Matching criterion")
        ax.tick_params(axis="x", rotation=35)
        if idx == 0:
            ax.set_ylabel("Supported EGG candidates")

    fig.suptitle("Supported candidate counts by database and coverage threshold", y=1.02)
    fig.tight_layout()

    fig.savefig(OUT / "database_stratified_supported_candidate_counts.pdf")
    fig.savefig(OUT / "database_stratified_supported_candidate_counts.png", dpi=300)
    fig.savefig(OUT / "database_stratified_supported_candidate_counts.svg")
    plt.close(fig)


def main():
    if not CANDIDATES_IN.exists():
        raise FileNotFoundError(f"Missing candidate file: {CANDIDATES_IN}")
    if not MATCHES_IN.exists():
        raise FileNotFoundError(f"Missing match file: {MATCHES_IN}")

    cand = pd.read_csv(CANDIDATES_IN, sep="\t", compression="gzip")
    matches = pd.read_csv(MATCHES_IN, sep="\t", compression="gzip")

    cand = make_ranks(cand)
    matches = prepare_matches(matches)
    cand = add_support_flags(cand, matches)

    counts = supported_counts(cand)
    counts.to_csv(
        OUT / "database_stratified_supported_candidate_counts.tsv",
        sep="\t",
        index=False,
    )

    # Previous-style all-methods comparison.
    hit_rates_prev = evaluate_methods(
        cand=cand,
        criteria=["min8"],
        databases=DATABASE_ORDER,
        top_n_cutoffs=TOP_N_CUTOFFS_PREVIOUS_STYLE,
        top_percent_cutoffs=TOP_PERCENT_CUTOFFS_PREVIOUS_STYLE,
    )

    diff_prev = difference_vs_ic50(hit_rates_prev)
    summary_prev = summarize_diff(diff_prev)

    hit_rates_prev.to_csv(
        OUT / "database_stratified_min8_all_methods_hit_rates.tsv",
        sep="\t",
        index=False,
    )
    diff_prev.to_csv(
        OUT / "database_stratified_min8_all_methods_difference_vs_ic50.tsv",
        sep="\t",
        index=False,
    )
    summary_prev.to_csv(
        OUT / "database_stratified_min8_all_methods_summary.tsv",
        sep="\t",
        index=False,
    )

    plot_previous_style_database_stratified(summary_prev)

    # Coverage-sensitivity comparison, Borda and all methods calculated but Borda plotted.
    hit_rates_cov = evaluate_methods(
        cand=cand,
        criteria=CRITERIA,
        databases=DATABASE_ORDER,
        top_n_cutoffs=TOP_N_CUTOFFS_COVERAGE,
        top_percent_cutoffs=TOP_PERCENT_CUTOFFS_COVERAGE,
    )

    diff_cov = difference_vs_ic50(hit_rates_cov)
    summary_cov = summarize_diff(diff_cov)

    hit_rates_cov.to_csv(
        OUT / "database_stratified_coverage_all_methods_hit_rates.tsv",
        sep="\t",
        index=False,
    )
    diff_cov.to_csv(
        OUT / "database_stratified_coverage_all_methods_difference_vs_ic50.tsv",
        sep="\t",
        index=False,
    )
    summary_cov.to_csv(
        OUT / "database_stratified_coverage_all_methods_summary.tsv",
        sep="\t",
        index=False,
    )

    plot_coverage_database_stratified(summary_cov, counts)

    # Save final candidate table with support flags for audit.
    cand.to_csv(
        OUT / "egg_candidates.database_stratified_support_flags.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    print("Done.")
    print(f"Output directory: {OUT.resolve()}")
    print()
    print("Supported candidate counts:")
    print(counts.to_string(index=False))


if __name__ == "__main__":
    main()
