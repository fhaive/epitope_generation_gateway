#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path.cwd()

INPUT_DIR = BASE / "reviewer2_coverage_sensitivity" / "patched_builder_outputs"
OUT = BASE / "reviewer2_coverage_sensitivity" / "final_s7_counts_and_borda_effect"
OUT.mkdir(parents=True, exist_ok=True)

CANDIDATES_IN = INPUT_DIR / "egg_candidates.exact_relaxed_evidence_flags.tsv.gz"
MATCHES_IN = INPUT_DIR / "all_matches.exact_relaxed.with_coverage.tsv.gz"

SUPPORTED_EVIDENCE_STRENGTHS = {
    "presentation",
    "strong_immunogenicity",
    "strong_immunogenicity_and_presentation",
}

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
    "min8_cov50": "Min-8 + ≥50%",
    "min8_cov70": "Min-8 + ≥70%",
    "min8_cov80": "Min-8 + ≥80%",
}

TOP_N_CUTOFFS = [10, 20, 30, 40, 50]
TOP_PERCENT_CUTOFFS = [10, 20, 30, 40, 50, 60]


def safe_numeric(x):
    return pd.to_numeric(x, errors="coerce")


def prepare_matches(matches):
    matches = matches.copy()

    matches["overlap_len"] = safe_numeric(matches["overlap_len"])
    matches["longer_peptide_coverage"] = safe_numeric(matches["longer_peptide_coverage"])

    matches = matches[
        matches["evidence_strength"].isin(SUPPORTED_EVIDENCE_STRENGTHS)
        & matches["match_scope"].eq("relaxed_min8_peptide")
        & (matches["overlap_len"] >= 8)
    ].copy()

    # Remove exact repeated evidence rows, while preserving distinct reference-peptide matches.
    # This prevents row-level duplication but does not discard alternative peptide matches
    # that may have different coverage values.
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
            cand[f"supported__{database}__{criterion}"] = cand["candidate_id"].isin(ids)

    return cand


def make_ranks(cand):
    cand = cand.copy()

    cand["Borda_Rank"] = safe_numeric(cand["Borda_Rank"])
    cand["Median.MT.IC50.Score"] = safe_numeric(cand["Median.MT.IC50.Score"])

    cand["IC50_Rank"] = (
        cand
        .groupby("sample_id")["Median.MT.IC50.Score"]
        .rank(method="first", ascending=True, na_option="bottom")
    )

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


def evaluate_borda_vs_ic50(cand):
    rows = []

    for sample_id, sdf in cand.groupby("sample_id"):
        sdf = sdf.copy()

        for criterion in CRITERIA:
            hit_col = f"supported__Any__{criterion}"

            for n in TOP_N_CUTOFFS:
                borda_hr = hit_rate_top_n(sdf, "Borda_Rank", hit_col, n)
                ic50_hr = hit_rate_top_n(sdf, "IC50_Rank", hit_col, n)

                rows.append({
                    "sample_id": sample_id,
                    "criterion": criterion,
                    "criterion_label": CRITERION_LABELS[criterion],
                    "cutoff_type": "Top-N",
                    "cutoff": n,
                    "borda_hit_rate": borda_hr,
                    "ic50_hit_rate": ic50_hr,
                    "difference": borda_hr - ic50_hr,
                    "difference_percentage_points": 100 * (borda_hr - ic50_hr),
                })

            for pct in TOP_PERCENT_CUTOFFS:
                borda_hr = hit_rate_top_percent(sdf, "Borda_Rank", hit_col, pct)
                ic50_hr = hit_rate_top_percent(sdf, "IC50_Rank", hit_col, pct)

                rows.append({
                    "sample_id": sample_id,
                    "criterion": criterion,
                    "criterion_label": CRITERION_LABELS[criterion],
                    "cutoff_type": "Top-percentile",
                    "cutoff": pct,
                    "borda_hit_rate": borda_hr,
                    "ic50_hit_rate": ic50_hr,
                    "difference": borda_hr - ic50_hr,
                    "difference_percentage_points": 100 * (borda_hr - ic50_hr),
                })

    return pd.DataFrame(rows)


def summarize_borda_effect(df):
    return (
        df
        .groupby(["criterion", "criterion_label", "cutoff_type", "cutoff"], as_index=False)
        .agg(
            median_difference=("difference", "median"),
            median_difference_percentage_points=("difference_percentage_points", "median"),
            q1_difference_percentage_points=("difference_percentage_points", lambda x: np.nanpercentile(x, 25)),
            q3_difference_percentage_points=("difference_percentage_points", lambda x: np.nanpercentile(x, 75)),
            n_samples=("difference", lambda x: x.notna().sum()),
        )
    )


def make_counts_table(cand):
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


def plot_counts_panel(ax, counts):
    counts = counts.copy()
    counts["criterion"] = pd.Categorical(counts["criterion"], categories=CRITERIA, ordered=True)
    counts["database"] = pd.Categorical(counts["database"], categories=DATABASE_ORDER, ordered=True)
    counts = counts.sort_values(["criterion", "database"])

    x = np.arange(len(CRITERIA))
    width = 0.20

    for i, database in enumerate(DATABASE_ORDER):
        sub = counts[counts["database"].eq(database)].copy()
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
                fontsize=7,
            )

    ax.set_xticks(x)
    ax.set_xticklabels([CRITERION_LABELS[c] for c in CRITERIA], rotation=20, ha="right")
    ax.set_ylabel("Supported EGG candidates")
    ax.set_xlabel("Matching criterion")
    ax.set_title("A. Candidate support by database")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)


def plot_borda_panel(ax, summary, cutoff_type, title):
    d = summary[summary["cutoff_type"].eq(cutoff_type)].copy()

    for criterion in CRITERIA:
        sub = d[d["criterion"].eq(criterion)].sort_values("cutoff")
        ax.plot(
            sub["cutoff"],
            sub["median_difference_percentage_points"],
            marker="o",
            label=CRITERION_LABELS[criterion],
        )

    ax.axhline(0, linestyle="--", linewidth=1)
    ax.grid(True, alpha=0.25)
    ax.set_title(title)
    ax.set_xlabel(
        "Top-N candidates within each patient"
        if cutoff_type == "Top-N"
        else "Top-percent cutoff within each patient"
    )
    ax.set_ylabel("Median delta hit rate vs IC50 ranking\npercentage points")
    ax.legend(frameon=False, fontsize=8)


def main():
    cand = pd.read_csv(CANDIDATES_IN, sep="\t", compression="gzip")
    matches = pd.read_csv(MATCHES_IN, sep="\t", compression="gzip")

    matches = prepare_matches(matches)
    cand = make_ranks(cand)
    cand = add_support_flags(cand, matches)

    counts = make_counts_table(cand)
    borda_effect = evaluate_borda_vs_ic50(cand)
    borda_summary = summarize_borda_effect(borda_effect)

    counts.to_csv(OUT / "s7_supported_candidate_counts_by_database_and_coverage_threshold.tsv", sep="\t", index=False)
    borda_effect.to_csv(OUT / "s7_borda_vs_ic50_by_sample.tsv", sep="\t", index=False)
    borda_summary.to_csv(OUT / "s7_borda_vs_ic50_summary.tsv", sep="\t", index=False)

    # Separate count-only bar plot.
    fig, ax = plt.subplots(figsize=(9, 5))
    plot_counts_panel(ax, counts)
    fig.tight_layout()
    fig.savefig(OUT / "S7A_supported_candidate_counts_by_database_and_coverage_threshold.pdf")
    fig.savefig(OUT / "S7A_supported_candidate_counts_by_database_and_coverage_threshold.png", dpi=300)
    fig.savefig(OUT / "S7A_supported_candidate_counts_by_database_and_coverage_threshold.svg")
    plt.close(fig)

    # Separate Borda-effect plot, previous-style two-panel layout.
    fig, axes = plt.subplots(2, 1, figsize=(9, 9), sharey=False)
    plot_borda_panel(axes[0], borda_summary, "Top-N", "A. Borda recovery versus IC50-only; top-N cutoffs")
    plot_borda_panel(axes[1], borda_summary, "Top-percentile", "B. Borda recovery versus IC50-only; top-percentile cutoffs")
    fig.suptitle("Effect of coverage thresholds on external evidence recovery; all databases pooled", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT / "S7B_C_borda_vs_ic50_coverage_thresholds_previous_style.pdf")
    fig.savefig(OUT / "S7B_C_borda_vs_ic50_coverage_thresholds_previous_style.png", dpi=300)
    fig.savefig(OUT / "S7B_C_borda_vs_ic50_coverage_thresholds_previous_style.svg")
    plt.close(fig)

    # Combined final S7: counts + Borda effect.
    fig, axes = plt.subplots(3, 1, figsize=(9, 13))

    plot_counts_panel(axes[0], counts)
    plot_borda_panel(
        axes[1],
        borda_summary,
        "Top-N",
        "B. Borda recovery versus IC50-only; top-N cutoffs"
    )
    plot_borda_panel(
        axes[2],
        borda_summary,
        "Top-percentile",
        "C. Borda recovery versus IC50-only; top-percentile cutoffs"
    )

    fig.suptitle(
        "Sensitivity of external evidence support to peptide-length-normalised coverage thresholds",
        y=0.995,
        fontsize=13,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.98])

    fig.savefig(OUT / "Supplementary_Figure_S7_counts_and_borda_coverage_effect.pdf")
    fig.savefig(OUT / "Supplementary_Figure_S7_counts_and_borda_coverage_effect.png", dpi=300)
    fig.savefig(OUT / "Supplementary_Figure_S7_counts_and_borda_coverage_effect.svg")
    plt.close(fig)

    print("Done.")
    print(f"Output directory: {OUT.resolve()}")
    print()
    print("Candidate counts:")
    print(counts.to_string(index=False))
    print()
    print("Borda summary:")
    print(borda_summary.to_string(index=False))


if __name__ == "__main__":
    main()
