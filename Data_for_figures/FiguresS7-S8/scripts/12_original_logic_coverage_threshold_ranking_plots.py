#!/usr/bin/env python3

from pathlib import Path
import math
import re
import tarfile
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(".").resolve()

CANDIDATES = ROOT / "egg_candidates.exact_relaxed_evidence_flags.tsv.gz"
MATCHES = ROOT / "reviewer2_coverage_sensitivity" / "patched_builder_outputs" / "all_matches.exact_relaxed.with_coverage.tsv.gz"

OLD_SUMMARY = ROOT / "ranking_vs_ic50_exact_relaxed_outputs" / "summary_ranking_delta_vs_ic50.tsv"

OUTDIR = ROOT / "reviewer2_coverage_sensitivity" / "original_logic_coverage_threshold_plots"
FIGDIR = OUTDIR / "figures"
OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

TOP_PCTS = [10, 20, 30, 40, 50, 60]
TOP_NS = [10, 20, 30, 40, 50]

ENDPOINT = "strict_positive_experimental"
ENDPOINT_LABEL = "Strict positive experimental evidence"

STRICT_STRENGTHS = {
    "strong_immunogenicity_and_presentation",
    "strong_immunogenicity",
    "presentation",
}

THRESHOLDS = [
    ("relaxed_min8_peptide", "Relaxed min-8 peptide match", None),
    ("relaxed_min8_peptide_cov50", "Relaxed min-8 + ≥50% longer-peptide coverage", 0.50),
    ("relaxed_min8_peptide_cov70", "Relaxed min-8 + ≥70% longer-peptide coverage", 0.70),
    ("relaxed_min8_peptide_cov80", "Relaxed min-8 + ≥80% longer-peptide coverage", 0.80),
]

THRESHOLD_LABELS = {
    scope: label for scope, label, _ in THRESHOLDS
}

DATABASE_ORDER = ["Any", "CEDAR", "IEDB", "TSNAdb_v2"]

DATABASE_LABELS = {
    "Any": "Any database",
    "CEDAR": "CEDAR",
    "IEDB": "IEDB",
    "TSNAdb_v2": "TSNAdb v2.0",
}

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


def norm_colname(x):
    return re.sub(r"[^a-z0-9]+", "_", str(x).strip().lower()).strip("_")


def first_existing(df, candidates):
    norm_to_real = {norm_colname(c): c for c in df.columns}
    for c in candidates:
        key = norm_colname(c)
        if key in norm_to_real:
            return norm_to_real[key]
    return None


def bool_series(df, cols):
    if not cols:
        return pd.Series(False, index=df.index)

    out = pd.Series(False, index=df.index)

    for c in cols:
        s = df[c]

        if s.dtype == bool:
            b = s.fillna(False)
        else:
            b = (
                s.astype(str)
                .str.strip()
                .str.lower()
                .isin(["1", "true", "yes", "y"])
            )

            numeric = pd.to_numeric(s, errors="coerce").fillna(0)
            b = b | (numeric > 0)

        out = out | b

    return out


def scope_matches(col_lower, scope):
    has_exact = "exact" in col_lower
    has_relaxed = ("relaxed" in col_lower) or ("min8" in col_lower) or ("contain" in col_lower)
    has_hla = "hla" in col_lower
    has_peptide = "peptide" in col_lower

    looks_old_exact_peptide = (
        "__peptide_match" in col_lower
        and not has_relaxed
        and not has_hla
    )

    looks_old_exact_peptide_hla = (
        "__peptide_hla_match" in col_lower
        and not has_relaxed
        and has_hla
    )

    if scope == "exact_peptide":
        return has_peptide and not has_hla and (has_exact or looks_old_exact_peptide)

    if scope == "exact_peptide_hla":
        return has_peptide and has_hla and (has_exact or looks_old_exact_peptide_hla)

    if scope == "relaxed_min8_peptide":
        return has_peptide and not has_hla and has_relaxed

    if scope == "relaxed_min8_peptide_hla":
        return has_peptide and has_hla and has_relaxed

    return False


def evidence_columns(df, evidence_tokens, scope):
    cols = []

    for c in df.columns:
        cl = c.lower()

        if "negative" in cl:
            continue

        if not scope_matches(cl, scope):
            continue

        if any(tok in cl for tok in evidence_tokens):
            cols.append(c)

    return cols


def make_original_min8_endpoint(df):
    """
    Recreate the original strict-positive relaxed-min8 endpoint exactly from the
    original candidate-level evidence flag columns.
    """

    cols = evidence_columns(df, STRICT_STRENGTHS, "relaxed_min8_peptide")
    endpoint_col = f"endpoint__{ENDPOINT}__relaxed_min8_peptide"

    df[endpoint_col] = bool_series(df, cols).astype(int)

    inventory = pd.DataFrame([
        {
            "endpoint": ENDPOINT,
            "endpoint_label": ENDPOINT_LABEL,
            "match_scope": "relaxed_min8_peptide",
            "match_scope_label": "Relaxed min-8 peptide",
            "n_columns_used": len(cols),
            "columns_used": ";".join(cols),
            "n_positive_candidates": int(df[endpoint_col].sum()),
            "positive_fraction": float(df[endpoint_col].mean()) if len(df) else np.nan,
        }
    ])

    return df, inventory


def prepare_strict_matches(matches):
    m = matches.copy()

    required = [
        "candidate_id",
        "match_scope",
        "database",
        "evidence_strength",
        "overlap_len",
        "longer_peptide_coverage",
    ]

    missing = [c for c in required if c not in m.columns]
    if missing:
        raise ValueError(f"Missing required columns in match table: {missing}")

    m["candidate_id"] = m["candidate_id"].astype(str)
    m["overlap_len"] = pd.to_numeric(m["overlap_len"], errors="coerce")
    m["longer_peptide_coverage"] = pd.to_numeric(m["longer_peptide_coverage"], errors="coerce")

    m = m[
        m["match_scope"].eq("relaxed_min8_peptide")
        & m["evidence_strength"].isin(STRICT_STRENGTHS)
        & (m["overlap_len"] >= 8)
    ].copy()

    return m


def add_coverage_endpoint_flags(df, strict_matches):
    df = df.copy()
    df["candidate_id"] = df["candidate_id"].astype(str)

    for scope, label, threshold in THRESHOLDS:
        endpoint_col = f"endpoint__{ENDPOINT}__{scope}"

        if threshold is None:
            # Important: keep original min8 endpoint from candidate-level flags,
            # not rebuilt from the match table. This guarantees the old S6 values match.
            continue

        ids = set(
            strict_matches.loc[
                strict_matches["longer_peptide_coverage"] >= threshold,
                "candidate_id",
            ].astype(str)
        )

        df[endpoint_col] = df["candidate_id"].isin(ids).astype(int)

    return df


def build_database_count_table(df, strict_matches):
    all_candidate_ids = set(df["candidate_id"].astype(str))

    original_min8_col = f"endpoint__{ENDPOINT}__relaxed_min8_peptide"
    original_min8_ids = set(df.loc[df[original_min8_col] == 1, "candidate_id"].astype(str))

    rows = []

    for database in DATABASE_ORDER:
        if database == "Any":
            db_matches = strict_matches.copy()
        else:
            db_matches = strict_matches[strict_matches["database"].eq(database)].copy()

        for scope, label, threshold in THRESHOLDS:
            if threshold is None:
                if database == "Any":
                    ids = original_min8_ids
                else:
                    ids = set(db_matches["candidate_id"].astype(str)) & all_candidate_ids
            else:
                ids = set(
                    db_matches.loc[
                        db_matches["longer_peptide_coverage"] >= threshold,
                        "candidate_id",
                    ].astype(str)
                ) & all_candidate_ids

            rows.append({
                "database": database,
                "database_label": DATABASE_LABELS[database],
                "match_scope": scope,
                "criterion_label": label,
                "supported_candidates": len(ids),
                "total_candidates": len(df),
                "supported_fraction": len(ids) / len(df) if len(df) else np.nan,
            })

    counts = pd.DataFrame(rows)
    counts.to_csv(OUTDIR / "supported_candidate_counts_by_database_and_coverage_threshold.tsv", sep="\t", index=False)

    return counts


def build_ranking_methods(df):
    methods = []

    ic50_col = first_existing(df, ["Median.MT.IC50.Score", "Median_MT_IC50_Score"])
    if ic50_col is None:
        raise ValueError("Could not find Median.MT.IC50.Score column.")

    methods.append(("ic50", "IC50", ic50_col, "lower_better"))

    borda_rank = first_existing(df, ["Borda_Rank", "borda_rank_numeric"])
    borda_score = first_existing(df, ["Borda_Score"])

    if borda_rank is not None:
        methods.append(("borda", "Borda", borda_rank, "lower_better"))
    elif borda_score is not None:
        methods.append(("borda", "Borda", borda_score, "higher_better"))

    method_specs = [
        ("depmap_survivability", "DepMap survivability", ["Depmap_survivability_score", "DepMap_survivability_score"], "lower_better"),
        ("net_betweenness", "Network betweenness", ["Net_Betweenness", "Network_Betweenness"], "higher_better"),
        ("net_degree", "Network degree", ["Net_Degree", "Network_Degree"], "higher_better"),
        ("net_impact", "Network impact", ["Net_Impact", "Network_Impact"], "higher_better"),
        ("net_strength", "Network strength", ["Net_Strength", "Network_Strength"], "higher_better"),
        ("net_wci", "Network WCI", ["Net_WCI", "Network_WCI"], "higher_better"),
    ]

    vaf_col = first_existing(df, ["Tumor_DNA_VAF", "Tumor DNA VAF", "DNA_VAF", "tumor_vaf"])
    if vaf_col is not None:
        method_specs.append(("tumor_dna_vaf", "Tumor DNA VAF", [vaf_col], "higher_better"))

    for method_id, label, candidates, direction in method_specs:
        col = first_existing(df, candidates)
        if col is not None:
            methods.append((method_id, label, col, direction))

    pd.DataFrame(
        [
            {
                "ranking_method": m,
                "ranking_label": lab,
                "ranking_column": col,
                "direction": direction,
            }
            for m, lab, col, direction in methods
        ]
    ).to_csv(OUTDIR / "ranking_methods_used.tsv", sep="\t", index=False)

    return methods


def top_indices_for_group(g, col, direction, top_type, cutoff):
    """
    Original top-set selection logic from 02_compare...py:
    sort by raw feature/rank value, push missing values to bottom,
    and use candidate_id as deterministic tie-breaker.
    """

    values = pd.to_numeric(g[col], errors="coerce")
    tmp = g.copy()
    tmp["_rank_value"] = values

    if direction == "lower_better":
        tmp["_rank_value_sort"] = tmp["_rank_value"].fillna(np.inf)
        ascending = True
    else:
        tmp["_rank_value_sort"] = tmp["_rank_value"].fillna(-np.inf)
        ascending = False

    tmp = tmp.sort_values(
        ["_rank_value_sort", "candidate_id"],
        ascending=[ascending, True],
        kind="mergesort",
    )

    n = len(tmp)

    if top_type == "percent":
        k = max(1, int(math.ceil(n * cutoff / 100.0)))
    elif top_type == "topn":
        k = min(int(cutoff), n)
    else:
        raise ValueError(top_type)

    return tmp.index[:k]


def compute_per_patient(df, methods):
    sample_col = first_existing(df, ["sample_id", "Sample", "sample"])
    if sample_col is None:
        raise ValueError("Could not find sample_id column.")

    method_map = {m[0]: m for m in methods}
    ic50 = method_map["ic50"]

    endpoint_cols = [
        f"endpoint__{ENDPOINT}__{scope}"
        for scope, _, _ in THRESHOLDS
    ]

    missing = [c for c in endpoint_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing endpoint columns: {missing}")

    rows = []

    for endpoint_col in endpoint_cols:
        _, endpoint, scope = endpoint_col.split("__", 2)

        for top_type, cutoffs in [("percent", TOP_PCTS), ("topn", TOP_NS)]:
            for cutoff in cutoffs:
                for sample_id, g in df.groupby(sample_col, sort=True):
                    ic50_top = top_indices_for_group(g, ic50[2], ic50[3], top_type, cutoff)
                    ic50_hits = int(df.loc[ic50_top, endpoint_col].sum())
                    ic50_rate = float(df.loc[ic50_top, endpoint_col].mean()) if len(ic50_top) else np.nan

                    for method_id, label, col, direction in methods:
                        if method_id == "ic50":
                            continue

                        method_top = top_indices_for_group(g, col, direction, top_type, cutoff)
                        method_hits = int(df.loc[method_top, endpoint_col].sum())
                        method_rate = float(df.loc[method_top, endpoint_col].mean()) if len(method_top) else np.nan

                        rows.append({
                            "endpoint": endpoint,
                            "endpoint_label": ENDPOINT_LABEL,
                            "match_scope": scope,
                            "match_scope_label": THRESHOLD_LABELS.get(scope, scope),
                            "top_type": top_type,
                            "cutoff": cutoff,
                            "sample_id": sample_id,
                            "n_candidates_sample": len(g),
                            "ranking_method": method_id,
                            "ranking_label": label,
                            "ranking_column": col,
                            "ranking_direction": direction,
                            "method_top_n_used": len(method_top),
                            "ic50_top_n_used": len(ic50_top),
                            "method_hits": method_hits,
                            "ic50_hits": ic50_hits,
                            "delta_hits_vs_ic50": method_hits - ic50_hits,
                            "method_hit_rate": method_rate,
                            "ic50_hit_rate": ic50_rate,
                            "delta_hit_rate_vs_ic50": method_rate - ic50_rate,
                        })

    per_patient = pd.DataFrame(rows)
    per_patient.to_csv(OUTDIR / "per_patient_ranking_delta_vs_ic50_coverage_thresholds.tsv", sep="\t", index=False)

    return per_patient


def summarize(per_patient):
    def better_count(x):
        return int((x > 1e-12).sum())

    def worse_count(x):
        return int((x < -1e-12).sum())

    def same_count(x):
        return int((x.abs() <= 1e-12).sum())

    summary = (
        per_patient
        .groupby(
            [
                "endpoint",
                "endpoint_label",
                "match_scope",
                "match_scope_label",
                "top_type",
                "cutoff",
                "ranking_method",
                "ranking_label",
                "ranking_column",
                "ranking_direction",
            ],
            dropna=False,
        )
        .agg(
            n_patients=("sample_id", "nunique"),
            median_delta_hit_rate_vs_ic50=("delta_hit_rate_vs_ic50", "median"),
            mean_delta_hit_rate_vs_ic50=("delta_hit_rate_vs_ic50", "mean"),
            median_delta_hits_vs_ic50=("delta_hits_vs_ic50", "median"),
            sum_delta_hits_vs_ic50=("delta_hits_vs_ic50", "sum"),
            patients_better_than_ic50=("delta_hit_rate_vs_ic50", better_count),
            patients_same_as_ic50=("delta_hit_rate_vs_ic50", same_count),
            patients_worse_than_ic50=("delta_hit_rate_vs_ic50", worse_count),
            median_method_hit_rate=("method_hit_rate", "median"),
            median_ic50_hit_rate=("ic50_hit_rate", "median"),
        )
        .reset_index()
    )

    summary.to_csv(OUTDIR / "summary_ranking_delta_vs_ic50_coverage_thresholds.tsv", sep="\t", index=False)

    return summary


def savefig(path):
    plt.tight_layout()
    plt.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    plt.close()


def plot_previous_style_single_threshold(summary, scope, top_type, outstem):
    df = summary.copy()
    df = df[df["ranking_method"] != "tumor_dna_vaf"].copy()
    df = df[
        (df["endpoint"] == ENDPOINT)
        & (df["match_scope"] == scope)
        & (df["top_type"] == top_type)
    ].copy()

    if df.empty:
        print(f"[WARN] no data for {scope}, {top_type}")
        return

    df["median_delta_hit_rate_vs_ic50_pct_points"] = (
        100.0 * pd.to_numeric(df["median_delta_hit_rate_vs_ic50"], errors="coerce")
    )

    df[
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

    methods = [m for m in METHOD_ORDER if m in set(df["ranking_method"])]

    plt.figure(figsize=(9.5, 5.8))

    for method in methods:
        mdf = df[df["ranking_method"] == method].sort_values("cutoff")
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

    plt.xlabel(xlabel)
    plt.ylabel("Median delta hit rate vs IC50 ranking\npercentage points")
    plt.title(f"{ENDPOINT_LABEL}\n{THRESHOLD_LABELS[scope]}; {title_cutoff}")
    plt.grid(axis="y", alpha=0.25)
    plt.legend(loc="best", fontsize=8, frameon=True)

    savefig(FIGDIR / outstem)

    print(f"Wrote {FIGDIR / (outstem + '.png')}")
    print(f"Wrote {OUTDIR / (outstem + '.numbers.tsv')}")


def plot_combined_previous_style(summary):
    df = summary.copy()
    df = df[df["ranking_method"] != "tumor_dna_vaf"].copy()
    df["median_delta_hit_rate_vs_ic50_pct_points"] = (
        100.0 * pd.to_numeric(df["median_delta_hit_rate_vs_ic50"], errors="coerce")
    )

    fig, axes = plt.subplots(
        nrows=len(THRESHOLDS),
        ncols=2,
        figsize=(14, 4.0 * len(THRESHOLDS)),
        sharey=False,
    )

    for row_idx, (scope, label, _) in enumerate(THRESHOLDS):
        for col_idx, top_type in enumerate(["topn", "percent"]):
            ax = axes[row_idx, col_idx]

            sub = df[
                (df["endpoint"] == ENDPOINT)
                & (df["match_scope"] == scope)
                & (df["top_type"] == top_type)
            ].copy()

            for method in METHOD_ORDER:
                mdf = sub[sub["ranking_method"] == method].sort_values("cutoff")
                if mdf.empty:
                    continue

                x = pd.to_numeric(mdf["cutoff"], errors="coerce")
                y = pd.to_numeric(mdf["median_delta_hit_rate_vs_ic50_pct_points"], errors="coerce")

                ax.plot(
                    x,
                    y,
                    label=METHOD_LABELS.get(method, method),
                    marker=MARKERS.get(method, "o"),
                    linestyle=LINESTYLES.get(method, "-"),
                    linewidth=1.8,
                    markersize=5.5,
                    alpha=0.9,
                )

            ax.axhline(0, linestyle="--", linewidth=1.0)
            ax.grid(axis="y", alpha=0.25)

            if top_type == "topn":
                ax.set_xlabel("Top-N candidates within each patient")
                title_cutoff = "top-N cutoffs"
            else:
                ax.set_xlabel("Top-percent cutoff within each patient")
                title_cutoff = "top-percent cutoffs"

            ax.set_ylabel("Median delta hit rate vs IC50 ranking\npercentage points")
            ax.set_title(f"{label}; {title_cutoff}")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False, fontsize=9)

    fig.suptitle(
        "Strict positive experimental evidence recovery across coverage thresholds\nall databases pooled; original ranking-comparison logic",
        y=0.995,
        fontsize=14,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.955])

    fig.savefig(FIGDIR / "coverage_thresholds_all_methods_previous_style_original_logic.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR / "coverage_thresholds_all_methods_previous_style_original_logic.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / "coverage_thresholds_all_methods_previous_style_original_logic.svg", bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {FIGDIR / 'coverage_thresholds_all_methods_previous_style_original_logic.png'}")


def plot_count_barplot(counts):
    d = counts.copy()
    d["match_scope"] = pd.Categorical(
        d["match_scope"],
        categories=[scope for scope, _, _ in THRESHOLDS],
        ordered=True,
    )
    d["database"] = pd.Categorical(d["database"], categories=DATABASE_ORDER, ordered=True)
    d = d.sort_values(["match_scope", "database"])

    x = np.arange(len(THRESHOLDS))
    width = 0.20

    fig, ax = plt.subplots(figsize=(9.5, 5.2))

    for i, database in enumerate(DATABASE_ORDER):
        sub = d[d["database"].eq(database)].copy()

        values = [
            int(sub.loc[sub["match_scope"].eq(scope), "supported_candidates"].iloc[0])
            for scope, _, _ in THRESHOLDS
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
                fontsize=8,
            )

    ax.set_xticks(x)
    ax.set_xticklabels([label for _, label, _ in THRESHOLDS], rotation=15, ha="right")
    ax.set_ylabel("Supported EGG candidates")
    ax.set_xlabel("External-evidence matching criterion")
    ax.set_title("Supported EGG candidates retained under peptide-coverage thresholds")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(axis="y", alpha=0.25)

    fig.tight_layout()

    fig.savefig(FIGDIR / "supported_candidate_counts_by_database_and_coverage_threshold_original_logic.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR / "supported_candidate_counts_by_database_and_coverage_threshold_original_logic.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / "supported_candidate_counts_by_database_and_coverage_threshold_original_logic.svg", bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {FIGDIR / 'supported_candidate_counts_by_database_and_coverage_threshold_original_logic.png'}")


def verify_min8_against_old_summary(new_summary):
    if not OLD_SUMMARY.exists():
        print(f"[WARN] Old summary not found, skipping verification: {OLD_SUMMARY}")
        return

    old = pd.read_csv(OLD_SUMMARY, sep="\t")

    old_sub = old[
        (old["endpoint"] == ENDPOINT)
        & (old["match_scope"] == "relaxed_min8_peptide")
    ].copy()

    new_sub = new_summary[
        (new_summary["endpoint"] == ENDPOINT)
        & (new_summary["match_scope"] == "relaxed_min8_peptide")
    ].copy()

    keep = [
        "endpoint",
        "match_scope",
        "top_type",
        "cutoff",
        "ranking_method",
        "ranking_label",
        "median_delta_hit_rate_vs_ic50",
        "patients_better_than_ic50",
        "patients_same_as_ic50",
        "patients_worse_than_ic50",
    ]

    old_sub = old_sub[keep].rename(columns={
        "median_delta_hit_rate_vs_ic50": "old_median_delta_hit_rate_vs_ic50",
        "patients_better_than_ic50": "old_patients_better_than_ic50",
        "patients_same_as_ic50": "old_patients_same_as_ic50",
        "patients_worse_than_ic50": "old_patients_worse_than_ic50",
    })

    new_sub = new_sub[keep].rename(columns={
        "median_delta_hit_rate_vs_ic50": "new_median_delta_hit_rate_vs_ic50",
        "patients_better_than_ic50": "new_patients_better_than_ic50",
        "patients_same_as_ic50": "new_patients_same_as_ic50",
        "patients_worse_than_ic50": "new_patients_worse_than_ic50",
    })

    merged = old_sub.merge(
        new_sub,
        on=[
            "endpoint",
            "match_scope",
            "top_type",
            "cutoff",
            "ranking_method",
            "ranking_label",
        ],
        how="outer",
    )

    merged["delta_new_minus_old"] = (
        merged["new_median_delta_hit_rate_vs_ic50"]
        - merged["old_median_delta_hit_rate_vs_ic50"]
    )

    merged["delta_new_minus_old_pct_points"] = 100 * merged["delta_new_minus_old"]

    merged.to_csv(
        OUTDIR / "verification_original_min8_summary_vs_new_original_logic.tsv",
        sep="\t",
        index=False,
    )

    finite = merged["delta_new_minus_old"].dropna()

    if len(finite):
        max_abs = float(finite.abs().max())
    else:
        max_abs = np.nan

    print()
    print("Verification against old relaxed-min8 summary:")
    print(f"  Rows compared: {len(merged)}")
    print(f"  Max absolute delta difference: {max_abs}")

    mismatch = merged[
        merged["delta_new_minus_old"].notna()
        & (merged["delta_new_minus_old"].abs() > 1e-12)
    ].copy()

    if mismatch.empty:
        print("  PASS: old relaxed-min8 summary values are reproduced exactly.")
    else:
        print("  WARNING: mismatches found. See:")
        print(f"  {OUTDIR / 'verification_original_min8_summary_vs_new_original_logic.tsv'}")
        print(mismatch.to_string(index=False))


def main():
    if not CANDIDATES.exists():
        raise FileNotFoundError(f"Missing candidate file: {CANDIDATES}")

    if not MATCHES.exists():
        raise FileNotFoundError(f"Missing match file: {MATCHES}")

    print("Reading candidates:")
    print(CANDIDATES)
    df = pd.read_csv(CANDIDATES, sep="\t", low_memory=False)

    # Avoid accidentally reusing endpoint columns from a previous run.
    df = df.loc[:, [c for c in df.columns if not c.startswith("endpoint__")]].copy()

    if "candidate_id" not in df.columns:
        df = df.reset_index(drop=True)
        df["candidate_id"] = np.arange(1, len(df) + 1).astype(str)

    df["candidate_id"] = df["candidate_id"].astype(str)

    print("Reading full match table with coverage:")
    print(MATCHES)
    matches = pd.read_csv(MATCHES, sep="\t", low_memory=False)

    strict_matches = prepare_strict_matches(matches)

    print("Recreating original strict-positive relaxed-min8 endpoint from original evidence flags...")
    df, endpoint_inventory = make_original_min8_endpoint(df)

    print("Adding coverage-threshold endpoints from full match-level table...")
    df = add_coverage_endpoint_flags(df, strict_matches)

    print("Writing endpoint inventory and candidate support table...")
    endpoint_inventory.to_csv(OUTDIR / "endpoint_inventory_original_min8.tsv", sep="\t", index=False)

    support_cols = [
        "candidate_id",
        "sample_id",
    ] + [
        f"endpoint__{ENDPOINT}__{scope}"
        for scope, _, _ in THRESHOLDS
    ]

    support_cols = [c for c in support_cols if c in df.columns]

    df[support_cols].to_csv(
        OUTDIR / "candidate_support_flags_original_min8_plus_coverage_thresholds.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    counts = build_database_count_table(df, strict_matches)

    print("Building ranking methods using original directions...")
    methods = build_ranking_methods(df)

    print("Computing per-patient deltas using original top-set logic...")
    per_patient = compute_per_patient(df, methods)

    print("Summarising...")
    summary = summarize(per_patient)

    print("Verifying min8 results against old summary...")
    verify_min8_against_old_summary(summary)

    print("Plotting previous-style single-threshold plots...")
    for scope, _, _ in THRESHOLDS:
        for top_type in ["topn", "percent"]:
            outstem = f"{ENDPOINT}__{scope}__{top_type}_original_logic"
            plot_previous_style_single_threshold(summary, scope, top_type, outstem)

    print("Plotting combined previous-style coverage-threshold figure...")
    plot_combined_previous_style(summary)

    print("Plotting candidate-count bar plot...")
    plot_count_barplot(counts)

    tar_path = OUTDIR / "original_logic_coverage_threshold_plots.tar.gz"
    with tarfile.open(tar_path, "w:gz") as tar:
        tar.add(FIGDIR, arcname="figures")
        for f in [
            OUTDIR / "summary_ranking_delta_vs_ic50_coverage_thresholds.tsv",
            OUTDIR / "per_patient_ranking_delta_vs_ic50_coverage_thresholds.tsv",
            OUTDIR / "supported_candidate_counts_by_database_and_coverage_threshold.tsv",
            OUTDIR / "candidate_support_flags_original_min8_plus_coverage_thresholds.tsv.gz",
            OUTDIR / "verification_original_min8_summary_vs_new_original_logic.tsv",
            OUTDIR / "ranking_methods_used.tsv",
        ]:
            if f.exists():
                tar.add(f, arcname=f.name)

    print()
    print("Done.")
    print(f"Output directory: {OUTDIR.resolve()}")
    print(f"Figures: {FIGDIR.resolve()}")
    print(f"Package: {tar_path.resolve()}")
    print()
    print("Candidate counts:")
    print(counts.to_string(index=False))


if __name__ == "__main__":
    main()
