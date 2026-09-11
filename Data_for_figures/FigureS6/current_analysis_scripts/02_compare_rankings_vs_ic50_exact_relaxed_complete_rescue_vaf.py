#!/usr/bin/env python3

from pathlib import Path
import math
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(".").resolve()
INFILE = ROOT / "egg_candidates.exact_relaxed_evidence_flags.tsv.gz"
OUTDIR = ROOT / "ranking_vs_ic50_exact_relaxed_outputs"
FIGDIR = OUTDIR / "figures"
OUTDIR.mkdir(exist_ok=True)
FIGDIR.mkdir(exist_ok=True)

TOP_PCTS = [10, 20, 30, 40, 50, 60]
TOP_NS = [10, 20, 30, 40, 50]

# Evidence groups to benchmark.
# Strict = direct positive experimental/validated evidence.
# Broad = strict evidence plus prediction-concordance and binding evidence.
ENDPOINTS = {
    "strict_positive_experimental": [
        "strong_immunogenicity_and_presentation",
        "strong_immunogenicity",
        "presentation",
    ],
    "broad_supportive": [
        "strong_immunogenicity_and_presentation",
        "strong_immunogenicity",
        "presentation",
        "prediction_concordance",
        "binding_evidence",
    ],
}

ENDPOINT_LABELS = {
    "strict_positive_experimental": "Strict positive experimental evidence",
    "broad_supportive": "Broad supportive evidence",
}

# Main scopes for reviewer-facing plots.
# HLA-specific scopes are still summarized if columns are present, but peptide-only is used for the main plots.
PLOT_SCOPES = ["exact_peptide", "relaxed_min8_peptide"]

SCOPE_LABELS = {
    "exact_peptide": "Exact peptide",
    "exact_peptide_hla": "Exact peptide + HLA",
    "relaxed_min8_peptide": "Relaxed min-8 peptide",
    "relaxed_min8_peptide_hla": "Relaxed min-8 peptide + HLA",
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

    # Backward-compatible exact naming: "__peptide_match" without relaxed/min8.
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


def make_endpoint_flags(df):
    rows = []
    scopes = [
        "exact_peptide",
        "exact_peptide_hla",
        "relaxed_min8_peptide",
        "relaxed_min8_peptide_hla",
    ]

    for endpoint, tokens in ENDPOINTS.items():
        for scope in scopes:
            cols = evidence_columns(df, tokens, scope)
            flag_col = f"endpoint__{endpoint}__{scope}"
            df[flag_col] = bool_series(df, cols).astype(int)

            rows.append({
                "endpoint": endpoint,
                "endpoint_label": ENDPOINT_LABELS[endpoint],
                "match_scope": scope,
                "match_scope_label": SCOPE_LABELS[scope],
                "n_columns_used": len(cols),
                "columns_used": ";".join(cols),
                "n_positive_candidates": int(df[flag_col].sum()),
                "positive_fraction": float(df[flag_col].mean()) if len(df) else np.nan,
            })

    endpoint_inventory = pd.DataFrame(rows)
    endpoint_inventory.to_csv(OUTDIR / "endpoint_flag_inventory.tsv", sep="\t", index=False)
    return df, endpoint_inventory


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

    # Optional: include Tumor DNA VAF if it exists.
    vaf_col = first_existing(df, ["Tumor_DNA_VAF", "Tumor DNA VAF", "DNA_VAF", "tumor_vaf"])
    if vaf_col is not None:
        method_specs.append(("tumor_dna_vaf", "Tumor DNA VAF", [vaf_col], "higher_better"))

    for method_id, label, candidates, direction in method_specs:
        col = first_existing(df, candidates)
        if col is not None:
            methods.append((method_id, label, col, direction))

    pd.DataFrame(
        [
            {"ranking_method": m, "ranking_label": lab, "ranking_column": col, "direction": direction}
            for m, lab, col, direction in methods
        ]
    ).to_csv(OUTDIR / "ranking_methods_used.tsv", sep="\t", index=False)

    return methods


def top_indices_for_group(g, col, direction, top_type, cutoff):
    values = pd.to_numeric(g[col], errors="coerce")
    tmp = g.copy()
    tmp["_rank_value"] = values

    # Missing values go to bottom.
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
        c for c in df.columns
        if c.startswith("endpoint__")
    ]

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
                            "endpoint_label": ENDPOINT_LABELS.get(endpoint, endpoint),
                            "match_scope": scope,
                            "match_scope_label": SCOPE_LABELS.get(scope, scope),
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
    per_patient.to_csv(OUTDIR / "per_patient_ranking_delta_vs_ic50.tsv", sep="\t", index=False)
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

    summary.to_csv(OUTDIR / "summary_ranking_delta_vs_ic50.tsv", sep="\t", index=False)
    return summary


def savefig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_series(summary):
    for endpoint in ENDPOINTS:
        for scope in PLOT_SCOPES:
            for top_type in ["percent", "topn"]:
                sub = summary[
                    (summary["endpoint"] == endpoint)
                    & (summary["match_scope"] == scope)
                    & (summary["top_type"] == top_type)
                ].copy()

                if sub.empty:
                    continue

                plt.figure(figsize=(10, 6))

                for method, mdf in sub.groupby("ranking_method", sort=False):
                    mdf = mdf.sort_values("cutoff")
                    label = mdf["ranking_label"].iloc[0]
                    plt.plot(
                        mdf["cutoff"],
                        mdf["median_delta_hit_rate_vs_ic50"],
                        marker="o",
                        label=label,
                    )

                plt.axhline(0, linestyle="--", linewidth=1)
                if top_type == "percent":
                    xlabel = "Top-percent cutoff within each patient"
                    title_suffix = "top-percent"
                else:
                    xlabel = "Top-N candidates within each patient"
                    title_suffix = "top-N"

                plt.xlabel(xlabel)
                plt.ylabel("Median delta hit rate vs IC50")
                plt.title(f"All ranking methods vs IC50: {ENDPOINT_LABELS[endpoint]}\n{scope}, {title_suffix}")
                plt.legend()
                safe = f"{endpoint}__{scope}__{top_type}".replace("/", "_")
                savefig(FIGDIR / f"{safe}.png")


def plot_per_patient_borda(summary_per_patient):
    # A compact per-patient plot for the most useful reviewer-facing case:
    # Borda vs IC50, top 20%, relaxed min-8 peptide, strict endpoint.
    sub = summary_per_patient[
        (summary_per_patient["endpoint"] == "strict_positive_experimental")
        & (summary_per_patient["match_scope"] == "relaxed_min8_peptide")
        & (summary_per_patient["top_type"] == "percent")
        & (summary_per_patient["cutoff"] == 20)
        & (summary_per_patient["ranking_method"] == "borda")
    ].copy()

    if sub.empty:
        return

    sub = sub.sort_values("delta_hits_vs_ic50", ascending=False)

    plt.figure(figsize=(12, 5))
    plt.bar(sub["sample_id"], sub["delta_hits_vs_ic50"])
    plt.axhline(0, linestyle="--", linewidth=1)
    plt.xticks(rotation=70, ha="right")
    plt.ylabel("Delta top-hit count vs IC50")
    plt.title("Per-patient Borda improvement vs IC50, top 20%, relaxed min-8 strict endpoint")
    savefig(FIGDIR / "per_patient_borda_vs_ic50_top20pct_relaxed_min8_strict.png")


def write_markdown(summary, endpoint_inventory, methods):
    lines = []
    lines.append("# Complete-rescue-with-VAF ranking benchmark vs IC50")
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append("For each patient, candidates were ranked by IC50, Borda, or one individual feature. At each top cutoff, the hit rate for an evidence endpoint was calculated among the top-ranked candidates. The plotted value is:")
    lines.append("")
    lines.append("```text")
    lines.append("delta hit rate = hit rate from ranking method - hit rate from IC50 ranking")
    lines.append("```")
    lines.append("")
    lines.append("Positive values mean the method placed more evidence-supported candidates near the top than IC50 alone. Values near zero mean no improvement over IC50. Negative values mean fewer evidence-supported candidates than IC50.")
    lines.append("")
    lines.append("## Endpoint definitions")
    lines.append("")
    lines.append("- **Strict positive experimental evidence** includes `strong_immunogenicity_and_presentation`, `strong_immunogenicity`, and `presentation`.")
    lines.append("- **Broad supportive evidence** includes the strict endpoint plus `prediction_concordance` and `binding_evidence`.")
    lines.append("")
    lines.append("## Ranking methods used")
    lines.append("")
    for method_id, label, col, direction in methods:
        lines.append(f"- {label}: `{col}`, {direction}")
    lines.append("")
    lines.append("## Evidence flag counts")
    lines.append("")
    keep = endpoint_inventory[
        endpoint_inventory["match_scope"].isin([
            "exact_peptide",
            "exact_peptide_hla",
            "relaxed_min8_peptide",
            "relaxed_min8_peptide_hla",
        ])
    ].copy()

    lines.append("| Endpoint | Match scope | Positive candidates | Columns used |")
    lines.append("|---|---|---:|---:|")
    for _, r in keep.iterrows():
        lines.append(
            f"| {r['endpoint_label']} | {r['match_scope_label']} | {int(r['n_positive_candidates'])} | {int(r['n_columns_used'])} |"
        )

    lines.append("")
    lines.append("## Key output files")
    lines.append("")
    lines.append("- `endpoint_flag_inventory.tsv`: evidence endpoint definitions and counts.")
    lines.append("- `ranking_methods_used.tsv`: ranking columns and directions.")
    lines.append("- `per_patient_ranking_delta_vs_ic50.tsv`: per-patient deltas.")
    lines.append("- `summary_ranking_delta_vs_ic50.tsv`: median/mean delta summaries across patients.")
    lines.append("- `figures/`: top-percent and top-N plots.")

    (OUTDIR / "README_ranking_benchmark_summary.md").write_text("\n".join(lines))


def main():
    print("Reading:", INFILE)
    df = pd.read_csv(INFILE, sep="\t", low_memory=False)

    if "candidate_id" not in df.columns:
        df = df.reset_index(drop=True)
        df["candidate_id"] = np.arange(1, len(df) + 1)

    print("Rows:", len(df))
    print("Columns:", len(df.columns))

    df, endpoint_inventory = make_endpoint_flags(df)
    methods = build_ranking_methods(df)

    print("Ranking methods:")
    for method_id, label, col, direction in methods:
        print(f"  {method_id}: {col} ({direction})")

    print("Computing per-patient ranking deltas...")
    per_patient = compute_per_patient(df, methods)

    print("Summarizing...")
    summary = summarize(per_patient)

    print("Plotting...")
    plot_series(summary)
    plot_per_patient_borda(per_patient)

    write_markdown(summary, endpoint_inventory, methods)

    print("Done.")
    print("Output folder:", OUTDIR)
    print("Main summary:", OUTDIR / "summary_ranking_delta_vs_ic50.tsv")
    print("Figures:", FIGDIR)


if __name__ == "__main__":
    main()
