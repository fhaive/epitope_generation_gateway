#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


BASE = Path("TCGA_melanoma")
FINAL_DIR = BASE / "epitopes_prioritisation_complete_rescue" / "final_epitopes"

OUT_DIR = BASE / "rescue_final_analysis" / "reviewer_updates" / "reviewer_3_borda_sensitivity_v2"
PLOT_DIR = OUT_DIR / "plots"
TABLE_DIR = OUT_DIR / "tables"

PLOT_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)


NETWORK_COLS = [
    "Net_Betweenness",
    "Net_Degree",
    "Net_Impact",
    "Net_Strength",
    "Net_WCI",
]

IC50_CANDIDATES = [
    "Median.MT.IC50.Score",
    "Median MT IC50 Score",
    "Median_MT_IC50_Score",
    "Median IC50 Score",
]

DEPMAP_CANDIDATES = [
    "Depmap_survivability_score",
    "DepMap_survivability_score",
    "Depmap.Survivability.Score",
    "DepMap_Survivability_Score",
]

DEFAULT_RANK_CANDIDATES = [
    "Borda_Rank",
    "Borda Rank",
    "Borda.Rank",
    "Final_Rank",
    "Final Rank",
    "Rank",
]

DEFAULT_SCORE_CANDIDATES = [
    "Borda_Score",
    "Borda Score",
    "Borda.Score",
]

DRIVER_CANDIDATES = [
    "Intogen_Driver_role",
    "Intogen Driver role",
    "Intogen.Driver.role",
    "Cancer_Driver",
    "Driver",
]


def find_col(df, candidates, required=False):
    for c in candidates:
        if c in df.columns:
            return c
    if required:
        raise ValueError(
            f"Could not find any of {candidates}. Available columns:\n"
            + "\n".join(map(str, df.columns))
        )
    return None


def sample_from_path(path):
    return path.name.replace("_epitopes_final.csv", "")


def read_final_tables():
    files = sorted(FINAL_DIR.glob("*_epitopes_final.csv"))
    if not files:
        raise FileNotFoundError(f"No final tables found in {FINAL_DIR}")

    rows = []
    for f in files:
        df = pd.read_csv(f, dtype=str, low_memory=False).fillna("")
        df["sample"] = sample_from_path(f)
        df["source_file"] = str(f)
        rows.append(df)

    df = pd.concat(rows, ignore_index=True)

    ic50_col = find_col(df, IC50_CANDIDATES, required=True)
    depmap_col = find_col(df, DEPMAP_CANDIDATES, required=True)
    default_rank_col = find_col(df, DEFAULT_RANK_CANDIDATES, required=True)
    default_score_col = find_col(df, DEFAULT_SCORE_CANDIDATES, required=False)
    driver_col = find_col(df, DRIVER_CANDIDATES, required=False)

    missing_network = [c for c in NETWORK_COLS if c not in df.columns]
    if missing_network:
        raise ValueError(f"Missing expected network columns: {missing_network}")

    df["default_rank"] = pd.to_numeric(df[default_rank_col], errors="coerce")
    df["default_score"] = pd.to_numeric(df[default_score_col], errors="coerce") if default_score_col else np.nan

    df["ic50_raw"] = pd.to_numeric(df[ic50_col], errors="coerce")
    df["depmap_raw"] = pd.to_numeric(df[depmap_col], errors="coerce")

    for c in NETWORK_COLS:
        df[c + "_raw"] = pd.to_numeric(df[c], errors="coerce")

    if driver_col:
        x = df[driver_col].astype(str).str.strip()
        df["is_driver"] = ~x.isin(["", "NA", "NaN", "nan", "None", ".", "FALSE", "False", "false"])
    else:
        df["is_driver"] = False

    df = df[df["default_rank"].notna()].copy()

    meta = pd.DataFrame([{
        "ic50_col": ic50_col,
        "ic50_direction": "lower is better",
        "depmap_col": depmap_col,
        "depmap_direction": "lower/more negative is better",
        "default_rank_col": default_rank_col,
        "default_score_col": default_score_col if default_score_col else "",
        "driver_col": driver_col if driver_col else "",
        "network_cols": ";".join(NETWORK_COLS),
        "note": "HLA expression was excluded because it was not used in the Borda ranking."
    }])

    return df, meta


def within_sample_score(df, value_col, direction):
    """
    Returns score where larger is better.
    direction='low': smaller raw value is better.
    direction='high': larger raw value is better.
    """
    out = pd.Series(index=df.index, dtype=float)

    for sample, idx in df.groupby("sample").groups.items():
        x = pd.to_numeric(df.loc[idx, value_col], errors="coerce")
        valid = x.notna()

        if valid.sum() == 0:
            out.loc[idx] = np.nan
            continue

        if valid.sum() == 1:
            tmp = pd.Series(np.nan, index=idx, dtype=float)
            tmp.loc[x[valid].index] = 1.0
            out.loc[idx] = tmp
            continue

        if direction == "high":
            ranks = x.rank(method="average", ascending=True, na_option="keep")
        elif direction == "low":
            ranks = x.rank(method="average", ascending=False, na_option="keep")
        else:
            raise ValueError(direction)

        out.loc[idx] = (ranks - 1) / (valid.sum() - 1)

    return out


def compute_metric_scores(df):
    scored = df.copy()

    scored["score_mhc_ic50"] = within_sample_score(scored, "ic50_raw", direction="low")
    scored["score_depmap"] = within_sample_score(scored, "depmap_raw", direction="low")

    network_score_cols = []
    for c in NETWORK_COLS:
        sc = "score_" + c
        scored[sc] = within_sample_score(scored, c + "_raw", direction="high")
        network_score_cols.append(sc)

    score_cols = ["score_mhc_ic50", "score_depmap"] + network_score_cols

    for c in score_cols:
        scored[c] = scored.groupby("sample")[c].transform(lambda x: x.fillna(x.median()))
        scored[c] = scored[c].fillna(0.5)

    return scored, network_score_cols


def weighted_score(df, weights):
    total = sum(weights.values())
    if total <= 0:
        raise ValueError("Total weight must be positive")

    s = np.zeros(len(df), dtype=float)

    for col, w in weights.items():
        s += df[col].astype(float).values * float(w)

    return s / total


def make_scenarios(network_score_cols):
    scenarios = {}

    # Core sensitivity scenarios.

    scenarios["mhc_ic50_only"] = {
        "score_mhc_ic50": 1.0,
    }

    scenarios["depmap_only"] = {
        "score_depmap": 1.0,
    }

    scenarios["network_only"] = {
        **{c: 1.0 for c in network_score_cols},
    }

    scenarios["no_network_mhc_ic50_depmap"] = {
        "score_mhc_ic50": 0.5,
        "score_depmap": 0.5,
    }

    scenarios["no_mhc_ic50_depmap_network"] = {
        "score_depmap": 0.25,
        **{c: 0.10 for c in network_score_cols},
    }

    scenarios["no_depmap_mhc_ic50_network"] = {
        "score_mhc_ic50": 0.25,
        **{c: 0.10 for c in network_score_cols},
    }

    # Network-weight sweep.
    for total_network_weight in [0.0, 0.25, 0.50, 0.75, 1.0]:
        name = f"network_weight_{int(total_network_weight * 100)}pct"

        if total_network_weight < 1.0:
            non_network_weight = 1.0 - total_network_weight
            mhc_w = non_network_weight / 2.0
            depmap_w = non_network_weight / 2.0
        else:
            mhc_w = 0.0
            depmap_w = 0.0

        net_w_each = total_network_weight / len(network_score_cols)

        weights = {}
        if mhc_w > 0:
            weights["score_mhc_ic50"] = mhc_w
        if depmap_w > 0:
            weights["score_depmap"] = depmap_w
        for c in network_score_cols:
            if net_w_each > 0:
                weights[c] = net_w_each

        scenarios[name] = weights

    # Leave-one-network-metric-out.
    # Keep total network block weight at 0.50 and redistribute across the remaining network metrics.
    for omitted in network_score_cols:
        remaining = [c for c in network_score_cols if c != omitted]
        scenarios[f"leave_out_{omitted.replace('score_', '')}"] = {
            "score_mhc_ic50": 0.25,
            "score_depmap": 0.25,
            **{c: 0.50 / len(remaining) for c in remaining},
        }

    return scenarios


def compute_scenario_ranks(df, scenarios):
    out = df.copy()

    for scenario, weights in scenarios.items():
        score_col = f"{scenario}_score"
        rank_col = f"{scenario}_rank"

        out[score_col] = weighted_score(out, weights)

        out[rank_col] = (
            out.groupby("sample")[score_col]
            .rank(method="first", ascending=False)
        )

    return out


def compare_to_default(df, scenarios):
    rows = []

    for sample, g in df.groupby("sample", dropna=False):
        g = g.copy()

        default_rank = g["default_rank"]
        default_top10 = set(g.nsmallest(10, "default_rank").index)
        default_top20 = set(g.nsmallest(20, "default_rank").index)

        for scenario in scenarios:
            rank_col = f"{scenario}_rank"
            scenario_rank = g[rank_col]

            if len(g) >= 3 and default_rank.nunique() > 1 and scenario_rank.nunique() > 1:
                rho, p = spearmanr(default_rank, scenario_rank)
            else:
                rho, p = np.nan, np.nan

            scen_top10 = set(g.nsmallest(10, rank_col).index)
            scen_top20 = set(g.nsmallest(20, rank_col).index)

            rows.append({
                "sample": sample,
                "scenario": scenario,
                "n_candidates": int(len(g)),
                "spearman_rho_vs_default_rank": rho,
                "spearman_p_vs_default_rank": p,
                "top10_overlap_with_default_n": len(default_top10 & scen_top10),
                "top10_overlap_with_default_fraction": len(default_top10 & scen_top10) / max(1, len(default_top10)),
                "top20_overlap_with_default_n": len(default_top20 & scen_top20),
                "top20_overlap_with_default_fraction": len(default_top20 & scen_top20) / max(1, len(default_top20)),
                "scenario_top10_driver_rows": int(g.loc[list(scen_top10), "is_driver"].sum()) if len(scen_top10) else 0,
                "default_top10_driver_rows": int(g.loc[list(default_top10), "is_driver"].sum()) if len(default_top10) else 0,
            })

    return pd.DataFrame(rows)


def summarize_comparisons(comp):
    rows = []

    for scenario, g in comp.groupby("scenario", dropna=False):
        rows.append({
            "scenario": scenario,
            "n_samples": int(g["sample"].nunique()),
            "median_spearman_rho_vs_default": g["spearman_rho_vs_default_rank"].median(),
            "q1_spearman_rho_vs_default": g["spearman_rho_vs_default_rank"].quantile(0.25),
            "q3_spearman_rho_vs_default": g["spearman_rho_vs_default_rank"].quantile(0.75),
            "median_top10_overlap_fraction": g["top10_overlap_with_default_fraction"].median(),
            "q1_top10_overlap_fraction": g["top10_overlap_with_default_fraction"].quantile(0.25),
            "q3_top10_overlap_fraction": g["top10_overlap_with_default_fraction"].quantile(0.75),
            "median_top20_overlap_fraction": g["top20_overlap_with_default_fraction"].median(),
            "q1_top20_overlap_fraction": g["top20_overlap_with_default_fraction"].quantile(0.25),
            "q3_top20_overlap_fraction": g["top20_overlap_with_default_fraction"].quantile(0.75),
            "total_scenario_top10_driver_rows": int(g["scenario_top10_driver_rows"].sum()),
            "total_default_top10_driver_rows": int(g["default_top10_driver_rows"].sum()),
        })

    return pd.DataFrame(rows).sort_values("scenario")


def pairwise_feature_correlations(scored, score_cols):
    rows = []

    for sample, g in scored.groupby("sample", dropna=False):
        sub = g[score_cols].apply(pd.to_numeric, errors="coerce")
        corr = sub.corr(method="spearman")

        for i, a in enumerate(score_cols):
            for b in score_cols[i + 1:]:
                rows.append({
                    "sample": sample,
                    "metric_a": a.replace("score_", ""),
                    "metric_b": b.replace("score_", ""),
                    "spearman_rho": corr.loc[a, b],
                })

    pair = pd.DataFrame(rows)

    summary = (
        pair.groupby(["metric_a", "metric_b"], as_index=False)
        .agg(
            median_spearman_rho=("spearman_rho", "median"),
            q1_spearman_rho=("spearman_rho", lambda x: x.quantile(0.25)),
            q3_spearman_rho=("spearman_rho", lambda x: x.quantile(0.75)),
            median_abs_spearman_rho=("spearman_rho", lambda x: x.abs().median()),
        )
        .sort_values("median_abs_spearman_rho", ascending=False)
    )

    overall = pd.DataFrame([{
        "all_feature_pairwise_median_abs_spearman_rho": pair["spearman_rho"].abs().median(),
        "all_feature_pairwise_q1_abs_spearman_rho": pair["spearman_rho"].abs().quantile(0.25),
        "all_feature_pairwise_q3_abs_spearman_rho": pair["spearman_rho"].abs().quantile(0.75),
    }])

    return pair, summary, overall


def set_plot_style():
    mpl.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 300,
        "font.size": 11,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 9,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "grid.linestyle": "--",
        "grid.linewidth": 0.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def boxplot_metric(comp, scenarios, labels, value_col, ylabel, outbase):
    set_plot_style()

    data = [
        comp.loc[comp["scenario"] == s, value_col].dropna().values
        for s in scenarios
    ]

    fig = plt.figure(figsize=(11.5, 4.9))
    ax = plt.gca()

    bp = ax.boxplot(data, labels=labels, patch_artist=True, showfliers=True)

    for patch in bp["boxes"]:
        patch.set_facecolor("#1f77b4")
        patch.set_alpha(0.75)
        patch.set_edgecolor("black")
        patch.set_linewidth(0.9)

    for k in ["medians", "whiskers", "caps"]:
        for artist in bp[k]:
            artist.set_color("black")
            artist.set_linewidth(1.0)

    for flier in bp["fliers"]:
        flier.set_marker("o")
        flier.set_markersize(3)
        flier.set_alpha(0.5)
        flier.set_markeredgecolor("black")

    ax.set_ylabel(ylabel)
    ax.set_ylim(-0.02, 1.02)
    plt.setp(ax.get_xticklabels(), rotation=35, ha="right")

    fig.tight_layout()

    png = outbase.with_suffix(".png")
    svg = outbase.with_suffix(".svg")
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)

    print("Wrote:", png)
    print("Wrote:", svg)


def plot_network_weight_sweep(summary):
    set_plot_style()

    rows = []
    for pct in [0, 25, 50, 75, 100]:
        scenario = f"network_weight_{pct}pct"
        row = summary[summary["scenario"] == scenario].iloc[0].to_dict()
        row["network_weight_pct"] = pct
        rows.append(row)

    d = pd.DataFrame(rows)

    fig = plt.figure(figsize=(5.7, 4.4))
    ax = plt.gca()

    ax.plot(
        d["network_weight_pct"],
        d["median_top10_overlap_fraction"],
        marker="o",
        linewidth=1.6,
        label="Top-10 overlap"
    )

    ax.plot(
        d["network_weight_pct"],
        d["median_spearman_rho_vs_default"],
        marker="o",
        linewidth=1.6,
        label="Spearman rho"
    )

    ax.set_xlabel("Total network weight (%)")
    ax.set_ylabel("Similarity to existing default ranking")
    ax.set_ylim(0, 1.02)
    ax.set_title("Effect of changing total network weight")
    ax.legend(frameon=False)

    fig.tight_layout()

    png = PLOT_DIR / "reviewer3_network_weight_sweep_similarity.png"
    svg = PLOT_DIR / "reviewer3_network_weight_sweep_similarity.svg"

    fig.savefig(png, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)

    print("Wrote:", png)
    print("Wrote:", svg)


def plot_leave_one_network(summary):
    set_plot_style()

    leave = summary[summary["scenario"].str.startswith("leave_out_")].copy()
    leave["metric_removed"] = leave["scenario"].str.replace("leave_out_", "", regex=False)

    fig = plt.figure(figsize=(7.2, 4.5))
    ax = plt.gca()

    ax.bar(
        leave["metric_removed"],
        leave["median_top10_overlap_fraction"],
        edgecolor="black",
        linewidth=0.6,
    )

    ax.set_ylabel("Median top-10 overlap with default")
    ax.set_ylim(0, 1.02)
    ax.set_title("Leave-one-network-metric-out sensitivity")
    plt.setp(ax.get_xticklabels(), rotation=35, ha="right")

    fig.tight_layout()

    png = PLOT_DIR / "reviewer3_leave_one_network_metric_top10_overlap.png"
    svg = PLOT_DIR / "reviewer3_leave_one_network_metric_top10_overlap.svg"

    fig.savefig(png, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)

    print("Wrote:", png)
    print("Wrote:", svg)


def main():
    df, meta = read_final_tables()
    scored, network_score_cols = compute_metric_scores(df)

    scenarios = make_scenarios(network_score_cols)
    scenario_df = compute_scenario_ranks(scored, scenarios)

    comp = compare_to_default(scenario_df, scenarios)
    summary = summarize_comparisons(comp)

    feature_score_cols = ["score_mhc_ic50", "score_depmap"] + network_score_cols
    pair_corr, pair_corr_summary, pair_corr_overall = pairwise_feature_correlations(scored, feature_score_cols)

    scenario_weight_rows = []
    for scenario, weights in scenarios.items():
        for metric, weight in weights.items():
            scenario_weight_rows.append({
                "scenario": scenario,
                "metric_score_column": metric,
                "weight": weight,
            })

    pd.DataFrame(scenario_weight_rows).to_csv(TABLE_DIR / "reviewer3_weighting_scenarios.tsv", sep="\t", index=False)
    meta.to_csv(TABLE_DIR / "reviewer3_columns_used.tsv", sep="\t", index=False)
    scenario_df.to_csv(TABLE_DIR / "reviewer3_epitope_scores_and_sensitivity_ranks.tsv.gz", sep="\t", index=False)
    comp.to_csv(TABLE_DIR / "reviewer3_rank_stability_by_sample.tsv", sep="\t", index=False)
    summary.to_csv(TABLE_DIR / "reviewer3_rank_stability_summary.tsv", sep="\t", index=False)
    pair_corr.to_csv(TABLE_DIR / "reviewer3_feature_pairwise_correlations_by_sample.tsv", sep="\t", index=False)
    pair_corr_summary.to_csv(TABLE_DIR / "reviewer3_feature_pairwise_correlations_summary.tsv", sep="\t", index=False)
    pair_corr_overall.to_csv(TABLE_DIR / "reviewer3_feature_pairwise_correlations_overall.tsv", sep="\t", index=False)

    main_scenarios = [
        "mhc_ic50_only",
        "depmap_only",
        "network_only",
        "no_network_mhc_ic50_depmap",
        "no_mhc_ic50_depmap_network",
        "no_depmap_mhc_ic50_network",
    ]

    main_labels = [
        "MHC-IC50\nonly",
        "DepMap\nonly",
        "Network\nonly",
        "No network\nMHC+DepMap",
        "No MHC-IC50\nDepMap+network",
        "No DepMap\nMHC+network",
    ]

    boxplot_metric(
        comp,
        main_scenarios,
        main_labels,
        value_col="top10_overlap_with_default_fraction",
        ylabel="Top-10 overlap with existing default ranking",
        outbase=PLOT_DIR / "reviewer3_main_scenarios_top10_overlap",
    )

    boxplot_metric(
        comp,
        main_scenarios,
        main_labels,
        value_col="spearman_rho_vs_default_rank",
        ylabel="Spearman correlation with existing default rank",
        outbase=PLOT_DIR / "reviewer3_main_scenarios_spearman_vs_default",
    )

    plot_network_weight_sweep(summary)
    plot_leave_one_network(summary)

    print()
    print("Columns used:")
    print(meta.to_string(index=False))

    print()
    print("Rank stability summary:")
    print(summary.to_string(index=False))

    print()
    print("Feature correlation overall:")
    print(pair_corr_overall.to_string(index=False))

    print()
    print("Feature pairwise correlation summary:")
    print(pair_corr_summary.to_string(index=False))


if __name__ == "__main__":
    main()
