#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


BASE = Path("TCGA_melanoma")
R3 = BASE / "rescue_final_analysis" / "reviewer_updates" / "reviewer_3_borda_sensitivity_v2"

SUMMARY = R3 / "tables" / "reviewer3_rank_stability_summary.tsv"
PLOT_DIR = R3 / "plots"

PLOT_DIR.mkdir(parents=True, exist_ok=True)


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


def clean_metric_name(s):
    return (
        s.replace("leave_out_", "")
         .replace("Net_", "")
         .replace("_", " ")
    )


def get_summary():
    df = pd.read_csv(SUMMARY, sep="\t")
    df = df[df["scenario"] != "recomputed_default"].copy()
    return df


def plot_main_plus_leave_one_spearman(df):
    """
    This plot shows Spearman correlation with the default ranking over the full ranked candidate list.
    It now includes the leave-one-network-metric-out scenarios.
    """
    set_plot_style()

    main_order = [
        "mhc_ic50_only",
        "depmap_only",
        "network_only",
        "no_network_mhc_ic50_depmap",
        "no_mhc_ic50_depmap_network",
        "no_depmap_mhc_ic50_network",
    ]

    leave_order = [
        "leave_out_Net_Betweenness",
        "leave_out_Net_Degree",
        "leave_out_Net_Impact",
        "leave_out_Net_Strength",
        "leave_out_Net_WCI",
    ]

    order = main_order + leave_order

    labels = {
        "mhc_ic50_only": "MHC-IC50\nonly",
        "depmap_only": "DepMap\nonly",
        "network_only": "Network\nonly",
        "no_network_mhc_ic50_depmap": "No network\nMHC+DepMap",
        "no_mhc_ic50_depmap_network": "No MHC-IC50\nDepMap+network",
        "no_depmap_mhc_ic50_network": "No DepMap\nMHC+network",
        "leave_out_Net_Betweenness": "Remove\nBetweenness",
        "leave_out_Net_Degree": "Remove\nDegree",
        "leave_out_Net_Impact": "Remove\nImpact",
        "leave_out_Net_Strength": "Remove\nStrength",
        "leave_out_Net_WCI": "Remove\nWCI",
    }

    plot = df.set_index("scenario").loc[order].reset_index()
    x = np.arange(len(plot))

    colors = ["#1f77b4"] * len(main_order) + ["#ff7f0e"] * len(leave_order)

    fig = plt.figure(figsize=(12.2, 4.9))
    ax = plt.gca()

    bars = ax.bar(
        x,
        plot["median_spearman_rho_vs_default"],
        color=colors,
        edgecolor="black",
        linewidth=0.6,
    )

    ax.errorbar(
        x=x,
        y=plot["median_spearman_rho_vs_default"],
        yerr=[
            plot["median_spearman_rho_vs_default"] - plot["q1_spearman_rho_vs_default"],
            plot["q3_spearman_rho_vs_default"] - plot["median_spearman_rho_vs_default"],
        ],
        fmt="none",
        ecolor="black",
        elinewidth=1.0,
        capsize=3,
    )

    for i, row in plot.iterrows():
        y = row["median_spearman_rho_vs_default"]
        ax.text(
            i,
            min(y + 0.035, 1.02),
            f"{y:.2f}",
            ha="center",
            va="bottom",
            fontsize=8,
            rotation=0,
        )

    ax.axvline(len(main_order) - 0.5, color="black", linestyle=":", linewidth=1.0)

    ax.text(
        (len(main_order) - 1) / 2,
        1.06,
        "Feature/block sensitivity",
        ha="center",
        va="bottom",
        fontsize=10,
        transform=ax.get_xaxis_transform(),
    )

    ax.text(
        len(main_order) + (len(leave_order) - 1) / 2,
        1.06,
        "Leave-one-network-metric-out",
        ha="center",
        va="bottom",
        fontsize=10,
        transform=ax.get_xaxis_transform(),
    )

    ax.set_xticks(x)
    ax.set_xticklabels([labels[s] for s in order])
    ax.set_ylabel("Median Spearman correlation\nwith default rank")
    ax.set_ylim(0, 1.10)
    ax.set_title("Full-rank similarity under alternative Borda weighting scenarios")

    plt.setp(ax.get_xticklabels(), rotation=35, ha="right")

    fig.tight_layout()

    png = PLOT_DIR / "reviewer3_main_scenarios_spearman_vs_default.png"
    svg = PLOT_DIR / "reviewer3_main_scenarios_spearman_vs_default.svg"

    fig.savefig(png, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)

    print("Wrote:", png)
    print("Wrote:", svg)


def plot_leave_one_top10_top20_overlap(df):
    set_plot_style()

    leave_order = [
        "leave_out_Net_Betweenness",
        "leave_out_Net_Degree",
        "leave_out_Net_Impact",
        "leave_out_Net_Strength",
        "leave_out_Net_WCI",
    ]

    leave = df.set_index("scenario").loc[leave_order].reset_index()
    leave["metric_removed"] = leave["scenario"].map(clean_metric_name)

    x = np.arange(len(leave))
    width = 0.38

    fig = plt.figure(figsize=(7.8, 4.7))
    ax = plt.gca()

    ax.bar(
        x - width / 2,
        leave["median_top10_overlap_fraction"],
        width,
        label="Top-10 overlap",
        edgecolor="black",
        linewidth=0.6,
    )

    ax.bar(
        x + width / 2,
        leave["median_top20_overlap_fraction"],
        width,
        label="Top-20 overlap",
        edgecolor="black",
        linewidth=0.6,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(leave["metric_removed"])
    ax.set_ylabel("Median overlap with default ranking")
    ax.set_ylim(0, 1.05)
    ax.set_title("Leave-one-network-metric-out overlap")

    ax.legend(
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        title="Overlap set",
    )

    plt.setp(ax.get_xticklabels(), rotation=35, ha="right")

    fig.tight_layout(rect=[0, 0, 0.82, 1])

    png = PLOT_DIR / "reviewer3_leave_one_network_metric_top10_top20_overlap.png"
    svg = PLOT_DIR / "reviewer3_leave_one_network_metric_top10_top20_overlap.svg"

    fig.savefig(png, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)

    print("Wrote:", png)
    print("Wrote:", svg)


def plot_leave_one_spearman(df):
    set_plot_style()

    leave_order = [
        "leave_out_Net_Betweenness",
        "leave_out_Net_Degree",
        "leave_out_Net_Impact",
        "leave_out_Net_Strength",
        "leave_out_Net_WCI",
    ]

    leave = df.set_index("scenario").loc[leave_order].reset_index()
    leave["metric_removed"] = leave["scenario"].map(clean_metric_name)

    x = np.arange(len(leave))

    fig = plt.figure(figsize=(7.4, 4.7))
    ax = plt.gca()

    ax.bar(
        x,
        leave["median_spearman_rho_vs_default"],
        edgecolor="black",
        linewidth=0.6,
    )

    ax.errorbar(
        x=x,
        y=leave["median_spearman_rho_vs_default"],
        yerr=[
            leave["median_spearman_rho_vs_default"] - leave["q1_spearman_rho_vs_default"],
            leave["q3_spearman_rho_vs_default"] - leave["median_spearman_rho_vs_default"],
        ],
        fmt="none",
        ecolor="black",
        elinewidth=1.0,
        capsize=3,
    )

    for i, row in leave.iterrows():
        y = row["median_spearman_rho_vs_default"]
        ax.text(
            i,
            min(y + 0.025, 1.02),
            f"{y:.3f}",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(leave["metric_removed"])
    ax.set_ylabel("Median Spearman correlation\nwith default rank")
    ax.set_ylim(0, 1.08)
    ax.set_title("Leave-one-network-metric-out full-rank similarity")

    plt.setp(ax.get_xticklabels(), rotation=35, ha="right")

    fig.tight_layout()

    png = PLOT_DIR / "reviewer3_leave_one_network_metric_spearman_vs_default.png"
    svg = PLOT_DIR / "reviewer3_leave_one_network_metric_spearman_vs_default.svg"

    fig.savefig(png, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)

    print("Wrote:", png)
    print("Wrote:", svg)


def main():
    df = get_summary()

    # Export a focused leave-one-network summary table.
    leave = df[df["scenario"].str.startswith("leave_out_")].copy()
    leave["metric_removed"] = leave["scenario"].map(clean_metric_name)
    leave.to_csv(R3 / "tables" / "reviewer3_leave_one_network_metric_summary.tsv", sep="\t", index=False)

    plot_main_plus_leave_one_spearman(df)
    plot_leave_one_top10_top20_overlap(df)
    plot_leave_one_spearman(df)


if __name__ == "__main__":
    main()
