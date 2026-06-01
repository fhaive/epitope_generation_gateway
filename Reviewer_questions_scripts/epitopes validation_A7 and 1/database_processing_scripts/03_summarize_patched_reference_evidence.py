#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict, Counter
import gzip
import csv
import math

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


PATCHED_REF_DIR = Path("normalized_reference_tables_hla_patched").resolve()
OUTDIR = Path("patched_reference_evidence_summary").resolve()
FIGDIR = OUTDIR / "figures"

OUTDIR.mkdir(exist_ok=True)
FIGDIR.mkdir(exist_ok=True)


STRICT_HLA_LEVELS = {"preexisting_normalized", "exact_allele"}


EVIDENCE_STRENGTH_TO_FOCUS = {
    # Highest-value evidence for reviewer response
    "strong_immunogenicity_and_presentation": "1_functional_and_presented",
    "strong_immunogenicity": "2_functional_tcell",
    "presentation": "3_experimental_presentation",

    # Useful but weaker / secondary
    "prediction_concordance": "4_prediction_concordance",
    "binding_evidence": "5_binding_only",

    # Useful as controls / exclusions, not positive evidence
    "negative_tcell_assay": "6_negative_evidence",
    "negative_presentation_assay": "6_negative_evidence",
    "negative_binding_assay": "6_negative_evidence",

    # Not a primary benchmark endpoint
    "uncategorized_mhc_ligand": "7_uncategorized_or_other",
}


FOCUS_DESCRIPTION = {
    "1_functional_and_presented": "Highest confidence: reported as immunogenic and supported by presentation/validation tier.",
    "2_functional_tcell": "Functional immunogenicity evidence, usually positive T-cell assay or validated neoantigen evidence.",
    "3_experimental_presentation": "Experimental MHC ligand / immunopeptidomics presentation evidence.",
    "4_prediction_concordance": "Predicted/shared TSNAdb evidence; useful as secondary support, not experimental validation.",
    "5_binding_only": "Binding assay evidence only; weaker than presentation or T-cell activity.",
    "6_negative_evidence": "Negative assay evidence; useful for controls, not positive benchmarking.",
    "7_uncategorized_or_other": "Other/uncategorized evidence; not used as primary positive endpoint.",
}


FOCUS_USE = {
    "1_functional_and_presented": "primary_positive",
    "2_functional_tcell": "primary_positive",
    "3_experimental_presentation": "primary_positive",
    "4_prediction_concordance": "secondary_supportive",
    "5_binding_only": "contextual_only",
    "6_negative_evidence": "negative_control_or_exclude",
    "7_uncategorized_or_other": "exclude_or_sensitivity",
}


DATABASE_ORDER = ["IEDB", "CEDAR", "TSNAdb_v2"]
FOCUS_ORDER = [
    "1_functional_and_presented",
    "2_functional_tcell",
    "3_experimental_presentation",
    "4_prediction_concordance",
    "5_binding_only",
    "6_negative_evidence",
    "7_uncategorized_or_other",
]


def clean(x):
    if x is None:
        return ""
    x = str(x).strip()
    if x.upper() in {"", "NA", "NAN", "NONE", "NULL"}:
        return ""
    return x


def peptide_key(x):
    return clean(x).upper()


def peptide_hla_key(peptide, hla):
    p = peptide_key(peptide)
    h = clean(hla)
    if not p or not h:
        return ""
    return f"{p}|{h}"


def evidence_focus(evidence_strength):
    e = clean(evidence_strength)
    return EVIDENCE_STRENGTH_TO_FOCUS.get(e, "7_uncategorized_or_other")


def strict_hla_eligible(row):
    hla = clean(row.get("hla_normalized", ""))
    level = clean(row.get("hla_match_level", ""))

    if not hla:
        return False

    if level in STRICT_HLA_LEVELS:
        return True

    # Fallback for older rows if hla_match_level was missing.
    if not level:
        return True

    return False


def savefig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()


def init_stats():
    return {
        "n_rows": 0,
        "rows_with_peptide": 0,
        "rows_with_hla_original": 0,
        "rows_with_hla_normalized": 0,
        "rows_strict_hla_eligible": 0,
        "rows_with_peptide_and_strict_hla": 0,
        "unique_peptides": set(),
        "unique_peptide_hla_pairs": set(),
        "unique_references": set(),
    }


def update_stats(stats, row):
    stats["n_rows"] += 1

    peptide = peptide_key(row.get("peptide", ""))
    hla_original = clean(row.get("hla_original", ""))
    hla_normalized = clean(row.get("hla_normalized", ""))
    pmid = clean(row.get("pmid", ""))
    reference_title = clean(row.get("reference_title", ""))

    is_strict = strict_hla_eligible(row)

    if peptide:
        stats["rows_with_peptide"] += 1
        stats["unique_peptides"].add(peptide)

    if hla_original:
        stats["rows_with_hla_original"] += 1

    if hla_normalized:
        stats["rows_with_hla_normalized"] += 1

    if is_strict:
        stats["rows_strict_hla_eligible"] += 1

    ph = peptide_hla_key(peptide, hla_normalized)
    if peptide and is_strict and ph:
        stats["rows_with_peptide_and_strict_hla"] += 1
        stats["unique_peptide_hla_pairs"].add(ph)

    ref = pmid if pmid else reference_title
    if ref:
        stats["unique_references"].add(ref)


def finalize_stats(stats):
    out = dict(stats)
    out["n_unique_peptides"] = len(stats["unique_peptides"])
    out["n_unique_peptide_hla_pairs_strict"] = len(stats["unique_peptide_hla_pairs"])
    out["n_unique_references"] = len(stats["unique_references"])

    del out["unique_peptides"]
    del out["unique_peptide_hla_pairs"]
    del out["unique_references"]

    return out


def plot_grouped_bar(df, x_col, y_col, group_col, title, ylabel, out_png, order_x=None, order_group=None):
    if df.empty:
        return

    plot_df = df.copy()

    if order_x is None:
        order_x = list(plot_df[x_col].drop_duplicates())
    if order_group is None:
        order_group = list(plot_df[group_col].drop_duplicates())

    matrix = []
    for group in order_group:
        vals = []
        for x in order_x:
            sub = plot_df[(plot_df[x_col] == x) & (plot_df[group_col] == group)]
            vals.append(float(sub[y_col].sum()) if not sub.empty else 0.0)
        matrix.append(vals)

    x_positions = list(range(len(order_x)))
    width = 0.8 / max(len(order_group), 1)

    plt.figure(figsize=(11, 5.5))

    for i, group in enumerate(order_group):
        offsets = [x + (i - (len(order_group) - 1) / 2) * width for x in x_positions]
        plt.bar(offsets, matrix[i], width=width, label=group)

    plt.xticks(x_positions, order_x, rotation=25, ha="right")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend(fontsize=8)
    savefig(out_png)


def plot_stacked_bar(df, x_col, stack_col, y_col, title, ylabel, out_png, order_x=None, order_stack=None):
    if df.empty:
        return

    plot_df = df.copy()

    if order_x is None:
        order_x = list(plot_df[x_col].drop_duplicates())
    if order_stack is None:
        order_stack = list(plot_df[stack_col].drop_duplicates())

    bottom = [0.0] * len(order_x)
    x_positions = list(range(len(order_x)))

    plt.figure(figsize=(10, 5.5))

    for stack in order_stack:
        vals = []
        for x in order_x:
            sub = plot_df[(plot_df[x_col] == x) & (plot_df[stack_col] == stack)]
            vals.append(float(sub[y_col].sum()) if not sub.empty else 0.0)

        plt.bar(x_positions, vals, bottom=bottom, label=stack)
        bottom = [b + v for b, v in zip(bottom, vals)]

    plt.xticks(x_positions, order_x, rotation=25, ha="right")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend(fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    savefig(out_png)


def main():
    files = sorted(PATCHED_REF_DIR.glob("*.normalized.tsv.gz"))

    if not files:
        raise FileNotFoundError(f"No patched normalized files found in: {PATCHED_REF_DIR}")

    database_stats = defaultdict(init_stats)
    focus_stats = defaultdict(init_stats)
    hla_match_stats = defaultdict(Counter)
    evidence_class_stats = defaultdict(init_stats)

    for f in files:
        print(f"Reading {f.name}")

        with gzip.open(f, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")

            required = {
                "database",
                "evidence_group",
                "evidence_class",
                "evidence_strength",
                "peptide",
                "hla_original",
                "hla_normalized",
                "hla_class",
            }

            missing = required - set(reader.fieldnames or [])
            if missing:
                raise RuntimeError(f"{f.name} missing required columns: {missing}")

            for row in reader:
                database = clean(row.get("database", "unknown"))
                evidence_strength = clean(row.get("evidence_strength", ""))
                evidence_class = clean(row.get("evidence_class", ""))
                match_level = clean(row.get("hla_match_level", "")) or "missing_hla_match_level"
                focus = evidence_focus(evidence_strength)

                update_stats(database_stats[database], row)
                update_stats(focus_stats[(database, focus)], row)
                update_stats(evidence_class_stats[(database, evidence_class, evidence_strength, focus)], row)

                hla_match_stats[(database, match_level)]["n_rows"] += 1

                if clean(row.get("hla_normalized", "")):
                    hla_match_stats[(database, match_level)]["rows_with_hla_normalized"] += 1

                if strict_hla_eligible(row):
                    hla_match_stats[(database, match_level)]["rows_strict_hla_eligible"] += 1

    # -------------------------
    # Table 1: database overview
    # -------------------------
    rows = []
    for database, stats in database_stats.items():
        row = {"database": database}
        row.update(finalize_stats(stats))
        rows.append(row)

    overview = pd.DataFrame(rows)

    if not overview.empty:
        overview["strict_hla_eligible_fraction"] = (
            overview["rows_strict_hla_eligible"] / overview["n_rows"]
        )
        overview = overview.sort_values(
            "database",
            key=lambda s: s.map({d: i for i, d in enumerate(DATABASE_ORDER)}).fillna(999),
        )

    overview.to_csv(OUTDIR / "table1_database_overview.tsv", sep="\t", index=False)

    # -------------------------
    # Table 2: evidence focus
    # -------------------------
    rows = []
    for (database, focus), stats in focus_stats.items():
        row = {
            "database": database,
            "evidence_focus": focus,
            "use_for_benchmarking": FOCUS_USE.get(focus, "unknown"),
            "description": FOCUS_DESCRIPTION.get(focus, ""),
        }
        row.update(finalize_stats(stats))
        rows.append(row)

    focus_df = pd.DataFrame(rows)

    if not focus_df.empty:
        focus_df["fraction_strict_hla_eligible"] = (
            focus_df["rows_strict_hla_eligible"] / focus_df["n_rows"]
        )
        focus_df["fraction_peptide_and_strict_hla"] = (
            focus_df["rows_with_peptide_and_strict_hla"] / focus_df["n_rows"]
        )
        focus_df = focus_df.sort_values(
            ["database", "evidence_focus"],
            key=lambda s: s.map({**{d: i for i, d in enumerate(DATABASE_ORDER)},
                                 **{e: i for i, e in enumerate(FOCUS_ORDER)}}).fillna(999)
            if s.name in {"database", "evidence_focus"} else s,
        )

    focus_df.to_csv(OUTDIR / "table2_evidence_focus_by_database.tsv", sep="\t", index=False)

    # -------------------------
    # Table 3: HLA matchability
    # -------------------------
    rows = []
    for (database, match_level), counter in hla_match_stats.items():
        rows.append({
            "database": database,
            "hla_match_level": match_level,
            "n_rows": counter["n_rows"],
            "rows_with_hla_normalized": counter["rows_with_hla_normalized"],
            "rows_strict_hla_eligible": counter["rows_strict_hla_eligible"],
        })

    hla_df = pd.DataFrame(rows)

    if not hla_df.empty:
        hla_df = hla_df.sort_values(["database", "n_rows"], ascending=[True, False])

    hla_df.to_csv(OUTDIR / "table3_hla_matchability_by_database.tsv", sep="\t", index=False)

    # -------------------------
    # Optional detail table: evidence class
    # -------------------------
    rows = []
    for (database, evidence_class, evidence_strength, focus), stats in evidence_class_stats.items():
        row = {
            "database": database,
            "evidence_class": evidence_class,
            "evidence_strength": evidence_strength,
            "evidence_focus": focus,
            "use_for_benchmarking": FOCUS_USE.get(focus, "unknown"),
        }
        row.update(finalize_stats(stats))
        rows.append(row)

    detail_df = pd.DataFrame(rows)

    if not detail_df.empty:
        detail_df = detail_df.sort_values(["database", "n_rows"], ascending=[True, False])

    detail_df.to_csv(OUTDIR / "supplement_evidence_class_detail.tsv", sep="\t", index=False)

    # -------------------------
    # Plot 1: rows by evidence focus/database
    # -------------------------
    plot1_df = focus_df.copy()
    plot1_df["evidence_focus_label"] = plot1_df["evidence_focus"].str.replace(r"^\d+_", "", regex=True)

    plot_stacked_bar(
        plot1_df,
        x_col="database",
        stack_col="evidence_focus_label",
        y_col="n_rows",
        title="Reference database rows by evidence focus",
        ylabel="Rows",
        out_png=FIGDIR / "plot1_rows_by_evidence_focus_and_database.png",
        order_x=[d for d in DATABASE_ORDER if d in set(plot1_df["database"])],
        order_stack=[
            x.replace(r"^\d+_", "")
            for x in []
        ],
    )

    # Redo plot1 with explicit stack order labels.
    stack_order = [
        "functional_and_presented",
        "functional_tcell",
        "experimental_presentation",
        "prediction_concordance",
        "binding_only",
        "negative_evidence",
        "uncategorized_or_other",
    ]

    plot_stacked_bar(
        plot1_df,
        x_col="database",
        stack_col="evidence_focus_label",
        y_col="n_rows",
        title="Reference database rows by evidence focus",
        ylabel="Rows",
        out_png=FIGDIR / "plot1_rows_by_evidence_focus_and_database.png",
        order_x=[d for d in DATABASE_ORDER if d in set(plot1_df["database"])],
        order_stack=stack_order,
    )

    # -------------------------
    # Plot 2: unique peptides by evidence focus/database
    # -------------------------
    plot_grouped_bar(
        plot1_df[
            plot1_df["evidence_focus"].isin([
                "1_functional_and_presented",
                "2_functional_tcell",
                "3_experimental_presentation",
                "4_prediction_concordance",
                "5_binding_only",
            ])
        ],
        x_col="evidence_focus_label",
        y_col="n_unique_peptides",
        group_col="database",
        title="Unique peptide coverage by evidence focus",
        ylabel="Unique peptides",
        out_png=FIGDIR / "plot2_unique_peptides_by_evidence_focus.png",
        order_x=[
            "functional_and_presented",
            "functional_tcell",
            "experimental_presentation",
            "prediction_concordance",
            "binding_only",
        ],
        order_group=[d for d in DATABASE_ORDER if d in set(plot1_df["database"])],
    )

    # -------------------------
    # Plot 3: HLA matchability
    # -------------------------
    hla_plot = hla_df.copy()

    keep_levels = [
        "preexisting_normalized",
        "exact_allele",
        "low_resolution",
        "broad_locus",
        "broad_or_ambiguous",
        "nonhuman_mhc",
        "multiple_or_complex",
        "unresolved",
        "missing",
        "missing_hla_match_level",
    ]

    hla_plot = hla_plot[hla_plot["hla_match_level"].isin(keep_levels)].copy()

    plot_stacked_bar(
        hla_plot,
        x_col="database",
        stack_col="hla_match_level",
        y_col="n_rows",
        title="HLA matchability after conservative normalization",
        ylabel="Rows",
        out_png=FIGDIR / "plot3_hla_matchability_by_database.png",
        order_x=[d for d in DATABASE_ORDER if d in set(hla_plot["database"])],
        order_stack=keep_levels,
    )

    # -------------------------
    # Markdown summary
    # -------------------------
    md = []
    md.append("# Patched reference evidence summary\n")
    md.append("## Evidence categories used for benchmarking\n")
    md.append("| Evidence focus | Use | Interpretation |")
    md.append("|---|---|---|")

    for focus in FOCUS_ORDER:
        md.append(
            f"| {focus} | {FOCUS_USE.get(focus, '')} | {FOCUS_DESCRIPTION.get(focus, '')} |"
        )

    md.append("\n## Recommended positive evidence endpoints\n")
    md.append(
        "- **Primary positive experimental/validated endpoint:** "
        "`1_functional_and_presented`, `2_functional_tcell`, and `3_experimental_presentation`."
    )
    md.append(
        "- **Secondary supportive endpoint:** `4_prediction_concordance`, mainly TSNAdb predicted/shared evidence."
    )
    md.append(
        "- **Not primary positive evidence:** binding-only, negative evidence, broad/uncategorized evidence."
    )
    md.append(
        "- **Strict peptide-HLA matching:** use only rows where `hla_match_level` is "
        "`preexisting_normalized` or `exact_allele`."
    )
    md.append(
        "- **Low-resolution HLA labels** such as HLA-DR1, HLA-DQ8, or HLA-Cw6 are retained separately "
        "but should not be mixed with strict allele-level matches."
    )

    md.append("\n## Output files\n")
    md.append("- `table1_database_overview.tsv`")
    md.append("- `table2_evidence_focus_by_database.tsv`")
    md.append("- `table3_hla_matchability_by_database.tsv`")
    md.append("- `supplement_evidence_class_detail.tsv`")
    md.append("- `figures/plot1_rows_by_evidence_focus_and_database.png`")
    md.append("- `figures/plot2_unique_peptides_by_evidence_focus.png`")
    md.append("- `figures/plot3_hla_matchability_by_database.png`")

    with open(OUTDIR / "README_summary.md", "w") as f:
        f.write("\n".join(md) + "\n")

    print()
    print(f"Wrote summary folder: {OUTDIR}")
    print()
    print("Quick overview:")
    print(overview.to_string(index=False))


if __name__ == "__main__":
    main()
