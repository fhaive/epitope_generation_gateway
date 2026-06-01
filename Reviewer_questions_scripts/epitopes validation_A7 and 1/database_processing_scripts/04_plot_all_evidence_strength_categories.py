#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict, Counter
import gzip
import csv

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


PATCHED_REF_DIR = Path("normalized_reference_tables_hla_patched").resolve()
OUTDIR = Path("patched_reference_evidence_summary/all_evidence_strength_categories").resolve()
FIGDIR = OUTDIR / "figures"

OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)


STRICT_HLA_LEVELS = {"preexisting_normalized", "exact_allele"}

EVIDENCE_ORDER = [
    "strong_immunogenicity_and_presentation",
    "strong_immunogenicity",
    "presentation",
    "prediction_concordance",
    "binding_evidence",
    "negative_tcell_assay",
    "negative_presentation_assay",
    "negative_binding_assay",
    "uncategorized_mhc_ligand",
]

EVIDENCE_MEANING = {
    "strong_immunogenicity_and_presentation": "Highest-confidence validated neoantigen / immunogenic and presentation-supported evidence.",
    "strong_immunogenicity": "Functional immunogenicity evidence, mainly positive T-cell assays or TSNAdb tier-2 validated neoantigens.",
    "presentation": "Experimental MHC ligand / immunopeptidomics presentation evidence.",
    "prediction_concordance": "Externally predicted/shared TSNAdb neoantigen evidence, not direct experimental validation.",
    "binding_evidence": "Positive MHC binding evidence only.",
    "negative_tcell_assay": "Negative T-cell assay evidence.",
    "negative_presentation_assay": "Negative MHC presentation evidence.",
    "negative_binding_assay": "Negative MHC binding evidence.",
    "uncategorized_mhc_ligand": "MHC ligand records not cleanly classified as binding, presentation, or negative evidence.",
}

USE_FOR_BENCHMARKING = {
    "strong_immunogenicity_and_presentation": "primary_positive",
    "strong_immunogenicity": "primary_positive",
    "presentation": "primary_positive",
    "prediction_concordance": "secondary_supportive",
    "binding_evidence": "contextual_only",
    "negative_tcell_assay": "negative_control_or_exclude",
    "negative_presentation_assay": "negative_control_or_exclude",
    "negative_binding_assay": "negative_control_or_exclude",
    "uncategorized_mhc_ligand": "exclude_or_sensitivity",
}


def clean(x):
    if x is None:
        return ""
    x = str(x).strip()
    if x.upper() in {"", "NA", "NAN", "NONE", "NULL"}:
        return ""
    return x


def peptide_key(x):
    return clean(x).upper()


def strict_hla_eligible(row):
    hla = clean(row.get("hla_normalized", ""))
    level = clean(row.get("hla_match_level", ""))

    if not hla:
        return False

    if level in STRICT_HLA_LEVELS:
        return True

    # Backward compatibility if older normalized files lack hla_match_level.
    if not level:
        return True

    return False


def savefig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    files = sorted(PATCHED_REF_DIR.glob("*.normalized.tsv.gz"))
    if not files:
        raise FileNotFoundError(f"No .normalized.tsv.gz files found in {PATCHED_REF_DIR}")

    stats = defaultdict(lambda: {
        "n_rows": 0,
        "rows_with_peptide": 0,
        "rows_with_hla_original": 0,
        "rows_with_hla_normalized": 0,
        "rows_strict_hla_eligible": 0,
        "rows_with_peptide_and_strict_hla": 0,
        "peptides": set(),
        "peptide_hla_pairs_strict": set(),
        "references": set(),
        "evidence_classes": Counter(),
    })

    for f in files:
        print(f"Reading {f.name}")

        with gzip.open(f, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")

            for row in reader:
                database = clean(row.get("database", "unknown"))
                evidence_strength = clean(row.get("evidence_strength", "uncategorized_mhc_ligand"))
                evidence_class = clean(row.get("evidence_class", ""))
                peptide = peptide_key(row.get("peptide", ""))
                hla_original = clean(row.get("hla_original", ""))
                hla_normalized = clean(row.get("hla_normalized", ""))
                pmid = clean(row.get("pmid", ""))
                reference_title = clean(row.get("reference_title", ""))

                key = (database, evidence_strength)
                s = stats[key]

                s["n_rows"] += 1
                s["evidence_classes"][evidence_class] += 1

                if peptide:
                    s["rows_with_peptide"] += 1
                    s["peptides"].add(peptide)

                if hla_original:
                    s["rows_with_hla_original"] += 1

                if hla_normalized:
                    s["rows_with_hla_normalized"] += 1

                is_strict = strict_hla_eligible(row)

                if is_strict:
                    s["rows_strict_hla_eligible"] += 1

                if peptide and hla_normalized and is_strict:
                    s["rows_with_peptide_and_strict_hla"] += 1
                    s["peptide_hla_pairs_strict"].add(f"{peptide}|{hla_normalized}")

                ref = pmid if pmid else reference_title
                if ref:
                    s["references"].add(ref)

    rows = []

    for (database, evidence_strength), s in stats.items():
        rows.append({
            "database": database,
            "evidence_strength": evidence_strength,
            "use_for_benchmarking": USE_FOR_BENCHMARKING.get(evidence_strength, "unknown"),
            "meaning": EVIDENCE_MEANING.get(evidence_strength, ""),
            "n_rows": s["n_rows"],
            "rows_with_peptide": s["rows_with_peptide"],
            "rows_with_hla_original": s["rows_with_hla_original"],
            "rows_with_hla_normalized": s["rows_with_hla_normalized"],
            "rows_strict_hla_eligible": s["rows_strict_hla_eligible"],
            "rows_with_peptide_and_strict_hla": s["rows_with_peptide_and_strict_hla"],
            "n_unique_peptides": len(s["peptides"]),
            "n_unique_peptide_hla_pairs_strict": len(s["peptide_hla_pairs_strict"]),
            "n_unique_references": len(s["references"]),
            "evidence_classes_merged": ";".join(
                f"{k}:{v}" for k, v in s["evidence_classes"].most_common() if k
            ),
        })

    df = pd.DataFrame(rows)

    order_map = {x: i for i, x in enumerate(EVIDENCE_ORDER)}
    db_order = {"IEDB": 0, "CEDAR": 1, "TSNAdb_v2": 2}

    df["evidence_order"] = df["evidence_strength"].map(order_map).fillna(999).astype(int)
    df["database_order"] = df["database"].map(db_order).fillna(999).astype(int)
    df = df.sort_values(["evidence_order", "database_order", "database"])

    out_table = OUTDIR / "all_evidence_strength_by_database.tsv"
    df.drop(columns=["evidence_order", "database_order"]).to_csv(out_table, sep="\t", index=False)

    # Collapsed table across databases.
    collapsed_rows = []
    for evidence_strength, sub in df.groupby("evidence_strength", sort=False):
        collapsed_rows.append({
            "evidence_strength": evidence_strength,
            "use_for_benchmarking": USE_FOR_BENCHMARKING.get(evidence_strength, "unknown"),
            "meaning": EVIDENCE_MEANING.get(evidence_strength, ""),
            "n_rows": int(sub["n_rows"].sum()),
            "rows_with_peptide": int(sub["rows_with_peptide"].sum()),
            "rows_with_hla_normalized": int(sub["rows_with_hla_normalized"].sum()),
            "rows_strict_hla_eligible": int(sub["rows_strict_hla_eligible"].sum()),
            "rows_with_peptide_and_strict_hla": int(sub["rows_with_peptide_and_strict_hla"].sum()),
            "n_unique_peptides_sum_by_database": int(sub["n_unique_peptides"].sum()),
            "n_unique_peptide_hla_pairs_strict_sum_by_database": int(sub["n_unique_peptide_hla_pairs_strict"].sum()),
            "n_unique_references_sum_by_database": int(sub["n_unique_references"].sum()),
        })

    collapsed = pd.DataFrame(collapsed_rows)
    collapsed["evidence_order"] = collapsed["evidence_strength"].map(order_map).fillna(999).astype(int)
    collapsed = collapsed.sort_values("evidence_order").drop(columns=["evidence_order"])

    collapsed.to_csv(OUTDIR / "all_evidence_strength_collapsed.tsv", sep="\t", index=False)

    # Plot 1: rows by evidence strength and database, log scale.
    plot_df = df.copy()
    pivot_rows = plot_df.pivot_table(
        index="evidence_strength",
        columns="database",
        values="n_rows",
        aggfunc="sum",
        fill_value=0,
    )

    pivot_rows = pivot_rows.reindex([x for x in EVIDENCE_ORDER if x in pivot_rows.index])

    ax = pivot_rows.plot(kind="bar", figsize=(12, 5.5))
    ax.set_yscale("log")
    ax.set_ylabel("Rows, log scale")
    ax.set_xlabel("")
    ax.set_title("All evidence-strength categories by database")
    plt.xticks(rotation=35, ha="right")
    savefig(FIGDIR / "all_evidence_strength_rows_by_database_logscale.png")

    # Plot 2: unique peptides by evidence strength and database, log scale.
    pivot_pep = plot_df.pivot_table(
        index="evidence_strength",
        columns="database",
        values="n_unique_peptides",
        aggfunc="sum",
        fill_value=0,
    )

    pivot_pep = pivot_pep.reindex([x for x in EVIDENCE_ORDER if x in pivot_pep.index])

    ax = pivot_pep.plot(kind="bar", figsize=(12, 5.5))
    ax.set_yscale("log")
    ax.set_ylabel("Unique peptides, log scale")
    ax.set_xlabel("")
    ax.set_title("Unique peptide coverage for all evidence-strength categories")
    plt.xticks(rotation=35, ha="right")
    savefig(FIGDIR / "all_evidence_strength_unique_peptides_by_database_logscale.png")

    # Plot 3: strict peptide-HLA usable rows by evidence strength and database.
    pivot_strict = plot_df.pivot_table(
        index="evidence_strength",
        columns="database",
        values="rows_with_peptide_and_strict_hla",
        aggfunc="sum",
        fill_value=0,
    )

    pivot_strict = pivot_strict.reindex([x for x in EVIDENCE_ORDER if x in pivot_strict.index])

    ax = pivot_strict.plot(kind="bar", figsize=(12, 5.5))
    ax.set_yscale("log")
    ax.set_ylabel("Rows with peptide + strict HLA, log scale")
    ax.set_xlabel("")
    ax.set_title("Strict peptide-HLA matchable evidence by evidence-strength category")
    plt.xticks(rotation=35, ha="right")
    savefig(FIGDIR / "all_evidence_strength_strict_peptide_hla_rows_logscale.png")

    print()
    print(f"Wrote: {OUTDIR}")
    print()
    print("Main table:")
    print(out_table)
    print()
    print("Figures:")
    for p in sorted(FIGDIR.glob("*.png")):
        print(p)


if __name__ == "__main__":
    main()
