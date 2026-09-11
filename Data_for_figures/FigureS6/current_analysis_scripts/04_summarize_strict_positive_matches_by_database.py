#!/usr/bin/env python3

from pathlib import Path
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(".").resolve()
MATCHES = ROOT / "match_examples.exact_relaxed.tsv.gz"
FLAGS = ROOT / "egg_candidates.exact_relaxed_evidence_flags.tsv.gz"

OUTDIR = ROOT / "strict_positive_match_database_summary"
FIGDIR = OUTDIR / "figures"
OUTDIR.mkdir(exist_ok=True)
FIGDIR.mkdir(exist_ok=True)

STRICT_STRENGTHS = {
    "strong_immunogenicity_and_presentation",
    "strong_immunogenicity",
    "presentation",
}


def norm(x):
    return re.sub(r"[^a-z0-9]+", "_", str(x).strip().lower()).strip("_")


def first_existing(df, names):
    mapping = {norm(c): c for c in df.columns}
    for name in names:
        key = norm(name)
        if key in mapping:
            return mapping[key]
    return None


def truthy(s):
    if s.dtype == bool:
        return s.fillna(False)
    numeric = pd.to_numeric(s, errors="coerce").fillna(0)
    string = s.astype(str).str.strip().str.lower()
    return (numeric > 0) | string.isin(["true", "yes", "y", "1"])


def classify_scope(raw):
    s = str(raw).lower()

    if "relaxed" in s or "min8" in s or "contain" in s:
        prefix = "relaxed_min8"
    else:
        prefix = "exact"

    if "hla" in s:
        return f"{prefix}_peptide_hla"
    return f"{prefix}_peptide"


def build_long_match_table(m):
    database_col = first_existing(m, ["database", "reference_database", "db"])
    strength_col = first_existing(m, ["evidence_strength", "reference_evidence_strength"])
    class_col = first_existing(m, ["evidence_class", "reference_evidence_class"])
    scope_col = first_existing(m, ["match_scope", "match_type", "scope", "match_kind"])

    candidate_col = first_existing(m, ["candidate_id"])
    sample_col = first_existing(m, ["sample_id", "sample"])
    peptide_col = first_existing(m, ["MT.Epitope.Seq", "candidate_peptide", "peptide", "query_peptide"])
    hla_col = first_existing(m, ["HLA.Allele", "candidate_hla", "hla_normalized", "hla"])
    gene_col = first_existing(m, ["Gene.Name", "gene", "candidate_gene"])

    required = {
        "database_col": database_col,
        "strength_col": strength_col,
    }
    missing = [k for k, v in required.items() if v is None]
    if missing:
        print("Columns available in match table:")
        for c in m.columns:
            print(" ", c)
        raise ValueError(f"Could not identify required columns: {missing}")

    if candidate_col is None:
        # Fallback candidate key.
        key_cols = [c for c in [sample_col, peptide_col, hla_col, gene_col] if c is not None]
        if not key_cols:
            raise ValueError("No candidate_id and no usable fallback candidate key columns.")
        m = m.copy()
        m["candidate_id_fallback"] = m[key_cols].astype(str).agg("||".join, axis=1)
        candidate_col = "candidate_id_fallback"

    keep_cols = [candidate_col, database_col, strength_col]
    optional_cols = [class_col, sample_col, peptide_col, hla_col, gene_col]
    keep_cols += [c for c in optional_cols if c is not None]
    keep_cols = list(dict.fromkeys(keep_cols))

    if scope_col is not None:
        long = m[keep_cols + [scope_col]].copy()
        long["match_scope_clean"] = long[scope_col].map(classify_scope)
    else:
        # Fallback: unpivot boolean match columns.
        match_flag_cols = [
            c for c in m.columns
            if "match" in norm(c)
            and (
                "exact" in norm(c)
                or "relaxed" in norm(c)
                or "min8" in norm(c)
                or "contain" in norm(c)
            )
        ]

        if not match_flag_cols:
            print("Columns available in match table:")
            for c in m.columns:
                print(" ", c)
            raise ValueError("Could not find match_scope column or boolean match flag columns.")

        pieces = []
        for c in match_flag_cols:
            sub = m.loc[truthy(m[c]), keep_cols].copy()
            sub["match_scope_clean"] = classify_scope(c)
            pieces.append(sub)

        long = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame()

    rename = {
        candidate_col: "candidate_id",
        database_col: "database",
        strength_col: "evidence_strength",
    }
    if class_col is not None:
        rename[class_col] = "evidence_class"
    if sample_col is not None:
        rename[sample_col] = "sample_id"
    if peptide_col is not None:
        rename[peptide_col] = "candidate_peptide"
    if hla_col is not None:
        rename[hla_col] = "candidate_hla"
    if gene_col is not None:
        rename[gene_col] = "candidate_gene"

    long = long.rename(columns=rename)

    for c in ["evidence_class", "sample_id", "candidate_peptide", "candidate_hla", "candidate_gene"]:
        if c not in long.columns:
            long[c] = ""

    long["evidence_strength"] = long["evidence_strength"].astype(str)
    long["database"] = long["database"].astype(str)

    return long


def savefig(path):
    plt.tight_layout()
    plt.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()


def main():
    print("Reading matches:", MATCHES)
    m = pd.read_csv(MATCHES, sep="\t", low_memory=False)

    long = build_long_match_table(m)

    long.to_csv(OUTDIR / "all_matches_long_scope_table.tsv.gz", sep="\t", index=False, compression="gzip")

    strict = long[long["evidence_strength"].isin(STRICT_STRENGTHS)].copy()
    strict.to_csv(OUTDIR / "strict_positive_matches_long_scope_table.tsv.gz", sep="\t", index=False, compression="gzip")

    # Candidate-level scope flags per database and evidence strength.
    key_cols = ["database", "evidence_strength", "candidate_id"]
    cand_scope = (
        strict
        .drop_duplicates(key_cols + ["match_scope_clean"])
        .assign(value=1)
        .pivot_table(
            index=key_cols,
            columns="match_scope_clean",
            values="value",
            aggfunc="max",
            fill_value=0,
        )
        .reset_index()
    )

    for c in [
        "exact_peptide",
        "exact_peptide_hla",
        "relaxed_min8_peptide",
        "relaxed_min8_peptide_hla",
    ]:
        if c not in cand_scope.columns:
            cand_scope[c] = 0

    cand_scope["relaxed_min8_peptide_nonexact"] = (
        (cand_scope["relaxed_min8_peptide"] == 1)
        & (cand_scope["exact_peptide"] == 0)
    ).astype(int)

    cand_scope["relaxed_min8_peptide_hla_nonexact"] = (
        (cand_scope["relaxed_min8_peptide_hla"] == 1)
        & (cand_scope["exact_peptide_hla"] == 0)
    ).astype(int)

    cand_scope.to_csv(OUTDIR / "strict_positive_candidate_scope_flags_by_database.tsv", sep="\t", index=False)

    by_db_strength = (
        cand_scope
        .groupby(["database", "evidence_strength"], dropna=False)
        .agg(
            n_unique_candidates_any_strict_match=("candidate_id", "nunique"),
            exact_peptide_candidates=("exact_peptide", "sum"),
            exact_peptide_hla_candidates=("exact_peptide_hla", "sum"),
            relaxed_min8_peptide_candidates=("relaxed_min8_peptide", "sum"),
            relaxed_min8_peptide_hla_candidates=("relaxed_min8_peptide_hla", "sum"),
            relaxed_min8_peptide_nonexact_candidates=("relaxed_min8_peptide_nonexact", "sum"),
            relaxed_min8_peptide_hla_nonexact_candidates=("relaxed_min8_peptide_hla_nonexact", "sum"),
        )
        .reset_index()
        .sort_values(["database", "evidence_strength"])
    )

    by_db_strength.to_csv(OUTDIR / "strict_positive_match_counts_by_database_and_strength.tsv", sep="\t", index=False)

    by_db = (
        cand_scope
        .groupby(["database"], dropna=False)
        .agg(
            n_unique_candidates_any_strict_match=("candidate_id", "nunique"),
            exact_peptide_candidates=("exact_peptide", "sum"),
            exact_peptide_hla_candidates=("exact_peptide_hla", "sum"),
            relaxed_min8_peptide_candidates=("relaxed_min8_peptide", "sum"),
            relaxed_min8_peptide_hla_candidates=("relaxed_min8_peptide_hla", "sum"),
            relaxed_min8_peptide_nonexact_candidates=("relaxed_min8_peptide_nonexact", "sum"),
            relaxed_min8_peptide_hla_nonexact_candidates=("relaxed_min8_peptide_hla_nonexact", "sum"),
        )
        .reset_index()
        .sort_values("n_unique_candidates_any_strict_match", ascending=False)
    )

    by_db.to_csv(OUTDIR / "strict_positive_match_counts_by_database.tsv", sep="\t", index=False)

    # Reference-row level counts too.
    ref_rows = (
        strict
        .groupby(["database", "evidence_strength", "match_scope_clean"], dropna=False)
        .size()
        .reset_index(name="n_match_rows")
        .sort_values(["database", "evidence_strength", "match_scope_clean"])
    )
    ref_rows.to_csv(OUTDIR / "strict_positive_reference_match_rows_by_database_strength_scope.tsv", sep="\t", index=False)

    # Plot database counts.
    plot_df = by_db.copy()
    x = np.arange(len(plot_df))
    width = 0.35

    plt.figure(figsize=(8.5, 5))
    plt.bar(x - width / 2, plot_df["exact_peptide_candidates"], width=width, label="Exact peptide")
    plt.bar(x + width / 2, plot_df["relaxed_min8_peptide_nonexact_candidates"], width=width, label="Relaxed min-8 only")
    plt.xticks(x, plot_df["database"], rotation=25, ha="right")
    plt.ylabel("Unique candidate/database matches")
    plt.title("Strict positive evidence matches by database")
    plt.legend()
    plt.grid(axis="y", alpha=0.25)
    savefig(FIGDIR / "strict_positive_exact_vs_relaxed_nonexact_by_database")

    print()
    print("Wrote:", OUTDIR / "strict_positive_match_counts_by_database.tsv")
    print("Wrote:", OUTDIR / "strict_positive_match_counts_by_database_and_strength.tsv")
    print("Wrote:", FIGDIR / "strict_positive_exact_vs_relaxed_nonexact_by_database.png")

    print()
    print("Quick database-level summary:")
    print(by_db.to_string(index=False))


if __name__ == "__main__":
    main()
