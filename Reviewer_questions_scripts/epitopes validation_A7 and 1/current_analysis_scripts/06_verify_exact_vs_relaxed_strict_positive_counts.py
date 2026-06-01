#!/usr/bin/env python3

from pathlib import Path
import gzip
import csv
import re
from collections import defaultdict

import pandas as pd


ROOT = Path(".").resolve()
PROJECT_ROOT = Path("/data/fsluma/pipelines/Epitope_Generation_Gateway")

CANDIDATE_FLAGS = ROOT / "egg_candidates.exact_relaxed_evidence_flags.tsv.gz"
REFERENCE_DIR = PROJECT_ROOT / "Immunopeptidomics_Data/TCGA_analysis/normalized_reference_tables_hla_patched"

OUTDIR = ROOT / "verify_strict_positive_exact_vs_relaxed"
OUTDIR.mkdir(exist_ok=True)

STRICT_STRENGTHS = {
    "strong_immunogenicity_and_presentation",
    "strong_immunogenicity",
    "presentation",
}

MIN_OVERLAP = 8


def clean_peptide(x):
    if pd.isna(x):
        return ""
    s = str(x).strip().upper()
    s = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", s)
    return s


def containment_match(a, b, min_len=8):
    """
    Relaxed min-8 containment:
    True if one full peptide is contained inside the other,
    and the contained peptide is at least min_len aa.
    """
    a = clean_peptide(a)
    b = clean_peptide(b)
    if not a or not b:
        return False, "", 0, ""

    if a == b and len(a) >= min_len:
        return True, a, len(a), "exact"

    if len(a) >= min_len and a in b:
        return True, a, len(a), "candidate_in_reference"

    if len(b) >= min_len and b in a:
        return True, b, len(b), "reference_in_candidate"

    return False, "", 0, ""


def open_maybe_gz(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "rt", newline="")


def main():
    print("Reading EGG candidate table:")
    print(CANDIDATE_FLAGS)

    cand = pd.read_csv(CANDIDATE_FLAGS, sep="\t", low_memory=False)

    # Find candidate columns robustly.
    peptide_col = None
    for c in ["MT.Epitope.Seq", "mt_epitope_seq", "candidate_peptide", "peptide"]:
        if c in cand.columns:
            peptide_col = c
            break
    if peptide_col is None:
        raise ValueError(f"Could not find peptide column. Columns: {list(cand.columns)}")

    if "candidate_id" not in cand.columns:
        key_cols = []
        for c in ["sample_id", "Gene.Name", "HLA.Allele", peptide_col, "Borda_Rank"]:
            if c in cand.columns:
                key_cols.append(c)
        cand["candidate_id"] = cand[key_cols].astype(str).agg("||".join, axis=1)

    cand["candidate_peptide_clean"] = cand[peptide_col].map(clean_peptide)
    cand = cand[cand["candidate_peptide_clean"].str.len() > 0].copy()

    candidate_peptides = sorted(cand["candidate_peptide_clean"].unique())
    peptide_to_candidate_ids = (
        cand.groupby("candidate_peptide_clean")["candidate_id"]
        .apply(lambda x: sorted(set(x)))
        .to_dict()
    )

    print(f"Candidate rows: {len(cand)}")
    print(f"Unique candidate peptides: {len(candidate_peptides)}")

    ref_files = sorted(REFERENCE_DIR.glob("*.normalized.tsv.gz"))
    if not ref_files:
        raise FileNotFoundError(f"No reference files found in {REFERENCE_DIR}")

    database_summary = []
    strength_summary = []
    exact_examples = []
    relaxed_examples = []

    # Store per database and strength.
    stats = defaultdict(lambda: {
        "strict_reference_rows": 0,
        "unique_strict_reference_peptides": set(),
        "exact_candidate_ids": set(),
        "relaxed_candidate_ids": set(),
        "exact_candidate_peptides": set(),
        "relaxed_candidate_peptides": set(),
        "matched_reference_peptides_exact": set(),
        "matched_reference_peptides_relaxed": set(),
    })

    # Pre-index candidate peptides for exact lookup.
    candidate_peptide_set = set(candidate_peptides)

    for ref_path in ref_files:
        print(f"Scanning {ref_path.name}")

        with open_maybe_gz(ref_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")

            required_cols = {"database", "evidence_strength", "peptide"}
            missing = required_cols - set(reader.fieldnames or [])
            if missing:
                raise ValueError(f"{ref_path} missing columns: {missing}")

            for row in reader:
                db = row.get("database", "")
                strength = row.get("evidence_strength", "")
                if strength not in STRICT_STRENGTHS:
                    continue

                ref_pep = clean_peptide(row.get("peptide", ""))
                if len(ref_pep) < MIN_OVERLAP:
                    continue

                keys = [
                    (db, "ALL_STRICT"),
                    (db, strength),
                ]

                for key in keys:
                    stats[key]["strict_reference_rows"] += 1
                    stats[key]["unique_strict_reference_peptides"].add(ref_pep)

                # Exact peptide-only match, ignoring HLA.
                if ref_pep in candidate_peptide_set:
                    for cid in peptide_to_candidate_ids[ref_pep]:
                        for key in keys:
                            stats[key]["exact_candidate_ids"].add(cid)
                            stats[key]["exact_candidate_peptides"].add(ref_pep)
                            stats[key]["matched_reference_peptides_exact"].add(ref_pep)

                    if len(exact_examples) < 200:
                        exact_examples.append({
                            "database": db,
                            "evidence_strength": strength,
                            "candidate_peptide": ref_pep,
                            "reference_peptide": ref_pep,
                            "match_type": "exact",
                            "reference_hla_original": row.get("hla_original", ""),
                            "reference_hla_normalized": row.get("hla_normalized", ""),
                            "reference_title": row.get("reference_title", ""),
                            "pmid": row.get("pmid", ""),
                        })

                # Relaxed containment match, ignoring HLA.
                # This brute-force loop is okay for ~2k candidate unique peptides.
                for cand_pep in candidate_peptides:
                    ok, shared, shared_len, match_type = containment_match(cand_pep, ref_pep, MIN_OVERLAP)
                    if not ok:
                        continue

                    for cid in peptide_to_candidate_ids[cand_pep]:
                        for key in keys:
                            stats[key]["relaxed_candidate_ids"].add(cid)
                            stats[key]["relaxed_candidate_peptides"].add(cand_pep)
                            stats[key]["matched_reference_peptides_relaxed"].add(ref_pep)

                    if len(relaxed_examples) < 500:
                        relaxed_examples.append({
                            "database": db,
                            "evidence_strength": strength,
                            "candidate_peptide": cand_pep,
                            "reference_peptide": ref_pep,
                            "candidate_length": len(cand_pep),
                            "reference_length": len(ref_pep),
                            "shared_sequence": shared,
                            "shared_length": shared_len,
                            "match_type": match_type,
                            "reference_hla_original": row.get("hla_original", ""),
                            "reference_hla_normalized": row.get("hla_normalized", ""),
                            "reference_title": row.get("reference_title", ""),
                            "pmid": row.get("pmid", ""),
                        })

    rows = []
    for (db, strength), d in sorted(stats.items()):
        rows.append({
            "database": db,
            "evidence_strength": strength,
            "strict_reference_rows": d["strict_reference_rows"],
            "unique_strict_reference_peptides": len(d["unique_strict_reference_peptides"]),
            "egg_candidates_exact_peptide_match": len(d["exact_candidate_ids"]),
            "egg_candidate_peptides_exact_match": len(d["exact_candidate_peptides"]),
            "reference_peptides_exactly_matched": len(d["matched_reference_peptides_exact"]),
            "egg_candidates_relaxed_min8_match": len(d["relaxed_candidate_ids"]),
            "egg_candidate_peptides_relaxed_min8_match": len(d["relaxed_candidate_peptides"]),
            "reference_peptides_relaxed_min8_matched": len(d["matched_reference_peptides_relaxed"]),
        })

    summary = pd.DataFrame(rows)

    by_database = summary[summary["evidence_strength"] == "ALL_STRICT"].copy()
    by_strength = summary[summary["evidence_strength"] != "ALL_STRICT"].copy()

    by_database.to_csv(
        OUTDIR / "strict_positive_exact_vs_relaxed_by_database.tsv",
        sep="\t",
        index=False,
    )

    by_strength.to_csv(
        OUTDIR / "strict_positive_exact_vs_relaxed_by_database_and_strength.tsv",
        sep="\t",
        index=False,
    )

    pd.DataFrame(exact_examples).to_csv(
        OUTDIR / "strict_positive_exact_match_examples.tsv",
        sep="\t",
        index=False,
    )

    pd.DataFrame(relaxed_examples).to_csv(
        OUTDIR / "strict_positive_relaxed_min8_match_examples.tsv",
        sep="\t",
        index=False,
    )

    print()
    print("Wrote:")
    print(OUTDIR / "strict_positive_exact_vs_relaxed_by_database.tsv")
    print(OUTDIR / "strict_positive_exact_vs_relaxed_by_database_and_strength.tsv")
    print(OUTDIR / "strict_positive_exact_match_examples.tsv")
    print(OUTDIR / "strict_positive_relaxed_min8_match_examples.tsv")

    print()
    print("Database-level verification:")
    print(by_database.to_string(index=False))


if __name__ == "__main__":
    main()
