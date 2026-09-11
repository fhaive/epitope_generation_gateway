#!/usr/bin/env python3

from pathlib import Path
import gzip
import csv
import re
from collections import defaultdict
import pandas as pd

ROOT = Path(".").resolve()
PROJECT_ROOT = Path("/data/fsluma/pipelines/Epitope_Generation_Gateway")

CANDIDATES = ROOT / "egg_candidates.exact_relaxed_evidence_flags.tsv.gz"
REFERENCE_DIR = PROJECT_ROOT / "Immunopeptidomics_Data/TCGA_analysis/normalized_reference_tables_hla_patched"

OUTDIR = ROOT / "verify_strict_positive_exact_vs_relaxed_fast"
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
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", str(x).strip().upper())


def substrings_min8(seq):
    seq = clean_peptide(seq)
    n = len(seq)
    for i in range(n):
        for j in range(i + MIN_OVERLAP, n + 1):
            yield seq[i:j]


def open_gz(path):
    return gzip.open(path, "rt", newline="")


def main():
    print("Reading EGG candidates")
    cand = pd.read_csv(CANDIDATES, sep="\t", low_memory=False)

    peptide_col = None
    for c in ["MT.Epitope.Seq", "candidate_peptide", "peptide", "mt_epitope_seq"]:
        if c in cand.columns:
            peptide_col = c
            break
    if peptide_col is None:
        raise ValueError(f"Could not find peptide column. Columns: {list(cand.columns)}")

    if "candidate_id" not in cand.columns:
        key_cols = [c for c in ["sample_id", "Gene.Name", "HLA.Allele", peptide_col, "Borda_Rank"] if c in cand.columns]
        cand["candidate_id"] = cand[key_cols].astype(str).agg("||".join, axis=1)

    cand["pep_clean"] = cand[peptide_col].map(clean_peptide)
    cand = cand[cand["pep_clean"].str.len() >= MIN_OVERLAP].copy()

    candidate_peptides = set(cand["pep_clean"])
    peptide_to_candidate_ids = (
        cand.groupby("pep_clean")["candidate_id"]
        .apply(lambda x: set(x))
        .to_dict()
    )

    # Index every >=8 aa substring found inside EGG candidate peptides.
    # This makes reference-contained-in-candidate matching fast.
    candidate_substring_to_candidate_ids = defaultdict(set)
    for pep, ids in peptide_to_candidate_ids.items():
        for sub in substrings_min8(pep):
            candidate_substring_to_candidate_ids[sub].update(ids)

    print(f"Candidate rows: {len(cand)}")
    print(f"Unique candidate peptides: {len(candidate_peptides)}")
    print(f"Unique candidate substrings length >=8: {len(candidate_substring_to_candidate_ids)}")

    stats = defaultdict(lambda: {
        "strict_reference_rows": 0,
        "unique_strict_reference_peptides": set(),
        "exact_candidate_ids": set(),
        "relaxed_candidate_ids": set(),
        "exact_reference_peptides_matched": set(),
        "relaxed_reference_peptides_matched": set(),
    })

    examples_exact = []
    examples_relaxed = []

    reference_files = sorted(REFERENCE_DIR.glob("*.normalized.tsv.gz"))
    if not reference_files:
        raise FileNotFoundError(f"No reference files found in {REFERENCE_DIR}")

    for ref in reference_files:
        print(f"Scanning {ref.name}")
        with open_gz(ref) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                strength = row.get("evidence_strength", "")
                if strength not in STRICT_STRENGTHS:
                    continue

                db = row.get("database", "")
                ref_pep = clean_peptide(row.get("peptide", ""))
                if len(ref_pep) < MIN_OVERLAP:
                    continue

                keys = [(db, "ALL_STRICT"), (db, strength)]

                for key in keys:
                    stats[key]["strict_reference_rows"] += 1
                    stats[key]["unique_strict_reference_peptides"].add(ref_pep)

                # Exact peptide-only match.
                exact_ids = peptide_to_candidate_ids.get(ref_pep, set())
                if exact_ids:
                    for key in keys:
                        stats[key]["exact_candidate_ids"].update(exact_ids)
                        stats[key]["exact_reference_peptides_matched"].add(ref_pep)

                    if len(examples_exact) < 50:
                        examples_exact.append({
                            "database": db,
                            "evidence_strength": strength,
                            "candidate_peptide": ref_pep,
                            "reference_peptide": ref_pep,
                            "match_type": "exact",
                            "reference_hla_original": row.get("hla_original", ""),
                            "reference_hla_normalized": row.get("hla_normalized", ""),
                            "pmid": row.get("pmid", ""),
                            "reference_title": row.get("reference_title", ""),
                        })

                relaxed_ids = set()

                # Case 1: database peptide is contained in an EGG candidate peptide.
                relaxed_ids.update(candidate_substring_to_candidate_ids.get(ref_pep, set()))

                # Case 2: EGG candidate peptide is contained in the database peptide.
                # Check candidate peptides by generating substrings of the reference peptide.
                for sub in substrings_min8(ref_pep):
                    if sub in peptide_to_candidate_ids:
                        relaxed_ids.update(peptide_to_candidate_ids[sub])

                if relaxed_ids:
                    for key in keys:
                        stats[key]["relaxed_candidate_ids"].update(relaxed_ids)
                        stats[key]["relaxed_reference_peptides_matched"].add(ref_pep)

                    if len(examples_relaxed) < 100:
                        examples_relaxed.append({
                            "database": db,
                            "evidence_strength": strength,
                            "reference_peptide": ref_pep,
                            "reference_length": len(ref_pep),
                            "n_candidate_ids_matched": len(relaxed_ids),
                            "match_type": "relaxed_min8",
                            "reference_hla_original": row.get("hla_original", ""),
                            "reference_hla_normalized": row.get("hla_normalized", ""),
                            "pmid": row.get("pmid", ""),
                            "reference_title": row.get("reference_title", ""),
                        })

    rows = []
    for (db, strength), d in sorted(stats.items()):
        rows.append({
            "database": db,
            "evidence_strength": strength,
            "strict_reference_rows": d["strict_reference_rows"],
            "unique_strict_reference_peptides": len(d["unique_strict_reference_peptides"]),
            "egg_candidates_exact_peptide_match": len(d["exact_candidate_ids"]),
            "reference_peptides_exactly_matched": len(d["exact_reference_peptides_matched"]),
            "egg_candidates_relaxed_min8_match": len(d["relaxed_candidate_ids"]),
            "reference_peptides_relaxed_min8_matched": len(d["relaxed_reference_peptides_matched"]),
        })

    out = pd.DataFrame(rows)
    by_db = out[out["evidence_strength"] == "ALL_STRICT"].copy()
    by_strength = out[out["evidence_strength"] != "ALL_STRICT"].copy()

    by_db.to_csv(OUTDIR / "strict_positive_exact_vs_relaxed_by_database.tsv", sep="\t", index=False)
    by_strength.to_csv(OUTDIR / "strict_positive_exact_vs_relaxed_by_database_and_strength.tsv", sep="\t", index=False)
    pd.DataFrame(examples_exact).to_csv(OUTDIR / "exact_match_examples.tsv", sep="\t", index=False)
    pd.DataFrame(examples_relaxed).to_csv(OUTDIR / "relaxed_match_examples.tsv", sep="\t", index=False)

    print()
    print("Database-level summary:")
    print(by_db.to_string(index=False))
    print()
    print("Wrote:", OUTDIR)


if __name__ == "__main__":
    main()
