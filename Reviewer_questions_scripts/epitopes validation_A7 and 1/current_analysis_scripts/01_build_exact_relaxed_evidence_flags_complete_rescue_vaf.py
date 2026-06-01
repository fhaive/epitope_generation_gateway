#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import re
import numpy as np
import pandas as pd


PROJECT_ROOT = Path("/data/fsluma/pipelines/Epitope_Generation_Gateway")

CANDIDATE_BASE = PROJECT_ROOT / "TCGA_melanoma/epitopes_prioritisation_complete_rescue_with_vaf/final_epitopes"
REFERENCE_DIR = PROJECT_ROOT / "Immunopeptidomics_Data/TCGA_analysis/normalized_reference_tables_hla_patched"

OUT = Path(".").resolve()
OUT.mkdir(exist_ok=True)

MIN_RELAXED_OVERLAP = 8
CHUNKSIZE = 200_000

CANDIDATE_OUT = OUT / "egg_candidates.exact_relaxed_evidence_flags.tsv.gz"
MATCH_EXAMPLES_OUT = OUT / "match_examples.exact_relaxed.tsv.gz"
INVENTORY_OUT = OUT / "candidate_evidence_inventory.tsv"
REFERENCE_SCAN_OUT = OUT / "reference_scan_summary.tsv"

REQUIRED_CANDIDATE_COLUMNS = [
    "Gene.Name",
    "HLA.Allele",
    "MT.Epitope.Seq",
    "Median.MT.IC50.Score",
    "Borda_Rank",
]

REFERENCE_USECOLS = [
    "database",
    "evidence_group",
    "evidence_class",
    "evidence_strength",
    "peptide",
    "peptide_length",
    "hla_original",
    "hla_normalized",
    "hla_class",
    "gene",
    "mutation",
    "cancer_type",
    "mutation_type",
    "pmid",
    "reference_title",
    "journal",
    "year",
    "source_file",
]


def clean_peptide(x):
    if pd.isna(x):
        return ""
    s = str(x).strip().upper()
    if s in {"", "NA", "NAN", "NONE"}:
        return ""
    s = re.sub(r"[^A-Z]", "", s)
    return s


def norm_gene(x):
    if pd.isna(x):
        return ""
    s = str(x).strip().upper()
    if s in {"", "NA", "NAN", "NONE"}:
        return ""
    return re.sub(r"[^A-Z0-9\-]", "", s)


def normalize_hla_one(x):
    if pd.isna(x):
        return ""
    s = str(x).strip().upper()
    if s in {"", "NA", "NAN", "NONE"}:
        return ""

    s = s.replace(" ", "")
    s = s.replace("_", "-")
    s = re.sub(r"^HLA[-:]?", "", s)

    # A02:01, A*02:01, E*01:01, G*01:04, etc.
    m = re.match(r"^(A|B|C|E|F|G)\*?(\d{2})(?::?(\d{2,3}[A-Z]?))?$", s)
    if m:
        gene, g1, g2 = m.groups()
        if g2:
            return f"HLA-{gene}*{g1}:{g2}"
        return f"HLA-{gene}*{g1}"

    # DRB1*01:01, DQA1*05:01, DPB1*04:01, etc.
    m = re.match(r"^(DRB[1345]?|DQA1|DQB1|DPA1|DPB1)\*?(\d{2})(?::?(\d{2,3}[A-Z]?))?$", s)
    if m:
        gene, g1, g2 = m.groups()
        if g2:
            return f"HLA-{gene}*{g1}:{g2}"
        return f"HLA-{gene}*{g1}"

    # Already low-resolution patched values such as HLA-DRB1*01 may pass here after HLA removal.
    m = re.match(r"^(DRB[1345]?|DQA1|DQB1|DPA1|DPB1)\*(\d{2})$", s)
    if m:
        gene, g1 = m.groups()
        return f"HLA-{gene}*{g1}"

    return ""


def normalize_candidate_hla_tokens(x):
    if pd.isna(x):
        return []

    s = str(x).strip()
    if not s:
        return []

    # Split composite class-II alleles, for example:
    # DQA1*02:01-DQB1*03:01
    parts = re.split(r"[;,/|]+", s)
    expanded = []
    for part in parts:
        expanded.extend(part.split("-"))

    tokens = []
    for p in expanded:
        n = normalize_hla_one(p)
        if n:
            tokens.append(n)

    return sorted(set(tokens))


def hla_compatible(candidate_tokens, ref_hla):
    ref = normalize_hla_one(ref_hla)
    if not ref or not candidate_tokens:
        return False

    for ct in candidate_tokens:
        if ct == ref:
            return True

        # Allow low-resolution external annotation to match high-resolution patient allele.
        # Example: HLA-DRB1*01 matches HLA-DRB1*01:01.
        if ct.startswith(ref + ":") or ref.startswith(ct + ":"):
            return True

    return False


def kmers(seq, k=8):
    if len(seq) < k:
        return []
    return [seq[i:i+k] for i in range(0, len(seq) - k + 1)]


def is_relaxed_containment(candidate_pep, ref_pep, min_overlap=8):
    if not candidate_pep or not ref_pep:
        return False
    if min(len(candidate_pep), len(ref_pep)) < min_overlap:
        return False
    return candidate_pep in ref_pep or ref_pep in candidate_pep


def find_candidate_dir():
    if (CANDIDATE_BASE / "HLA").exists():
        hla_dir = CANDIDATE_BASE / "HLA"
        if list(hla_dir.glob("Sample_*_epitopes_final*.csv")):
            return hla_dir
    return CANDIDATE_BASE


def sample_id_from_file(path):
    m = re.search(r"Sample_(TCGA-[A-Z0-9\-]+)_", path.name)
    if m:
        return m.group(1)
    return path.stem.replace("Sample_", "").split("_epitopes")[0]


def load_candidates():
    candidate_dir = find_candidate_dir()
    csvs = sorted(candidate_dir.glob("Sample_*_epitopes_final*.csv"))

    if not csvs:
        raise FileNotFoundError(f"No Sample_*_epitopes_final*.csv files found in {candidate_dir}")

    frames = []

    for f in csvs:
        sample_id = sample_id_from_file(f)
        df = pd.read_csv(f, dtype=str)
        df.insert(0, "source_candidate_file", str(f))
        df.insert(1, "sample_id", sample_id)
        df.insert(2, "candidate_row_in_file", np.arange(1, len(df) + 1))
        frames.append(df)

    cand = pd.concat(frames, ignore_index=True)

    missing = [c for c in REQUIRED_CANDIDATE_COLUMNS if c not in cand.columns]
    if missing:
        raise RuntimeError(f"Missing required candidate columns: {missing}")

    cand.insert(0, "candidate_id", [
        f"{sid}__row{row}"
        for sid, row in zip(cand["sample_id"], cand["candidate_row_in_file"])
    ])

    cand["candidate_peptide_clean"] = cand["MT.Epitope.Seq"].map(clean_peptide)
    cand["candidate_gene_clean"] = cand["Gene.Name"].map(norm_gene)
    cand["candidate_hla_tokens"] = cand["HLA.Allele"].map(normalize_candidate_hla_tokens)
    cand["candidate_hla_normalized_joined"] = cand["candidate_hla_tokens"].map(lambda xs: ";".join(xs))

    before = len(cand)
    cand = cand[cand["candidate_peptide_clean"].str.len() > 0].copy()
    after = len(cand)

    if after < before:
        print(f"Removed {before - after} candidate rows with unusable peptide sequence.")

    return cand


def main():
    print("Loading updated complete-rescue-with-VAF candidate tables...")
    cand = load_candidates()
    cand = cand.reset_index(drop=True)

    print(f"Candidate rows: {len(cand)}")
    print(f"Samples: {cand['sample_id'].nunique()}")
    print(f"Unique peptides: {cand['candidate_peptide_clean'].nunique()}")

    n = len(cand)

    peptide_to_ids = defaultdict(list)
    kmer_to_ids = defaultdict(set)

    id_to_peptide = {}
    id_to_hla_tokens = {}
    id_to_gene = {}

    for i, row in cand.iterrows():
        pep = row["candidate_peptide_clean"]
        peptide_to_ids[pep].append(i)
        id_to_peptide[i] = pep
        id_to_hla_tokens[i] = row["candidate_hla_tokens"]
        id_to_gene[i] = row["candidate_gene_clean"]

        for kmer in set(kmers(pep, MIN_RELAXED_OVERLAP)):
            kmer_to_ids[kmer].add(i)

    flags = {}

    def ensure_flag(col):
        if col not in flags:
            flags[col] = np.zeros(n, dtype=np.uint8)
        return flags[col]

    examples = []
    seen_examples = set()
    scan_rows = []

    def mark_candidate(
        cid,
        scope,
        database,
        evidence_strength,
        evidence_class,
        ref_pep,
        ref_hla_norm,
        ref_gene,
        row_meta,
    ):
        evidence_strength = evidence_strength or "unknown"
        evidence_strength = str(evidence_strength).strip()
        database = database or "unknown_database"

        base_cols = [
            f"{scope}__{evidence_strength}",
            f"{scope}__{database}__{evidence_strength}",
        ]

        hla_match = hla_compatible(id_to_hla_tokens[cid], ref_hla_norm)
        gene_match = bool(ref_gene) and id_to_gene[cid] == norm_gene(ref_gene)

        if hla_match:
            base_cols.extend([
                f"{scope}_hla__{evidence_strength}",
                f"{scope}_hla__{database}__{evidence_strength}",
            ])

        if gene_match:
            base_cols.extend([
                f"{scope}_gene__{evidence_strength}",
                f"{scope}_gene__{database}__{evidence_strength}",
            ])

        for col in base_cols:
            ensure_flag(col)[cid] = 1

        ex_key = (cid, scope, database, evidence_strength, evidence_class)
        if ex_key not in seen_examples:
            seen_examples.add(ex_key)
            cpep = id_to_peptide[cid]
            examples.append({
                "candidate_id": cand.loc[cid, "candidate_id"],
                "sample_id": cand.loc[cid, "sample_id"],
                "candidate_gene": cand.loc[cid, "Gene.Name"],
                "candidate_hla_original": cand.loc[cid, "HLA.Allele"],
                "candidate_hla_normalized_joined": cand.loc[cid, "candidate_hla_normalized_joined"],
                "candidate_peptide": cpep,
                "candidate_peptide_length": len(cpep),
                "match_scope": scope,
                "database": database,
                "evidence_group": row_meta.get("evidence_group", ""),
                "evidence_class": evidence_class,
                "evidence_strength": evidence_strength,
                "reference_peptide": ref_pep,
                "reference_peptide_length": len(ref_pep) if ref_pep else 0,
                "length_difference": abs(len(cpep) - len(ref_pep)) if ref_pep else "",
                "candidate_contains_reference": int(bool(ref_pep and ref_pep in cpep)),
                "reference_contains_candidate": int(bool(ref_pep and cpep in ref_pep)),
                "reference_hla_original": row_meta.get("hla_original", ""),
                "reference_hla_normalized": ref_hla_norm,
                "hla_match": int(hla_match),
                "reference_gene": ref_gene,
                "gene_match": int(gene_match),
                "pmid": row_meta.get("pmid", ""),
                "reference_title": row_meta.get("reference_title", ""),
                "journal": row_meta.get("journal", ""),
                "year": row_meta.get("year", ""),
                "source_file": row_meta.get("source_file", ""),
            })

    reference_files = sorted(REFERENCE_DIR.glob("*.normalized.tsv.gz"))
    if not reference_files:
        raise FileNotFoundError(f"No *.normalized.tsv.gz files found in {REFERENCE_DIR}")

    for ref_file in reference_files:
        print(f"Scanning {ref_file.name}...")
        header = pd.read_csv(ref_file, sep="\t", nrows=0).columns.tolist()
        usecols = [c for c in REFERENCE_USECOLS if c in header]

        rows_seen = 0
        exact_row_hits = 0
        exact_hla_row_hits = 0
        relaxed_row_hits = 0
        relaxed_hla_row_hits = 0

        for chunk in pd.read_csv(
            ref_file,
            sep="\t",
            usecols=usecols,
            dtype=str,
            keep_default_na=False,
            chunksize=CHUNKSIZE,
        ):
            rows_seen += len(chunk)

            for c in REFERENCE_USECOLS:
                if c not in chunk.columns:
                    chunk[c] = ""

            cols = [
                "database",
                "evidence_group",
                "evidence_class",
                "evidence_strength",
                "peptide",
                "hla_original",
                "hla_normalized",
                "hla_class",
                "gene",
                "pmid",
                "reference_title",
                "journal",
                "year",
                "source_file",
            ]

            for row in chunk[cols].itertuples(index=False, name=None):
                row_meta = dict(zip(cols, row))

                ref_pep = clean_peptide(row_meta["peptide"])
                if not ref_pep:
                    continue

                database = row_meta["database"] or ref_file.name.replace(".normalized.tsv.gz", "")
                evidence_strength = row_meta["evidence_strength"] or "unknown"
                evidence_class = row_meta["evidence_class"] or "unknown"
                ref_hla_norm = row_meta["hla_normalized"]
                ref_gene = row_meta["gene"]

                exact_ids = peptide_to_ids.get(ref_pep, [])
                exact_any_hla = False

                if exact_ids:
                    exact_row_hits += 1
                    for cid in exact_ids:
                        if hla_compatible(id_to_hla_tokens[cid], ref_hla_norm):
                            exact_any_hla = True
                        mark_candidate(
                            cid=cid,
                            scope="exact_peptide",
                            database=database,
                            evidence_strength=evidence_strength,
                            evidence_class=evidence_class,
                            ref_pep=ref_pep,
                            ref_hla_norm=ref_hla_norm,
                            ref_gene=ref_gene,
                            row_meta=row_meta,
                        )

                if exact_any_hla:
                    exact_hla_row_hits += 1

                relaxed_candidate_ids = set()
                if len(ref_pep) >= MIN_RELAXED_OVERLAP:
                    for kmer in set(kmers(ref_pep, MIN_RELAXED_OVERLAP)):
                        relaxed_candidate_ids.update(kmer_to_ids.get(kmer, []))

                relaxed_ids = [
                    cid for cid in relaxed_candidate_ids
                    if is_relaxed_containment(id_to_peptide[cid], ref_pep, MIN_RELAXED_OVERLAP)
                ]

                relaxed_any_hla = False
                if relaxed_ids:
                    relaxed_row_hits += 1
                    for cid in relaxed_ids:
                        if hla_compatible(id_to_hla_tokens[cid], ref_hla_norm):
                            relaxed_any_hla = True
                        mark_candidate(
                            cid=cid,
                            scope="relaxed_min8_peptide",
                            database=database,
                            evidence_strength=evidence_strength,
                            evidence_class=evidence_class,
                            ref_pep=ref_pep,
                            ref_hla_norm=ref_hla_norm,
                            ref_gene=ref_gene,
                            row_meta=row_meta,
                        )

                if relaxed_any_hla:
                    relaxed_hla_row_hits += 1

        scan_rows.append({
            "reference_file": str(ref_file),
            "rows_seen": rows_seen,
            "exact_peptide_matched_reference_rows": exact_row_hits,
            "exact_peptide_hla_matched_reference_rows": exact_hla_row_hits,
            "relaxed_min8_peptide_matched_reference_rows": relaxed_row_hits,
            "relaxed_min8_peptide_hla_matched_reference_rows": relaxed_hla_row_hits,
        })

    print("Writing candidate flag table...")
    for col in sorted(flags):
        cand[col] = flags[col]

    cand.to_csv(CANDIDATE_OUT, sep="\t", index=False, compression="gzip")

    print("Writing match examples...")
    pd.DataFrame(examples).to_csv(MATCH_EXAMPLES_OUT, sep="\t", index=False, compression="gzip")

    print("Writing reference scan summary...")
    pd.DataFrame(scan_rows).to_csv(REFERENCE_SCAN_OUT, sep="\t", index=False)

    print("Writing inventory...")
    inv_rows = [
        {"metric": "n_candidates", "value": len(cand)},
        {"metric": "n_samples", "value": cand["sample_id"].nunique()},
        {"metric": "n_unique_peptides", "value": cand["candidate_peptide_clean"].nunique()},
        {
            "metric": "n_unique_peptide_hla_pairs",
            "value": cand[["candidate_peptide_clean", "candidate_hla_normalized_joined"]].drop_duplicates().shape[0],
        },
    ]

    for col in sorted(flags):
        inv_rows.append({
            "metric": f"{col}__candidate_count",
            "value": int(flags[col].sum()),
        })

    pd.DataFrame(inv_rows).to_csv(INVENTORY_OUT, sep="\t", index=False)

    print("Done.")
    print(f"Candidate evidence flags: {CANDIDATE_OUT}")
    print(f"Match examples: {MATCH_EXAMPLES_OUT}")
    print(f"Inventory: {INVENTORY_OUT}")
    print(f"Reference scan summary: {REFERENCE_SCAN_OUT}")


if __name__ == "__main__":
    main()
