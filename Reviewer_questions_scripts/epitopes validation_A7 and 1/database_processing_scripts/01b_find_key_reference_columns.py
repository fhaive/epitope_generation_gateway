#!/usr/bin/env python3

from pathlib import Path
import re
import pandas as pd

ROOT = Path("/data/fsluma/pipelines/Epitope_Generation_Gateway/Immunopeptidomics_Data")
OUTDIR = ROOT / "evidence_db_audit"
OUTDIR.mkdir(exist_ok=True)

FILES = [
    ("CEDAR", "tcell_full", ROOT / "cedar/tcell_full/tcell_full_v3.csv", ",", [0, 1]),
    ("CEDAR", "mhc_ligand", ROOT / "cedar/mhc_ligand_single_file/mhc_ligand_full.csv", ",", [0, 1]),
    ("IEDB", "tcell_full", ROOT / "iedb/tcell_full/tcell_full_v3.csv", ",", [0, 1]),
    ("IEDB", "mhc_ligand", ROOT / "iedb/mhc_ligand_single_file/mhc_ligand_full.csv", ",", [0, 1]),
    ("TSNAdb_v2", "validated_collected", ROOT / "tsnadb_v2/validated_tsnadb2_download/validated_tsnadb2_download.txt", "\t", 0),
    ("TSNAdb_v2", "predicted_snv", ROOT / "tsnadb_v2/SNV-derived/SNV-derived.txt", "\t", 0),
    ("TSNAdb_v2", "predicted_indel", ROOT / "tsnadb_v2/INDEL-derived/INDEL-derived.txt", "\t", 0),
    ("TSNAdb_v2", "predicted_fusion", ROOT / "tsnadb_v2/Fusion-derived/Fusion-derived.txt", "\t", 0),
]

KEY_PATTERNS = {
    "peptide": r"epitope|peptide|sequence|mutant peptide",
    "hla_mhc": r"hla|mhc|restriction|allele",
    "assay": r"assay|method|technique|elution|mass|ms",
    "outcome": r"qualitative|outcome|positive|negative|measurement|response",
    "disease_tumor": r"disease|tumou?r|cancer|melanoma|cell type",
    "gene_antigen": r"gene|antigen|protein|source|parent|organism",
    "reference": r"pmid|pubmed|reference|doi|article",
    "evidence_level": r"tier|level|evidence|validated",
    "tcga_or_sample": r"tcga|sample|barcode|case|frequency",
    "binding_prediction": r"rank|ic50|affinity|mhcf|netmhc|score",
}


def flatten_columns(columns):
    out = []
    for col in columns:
        if isinstance(col, tuple):
            parts = []
            for x in col:
                x = str(x).strip()
                if not x or x.lower() == "nan" or x.startswith("Unnamed:"):
                    continue
                parts.append(x)
            out.append(" | ".join(parts))
        else:
            out.append(str(col).strip())
    return out


def classify(col):
    hits = []
    x = col.lower()
    for label, pat in KEY_PATTERNS.items():
        if re.search(pat, x, flags=re.I):
            hits.append(label)
    return ";".join(hits)


rows = []
preview_rows = []

for db, group, path, sep, header in FILES:
    print(f"\n=== {db} / {group} ===")
    print(path)

    df = pd.read_csv(
        path,
        sep=sep,
        header=header,
        nrows=8,
        dtype=str,
        low_memory=False,
        on_bad_lines="skip",
    )
    df.columns = flatten_columns(df.columns)

    for i, col in enumerate(df.columns):
        label = classify(col)
        if label:
            values = (
                df[col]
                .dropna()
                .astype(str)
                .head(5)
                .tolist()
            )
            rows.append({
                "database": db,
                "evidence_group": group,
                "column_index_0based": i,
                "column": col,
                "candidate_for": label,
                "example_values": " || ".join(values),
                "file": str(path),
            })

    # save small full preview with readable columns
    preview_file = OUTDIR / f"{db}__{group}.readable_header_preview.tsv"
    df.to_csv(preview_file, sep="\t", index=False)

out = pd.DataFrame(rows)
out_file = OUTDIR / "key_reference_column_candidates.tsv"
out.to_csv(out_file, sep="\t", index=False)

print("\nWrote:")
print(out_file)
print("\nReadable previews written to:")
print(OUTDIR / "*readable_header_preview.tsv")
