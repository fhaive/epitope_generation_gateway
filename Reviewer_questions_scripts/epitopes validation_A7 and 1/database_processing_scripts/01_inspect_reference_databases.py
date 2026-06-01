#!/usr/bin/env python3

from pathlib import Path
import csv
import re
import subprocess
import pandas as pd


ROOT = Path("/data/fsluma/pipelines/Epitope_Generation_Gateway/Immunopeptidomics_Data")
OUTDIR = ROOT / "evidence_db_audit"
PREVIEW_DIR = OUTDIR / "previews"
RAW_HEAD_DIR = OUTDIR / "raw_heads"

CEDAR_ROOT = ROOT / "cedar"
IEDB_ROOT = ROOT / "iedb"
TSNADB_ROOT = ROOT / "tsnadb_v2"

OUTDIR.mkdir(exist_ok=True)
PREVIEW_DIR.mkdir(exist_ok=True)
RAW_HEAD_DIR.mkdir(exist_ok=True)


REFERENCE_FILES = [
    {
        "database": "CEDAR",
        "evidence_group": "tcell_full",
        "file": CEDAR_ROOT / "tcell_full" / "tcell_full_v3.csv",
    },
    {
        "database": "CEDAR",
        "evidence_group": "mhc_ligand_single",
        "file": CEDAR_ROOT / "mhc_ligand_single_file" / "mhc_ligand_full.csv",
    },
    {
        "database": "IEDB",
        "evidence_group": "tcell_full",
        "file": IEDB_ROOT / "tcell_full" / "tcell_full_v3.csv",
    },
    {
        "database": "IEDB",
        "evidence_group": "mhc_ligand_single",
        "file": IEDB_ROOT / "mhc_ligand_single_file" / "mhc_ligand_full.csv",
    },
    {
        "database": "TSNAdb_v2",
        "evidence_group": "validated_collected",
        "file": TSNADB_ROOT / "validated_tsnadb2_download" / "validated_tsnadb2_download.txt",
    },
    {
        "database": "TSNAdb_v2",
        "evidence_group": "predicted_snv",
        "file": TSNADB_ROOT / "SNV-derived" / "SNV-derived.txt",
    },
    {
        "database": "TSNAdb_v2",
        "evidence_group": "predicted_indel",
        "file": TSNADB_ROOT / "INDEL-derived" / "INDEL-derived.txt",
    },
    {
        "database": "TSNAdb_v2",
        "evidence_group": "predicted_fusion",
        "file": TSNADB_ROOT / "Fusion-derived" / "Fusion-derived.txt",
    },
]


KEYWORD_PATTERNS = {
    "peptide_sequence": [
        r"peptide",
        r"epitope",
        r"sequence",
        r"seq",
        r"mutant",
        r"neo",
    ],
    "hla_or_mhc": [
        r"hla",
        r"mhc",
        r"allele",
        r"restriction",
    ],
    "assay": [
        r"assay",
        r"method",
        r"technique",
        r"elution",
        r"mass.?spec",
        r"ms",
    ],
    "outcome": [
        r"outcome",
        r"response",
        r"positive",
        r"negative",
        r"qualitative",
        r"measurement",
    ],
    "cancer_or_disease": [
        r"cancer",
        r"tumou?r",
        r"disease",
        r"melanoma",
        r"carcinoma",
        r"neoplasm",
    ],
    "gene_or_antigen": [
        r"gene",
        r"protein",
        r"antigen",
        r"source",
        r"parent",
    ],
    "reference": [
        r"pmid",
        r"pubmed",
        r"reference",
        r"doi",
        r"article",
    ],
    "evidence_tier": [
        r"tier",
        r"level",
        r"evidence",
        r"validated",
    ],
    "tcga": [
        r"tcga",
        r"sample",
        r"barcode",
        r"case",
    ],
}


def safe_name(database, evidence_group):
    return f"{database}__{evidence_group}".replace("/", "_").replace(" ", "_")


def count_lines(path):
    if not path.exists():
        return None
    try:
        result = subprocess.run(
            ["wc", "-l", str(path)],
            check=True,
            capture_output=True,
            text=True,
        )
        return int(result.stdout.strip().split()[0])
    except Exception:
        return None


def get_file_size_mb(path):
    if not path.exists():
        return None
    return round(path.stat().st_size / (1024 ** 2), 3)


def read_raw_head(path, n=8):
    lines = []
    if not path.exists():
        return lines
    with path.open("r", errors="replace") as f:
        for _ in range(n):
            line = f.readline()
            if not line:
                break
            lines.append(line.rstrip("\n"))
    return lines


def detect_delimiter(path):
    raw = "\n".join(read_raw_head(path, n=5))
    if not raw:
        return ","

    try:
        dialect = csv.Sniffer().sniff(raw, delimiters=[",", "\t", ";", "|"])
        return dialect.delimiter
    except Exception:
        suffix = path.suffix.lower()
        if suffix == ".txt":
            return "\t"
        return ","


def flatten_columns(columns):
    flattened = []
    for col in columns:
        if isinstance(col, tuple):
            pieces = [
                str(x).strip()
                for x in col
                if str(x).strip()
                and not str(x).startswith("Unnamed:")
                and str(x).lower() != "nan"
            ]
            flattened.append(" | ".join(pieces))
        else:
            flattened.append(str(col).strip())
    return flattened


def read_preview(path, delimiter, nrows=20):
    """
    Read small preview robustly.

    IEDB/CEDAR files may sometimes have multi-row style headers.
    We save both a standard header parse and, when possible, a 2-row header parse.
    """
    previews = {}

    try:
        df0 = pd.read_csv(
            path,
            sep=delimiter,
            nrows=nrows,
            dtype=str,
            low_memory=False,
            on_bad_lines="skip",
        )
        df0.columns = flatten_columns(df0.columns)
        previews["header_0"] = df0
    except Exception as e:
        previews["header_0_error"] = str(e)

    try:
        df_multi = pd.read_csv(
            path,
            sep=delimiter,
            header=[0, 1],
            nrows=nrows,
            dtype=str,
            low_memory=False,
            on_bad_lines="skip",
        )
        df_multi.columns = flatten_columns(df_multi.columns)
        previews["header_0_1"] = df_multi
    except Exception as e:
        previews["header_0_1_error"] = str(e)

    return previews


def choose_best_preview(previews):
    """
    Prefer the parse with more informative column names.
    """
    candidates = []
    for key in ["header_0", "header_0_1"]:
        df = previews.get(key)
        if isinstance(df, pd.DataFrame):
            cols = list(df.columns)
            unnamed = sum(c.lower().startswith("unnamed") for c in cols)
            emptyish = sum(c.strip() == "" for c in cols)
            unique_cols = len(set(cols))
            score = unique_cols - unnamed - emptyish
            candidates.append((score, key, df))

    if not candidates:
        return None, None

    candidates.sort(reverse=True, key=lambda x: x[0])
    return candidates[0][1], candidates[0][2]


def classify_candidate_columns(columns):
    rows = []
    for col in columns:
        col_lower = col.lower()
        matched_categories = []
        for category, patterns in KEYWORD_PATTERNS.items():
            if any(re.search(p, col_lower) for p in patterns):
                matched_categories.append(category)

        if matched_categories:
            rows.append(
                {
                    "column": col,
                    "candidate_for": ";".join(matched_categories),
                }
            )
    return rows


def main():
    inventory_rows = []
    columns_rows = []
    candidate_rows = []

    for item in REFERENCE_FILES:
        db = item["database"]
        group = item["evidence_group"]
        path = item["file"]
        label = safe_name(db, group)

        print(f"\n=== {db} / {group} ===")
        print(path)

        exists = path.exists()
        delimiter = detect_delimiter(path) if exists else None
        line_count = count_lines(path) if exists else None
        size_mb = get_file_size_mb(path) if exists else None

        raw_head = read_raw_head(path, n=10) if exists else []
        raw_head_file = RAW_HEAD_DIR / f"{label}.raw_head.txt"
        raw_head_file.write_text("\n".join(raw_head), encoding="utf-8")

        best_parse = None
        n_columns = None
        preview_file = None
        best_df = None

        if exists:
            previews = read_preview(path, delimiter=delimiter, nrows=25)

            for parse_name, obj in previews.items():
                if isinstance(obj, pd.DataFrame):
                    out_preview = PREVIEW_DIR / f"{label}.{parse_name}.preview.tsv"
                    obj.head(25).to_csv(out_preview, sep="\t", index=False)

            best_parse, best_df = choose_best_preview(previews)

            if best_df is not None:
                n_columns = best_df.shape[1]
                preview_file = PREVIEW_DIR / f"{label}.{best_parse}.preview.tsv"

                for idx, col in enumerate(best_df.columns, start=1):
                    columns_rows.append(
                        {
                            "database": db,
                            "evidence_group": group,
                            "file": str(path),
                            "best_parse": best_parse,
                            "column_index": idx,
                            "column": col,
                        }
                    )

                for row in classify_candidate_columns(best_df.columns):
                    row.update(
                        {
                            "database": db,
                            "evidence_group": group,
                            "file": str(path),
                            "best_parse": best_parse,
                        }
                    )
                    candidate_rows.append(row)

        inventory_rows.append(
            {
                "database": db,
                "evidence_group": group,
                "file": str(path),
                "exists": exists,
                "delimiter": repr(delimiter),
                "line_count_including_header": line_count,
                "approx_rows_excluding_header": None if line_count is None else max(line_count - 1, 0),
                "size_mb": size_mb,
                "best_parse": best_parse,
                "n_columns_in_best_parse": n_columns,
                "raw_head_file": str(raw_head_file),
                "preview_file": str(preview_file) if preview_file else None,
            }
        )

        print(f"exists: {exists}")
        print(f"delimiter: {repr(delimiter)}")
        print(f"lines: {line_count}")
        print(f"size_mb: {size_mb}")
        print(f"best_parse: {best_parse}")
        print(f"n_columns: {n_columns}")

        if best_df is not None:
            print("first 15 columns:")
            for c in list(best_df.columns)[:15]:
                print(f"  - {c}")

    inventory = pd.DataFrame(inventory_rows)
    columns_long = pd.DataFrame(columns_rows)
    candidate_cols = pd.DataFrame(candidate_rows)

    inventory.to_csv(OUTDIR / "file_inventory.tsv", sep="\t", index=False)
    columns_long.to_csv(OUTDIR / "columns_long.tsv", sep="\t", index=False)
    candidate_cols.to_csv(OUTDIR / "candidate_columns.tsv", sep="\t", index=False)

    print("\n\nWrote audit files:")
    print(OUTDIR / "file_inventory.tsv")
    print(OUTDIR / "columns_long.tsv")
    print(OUTDIR / "candidate_columns.tsv")
    print(PREVIEW_DIR)
    print(RAW_HEAD_DIR)


if __name__ == "__main__":
    main()
