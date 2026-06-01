#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
from tqdm import tqdm

BASE = Path("gdc_tcga_skcm_star_counts")
MANIFEST = BASE / "gdc_tcga_skcm_star_counts_manifest.tsv"

OUT = Path("tcga_skcm_expression_matrix")
OUT.mkdir(exist_ok=True)

COUNT_MATRIX_SELECTED = OUT / "tcga_skcm_star_counts_unstranded_primary_tumor_one_sample_per_patient.tsv"
SELECTED_MANIFEST = OUT / "tcga_skcm_primary_tumor_selected_one_sample_per_patient_manifest.tsv"
GENE_ANNOTATION = OUT / "tcga_skcm_star_counts_gene_annotation.tsv"

ORIGINAL_30_SAMPLES = [
    "Sample_TCGA-BF-A1PV",
    "Sample_TCGA-BF-A3DJ",
    "Sample_TCGA-BF-A3DL",
    "Sample_TCGA-BF-A3DM",
    "Sample_TCGA-BF-A3DN",
    "Sample_TCGA-BF-A5EO",
    "Sample_TCGA-BF-A5EP",
    "Sample_TCGA-BF-A5EQ",
    "Sample_TCGA-BF-A5ER",
    "Sample_TCGA-BF-A5ES",
    "Sample_TCGA-BF-AAOU",
    "Sample_TCGA-BF-AAOX",
    "Sample_TCGA-BF-AAP1",
    "Sample_TCGA-BF-AAP2",
    "Sample_TCGA-BF-AAP4",
    "Sample_TCGA-BF-AAP6",
    "Sample_TCGA-BF-AAP7",
    "Sample_TCGA-BF-AAP8",
    "Sample_TCGA-D3-A5GT",
    "Sample_TCGA-D9-A3Z4",
    "Sample_TCGA-D9-A4Z2",
    "Sample_TCGA-D9-A4Z3",
    "Sample_TCGA-EB-A1NK",
    "Sample_TCGA-EB-A3HV",
    "Sample_TCGA-EB-A3XB",
    "Sample_TCGA-EB-A3XC",
    "Sample_TCGA-EB-A3XD",
    "Sample_TCGA-EB-A3XE",
    "Sample_TCGA-EB-A3XF",
    "Sample_TCGA-EB-A41A",
]

ORIGINAL_30_CASES = [x.replace("Sample_", "") for x in ORIGINAL_30_SAMPLES]


def read_star_counts(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)

    required = {"gene_id", "gene_name", "gene_type", "unstranded"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {missing}; columns={list(df.columns)}")

    df = df[~df["gene_id"].astype(str).str.startswith("__")].copy()
    df["gene_id_clean"] = df["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)

    out = df[["gene_id_clean", "gene_id", "gene_name", "gene_type", "unstranded"]].copy()
    out["unstranded"] = pd.to_numeric(out["unstranded"], errors="coerce").fillna(0).astype(int)

    return out


def sample_col_name(case_submitter_id: str) -> str:
    return f"Sample_{case_submitter_id}"


def priority(row) -> tuple:
    case_id = str(row.get("case_submitter_id", ""))

    # Keep original 30 cases first if they are available as primary tumors.
    original_priority = 0 if case_id in ORIGINAL_30_CASES else 1

    return (
        original_priority,
        str(row.get("sample_submitter_id", "")),
        str(row.get("file_name", "")),
    )


def main():
    manifest = pd.read_csv(MANIFEST, sep="\t")

    if "local_path" not in manifest.columns:
        raise RuntimeError("Manifest lacks local_path column. Did download script finish?")

    # PRIMARY TUMOR ONLY
    selected_type = "Primary Tumor"
    manifest_primary = manifest[manifest["sample_type"] == selected_type].copy()

    print("Files by sample_type in full manifest:")
    print(manifest["sample_type"].value_counts(dropna=False).to_string())

    print("\nPrimary Tumor files retained:")
    print(manifest_primary["sample_type"].value_counts(dropna=False).to_string())

    if manifest_primary.empty:
        raise RuntimeError("No Primary Tumor files found in manifest.")

    # Choose one primary tumor RNA-seq file per case.
    manifest_primary["priority_tuple"] = manifest_primary.apply(priority, axis=1)
    manifest_primary = manifest_primary.sort_values("priority_tuple")
    selected = manifest_primary.drop_duplicates("case_submitter_id", keep="first").copy()
    selected = selected.sort_values("case_submitter_id")
    selected["matrix_col"] = selected["case_submitter_id"].map(sample_col_name)

    selected.to_csv(SELECTED_MANIFEST, sep="\t", index=False)

    print(f"\nSelected one Primary Tumor sample per case: {selected.shape[0]}")

    missing_original = sorted(set(ORIGINAL_30_CASES) - set(selected["case_submitter_id"]))
    present_original = sorted(set(ORIGINAL_30_CASES) & set(selected["case_submitter_id"]))

    print(f"\nOriginal 30 target cases present as Primary Tumor: {len(present_original)} / 30")

    if present_original:
        print("Present original target cases:")
        for x in present_original:
            print("  ", x)

    if missing_original:
        print("\nOriginal target cases missing from Primary Tumor matrix:")
        for x in missing_original:
            print("  ", x)

    gene_annot = None
    count_series = []

    for _, row in tqdm(selected.iterrows(), total=selected.shape[0], desc="reading STAR count files"):
        path = row["local_path"]
        col = row["matrix_col"]

        x = read_star_counts(path)

        if gene_annot is None:
            gene_annot = x[["gene_id_clean", "gene_id", "gene_name", "gene_type"]].drop_duplicates("gene_id_clean")

        s = x.set_index("gene_id_clean")["unstranded"]
        s.name = col
        count_series.append(s)

    mat = pd.concat(count_series, axis=1)
    mat = mat.fillna(0).astype(int)
    mat = mat.sort_index()

    gene_annot = gene_annot.drop_duplicates("gene_id_clean").set_index("gene_id_clean")
    gene_annot = gene_annot.reindex(mat.index).reset_index()
    gene_annot.to_csv(GENE_ANNOTATION, sep="\t", index=False)

    mat.to_csv(COUNT_MATRIX_SELECTED, sep="\t", index_label="gene_id")

    print("\nCount matrix written:")
    print(COUNT_MATRIX_SELECTED)
    print("Matrix shape:", mat.shape[0], "genes x", mat.shape[1], "samples")
    print("Selected manifest:", SELECTED_MANIFEST)
    print("Gene annotation:", GENE_ANNOTATION)


if __name__ == "__main__":
    main()
