#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def savefig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()



NORMALIZED_DIR = Path("/data/fsluma/pipelines/Epitope_Generation_Gateway/Immunopeptidomics_Data/normalized_reference_tables")

OUTDIR = Path("normalized_db_evidence_summary").resolve()
FIGDIR = OUTDIR / "figures"
OUTDIR.mkdir(exist_ok=True)
FIGDIR.mkdir(exist_ok=True)

CHUNKSIZE = 250_000


def clean_series(s):
    s = s.fillna("").astype(str).str.strip()
    s = s.replace({"NA": "", "nan": "", "None": "", "NaN": ""})
    return s


def derive_database_and_group(filename):
    name = filename.replace(".normalized.tsv.gz", "")

    if name.startswith("cedar_"):
        return "CEDAR", name.replace("cedar_", "", 1)

    if name.startswith("iedb_"):
        return "IEDB", name.replace("iedb_", "", 1)

    if name.startswith("tsnadb_"):
        return "TSNAdb_v2", name.replace("tsnadb_", "", 1)

    return "unknown", name


def find_col(columns, candidates, fallback_index=None):
    for c in candidates:
        if c in columns:
            return c

    lower_map = {c.lower(): c for c in columns}
    for c in candidates:
        if c.lower() in lower_map:
            return lower_map[c.lower()]

    if fallback_index is not None and fallback_index < len(columns):
        return columns[fallback_index]

    return None


def detect_columns(path):
    header = pd.read_csv(path, sep="\t", nrows=0, compression="gzip")
    cols = list(header.columns)

    detected = {
        "database": find_col(cols, ["database", "Database"], fallback_index=0),
        "evidence_group": find_col(cols, ["evidence_group", "Evidence_Group"], fallback_index=1),
        "evidence_class": find_col(cols, ["evidence_class", "Evidence_Class"], fallback_index=2),
        "evidence_strength": find_col(cols, ["evidence_strength", "Evidence_Strength"], fallback_index=3),

        # Usually normalized peptide is column 6 in your files, after raw peptide.
        "peptide": find_col(
            cols,
            [
                "peptide_normalized",
                "normalized_peptide",
                "peptide",
                "Peptide",
                "peptide_sequence",
                "Mutant Peptide",
            ],
            fallback_index=5,
        ),

        # In your previous check, raw HLA was column 7 and normalized HLA was column 8.
        "hla": find_col(
            cols,
            [
                "hla_normalized",
                "normalized_hla",
                "hla",
                "HLA",
                "hla_or_mhc",
            ],
            fallback_index=7,
        ),

        "hla_class": find_col(
            cols,
            [
                "hla_class",
                "HLA.Class",
                "mhc_class",
                "MHC Class",
            ],
            fallback_index=None,
        ),

        "reference": find_col(
            cols,
            [
                "reference",
                "Reference",
                "pmid",
                "PMID",
                "Pubmed ID",
                "reference_id",
            ],
            fallback_index=None,
        ),
    }

    return cols, detected


def hash_values(series):
    if len(series) == 0:
        return set()
    return set(pd.util.hash_pandas_object(series.astype(str), index=False).astype("uint64").tolist())


def update_summary(summary, key, sub, peptide_col, hla_col, reference_col):
    rec = summary[key]

    rec["n_rows"] += len(sub)

    peptide = clean_series(sub[peptide_col]) if peptide_col in sub.columns else pd.Series([""] * len(sub))
    hla = clean_series(sub[hla_col]) if hla_col in sub.columns else pd.Series([""] * len(sub))

    peptide_mask = peptide != ""
    hla_mask = hla != ""
    peptide_hla_mask = peptide_mask & hla_mask

    rec["rows_with_peptide"] += int(peptide_mask.sum())
    rec["rows_with_hla"] += int(hla_mask.sum())
    rec["rows_with_peptide_hla"] += int(peptide_hla_mask.sum())

    rec["unique_peptide_hashes"].update(hash_values(peptide[peptide_mask]))
    rec["unique_peptide_hla_hashes"].update(hash_values(peptide[peptide_hla_mask] + "|" + hla[peptide_hla_mask]))

    if reference_col is not None and reference_col in sub.columns:
        ref = clean_series(sub[reference_col])
        ref_mask = ref != ""
        rec["unique_reference_hashes"].update(hash_values(ref[ref_mask]))


def summarize_file(path):
    filename = path.name
    fallback_db, fallback_group = derive_database_and_group(filename)

    cols, detected = detect_columns(path)

    required = ["evidence_class", "evidence_strength", "peptide", "hla"]
    missing = [k for k in required if detected[k] is None]

    if missing:
        raise RuntimeError(
            f"Could not detect required columns in {path}\n"
            f"Missing: {missing}\n"
            f"Columns: {cols}"
        )

    detected_rows.append({
        "file": str(path),
        "database_col": detected["database"],
        "evidence_group_col": detected["evidence_group"],
        "evidence_class_col": detected["evidence_class"],
        "evidence_strength_col": detected["evidence_strength"],
        "peptide_col": detected["peptide"],
        "hla_col": detected["hla"],
        "hla_class_col": detected["hla_class"],
        "reference_col": detected["reference"],
        "all_columns": " | ".join(cols),
    })

    usecols = []
    for k, c in detected.items():
        if c is not None and c not in usecols:
            usecols.append(c)

    for chunk in pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        dtype=str,
        chunksize=CHUNKSIZE,
        usecols=usecols,
        low_memory=False,
    ):
        if detected["database"] is not None and detected["database"] in chunk.columns:
            database = clean_series(chunk[detected["database"]])
            database = database.replace("", fallback_db)
        else:
            database = pd.Series([fallback_db] * len(chunk))

        if detected["evidence_group"] is not None and detected["evidence_group"] in chunk.columns:
            evidence_group = clean_series(chunk[detected["evidence_group"]])
            evidence_group = evidence_group.replace("", fallback_group)
        else:
            evidence_group = pd.Series([fallback_group] * len(chunk))

        evidence_class = clean_series(chunk[detected["evidence_class"]])
        evidence_strength = clean_series(chunk[detected["evidence_strength"]])

        if detected["hla_class"] is not None and detected["hla_class"] in chunk.columns:
            hla_class = clean_series(chunk[detected["hla_class"]])
            hla_class = hla_class.replace("", "unknown")
        else:
            hla_class = pd.Series(["unknown"] * len(chunk))

        tmp = chunk.copy()
        tmp["_database"] = database.values
        tmp["_evidence_group"] = evidence_group.values
        tmp["_evidence_class"] = evidence_class.values
        tmp["_evidence_strength"] = evidence_strength.values
        tmp["_hla_class"] = hla_class.values

        group_cols = [
            "_database",
            "_evidence_group",
            "_evidence_class",
            "_evidence_strength",
            "_hla_class",
        ]

        for key, sub in tmp.groupby(group_cols, dropna=False):
            update_summary(
                summary=summary,
                key=key,
                sub=sub,
                peptide_col=detected["peptide"],
                hla_col=detected["hla"],
                reference_col=detected["reference"],
            )


def make_summary_dataframe(summary):
    rows = []

    for key, rec in summary.items():
        database, evidence_group, evidence_class, evidence_strength, hla_class = key

        rows.append({
            "database": database,
            "evidence_group": evidence_group,
            "evidence_class": evidence_class,
            "evidence_strength": evidence_strength,
            "hla_class": hla_class,
            "n_rows": rec["n_rows"],
            "rows_with_peptide": rec["rows_with_peptide"],
            "rows_with_hla": rec["rows_with_hla"],
            "rows_with_peptide_hla": rec["rows_with_peptide_hla"],
            "n_unique_peptides": len(rec["unique_peptide_hashes"]),
            "n_unique_peptide_hla_pairs": len(rec["unique_peptide_hla_hashes"]),
            "n_unique_references": len(rec["unique_reference_hashes"]),
        })

    df = pd.DataFrame(rows)

    if not df.empty:
        df = df.sort_values(
            ["database", "evidence_group", "evidence_strength", "evidence_class", "hla_class"],
            ascending=True,
        )

    return df


def aggregate_summary(df, group_cols):
    out = (
        df
        .groupby(group_cols, dropna=False)
        .agg(
            n_rows=("n_rows", "sum"),
            rows_with_peptide=("rows_with_peptide", "sum"),
            rows_with_hla=("rows_with_hla", "sum"),
            rows_with_peptide_hla=("rows_with_peptide_hla", "sum"),
            n_unique_peptides_approx=("n_unique_peptides", "sum"),
            n_unique_peptide_hla_pairs_approx=("n_unique_peptide_hla_pairs", "sum"),
            n_unique_references_approx=("n_unique_references", "sum"),
        )
        .reset_index()
    )

    out = out.sort_values("n_rows", ascending=False)
    return out


def plot_stacked(df, index_col, column_col, value_col, out_png, title, ylabel):
    if df.empty:
        return

    mat = (
        df
        .pivot_table(index=index_col, columns=column_col, values=value_col, aggfunc="sum", fill_value=0)
    )

    mat = mat.loc[mat.sum(axis=1).sort_values(ascending=False).index]

    ax = mat.plot(kind="bar", stacked=True, figsize=(max(8, 0.8 * len(mat)), 5))
    ax.set_title(title)
    ax.set_xlabel(index_col)
    ax.set_ylabel(ylabel)
    plt.xticks(rotation=35, ha="right")
    plt.legend(title=column_col, bbox_to_anchor=(1.02, 1), loc="upper left")
    savefig(out_png)


def plot_top_classes(df):
    top = (
        df
        .groupby("evidence_class", dropna=False)["n_rows"]
        .sum()
        .sort_values(ascending=False)
        .head(25)
    )

    plt.figure(figsize=(10, 6))
    top.sort_values().plot(kind="barh")
    plt.xlabel("Rows")
    plt.ylabel("Evidence class")
    plt.title("Top evidence classes across normalized reference databases")
    savefig(FIGDIR / "top_evidence_classes_overall.png")


def write_readme(files, full_summary, by_database, by_strength):
    txt = f"""# Normalized reference database evidence summary

## Input directory

{NORMALIZED_DIR}

## Files summarized

{chr(10).join(str(f) for f in files)}

## Main outputs

- `normalized_reference_evidence_summary_full.tsv`
- `normalized_reference_summary_by_database.tsv`
- `normalized_reference_summary_by_database_and_evidence_strength.tsv`
- `normalized_reference_summary_by_database_and_evidence_class.tsv`
- `normalized_reference_summary_by_database_group_strength_class_hla.tsv`
- `detected_columns.tsv`

## Plot outputs

- `figures/rows_by_database_and_evidence_strength.png`
- `figures/rows_by_database_and_evidence_class.png`
- `figures/rows_by_database_and_hla_class.png`
- `figures/top_evidence_classes_overall.png`

## Notes

`n_unique_peptides` and `n_unique_peptide_hla_pairs` in the full table are exact within each row of the summary table.

In aggregated tables, columns ending in `_approx` are sums across groups, so the same peptide may be counted more than once if it appears in multiple evidence classes or files.
"""

    (OUTDIR / "README.md").write_text(txt)


summary = defaultdict(lambda: {
    "n_rows": 0,
    "rows_with_peptide": 0,
    "rows_with_hla": 0,
    "rows_with_peptide_hla": 0,
    "unique_peptide_hashes": set(),
    "unique_peptide_hla_hashes": set(),
    "unique_reference_hashes": set(),
})

detected_rows = []


def main():
    files = sorted(NORMALIZED_DIR.glob("*.normalized.tsv.gz"))

    if not files:
        raise FileNotFoundError(f"No *.normalized.tsv.gz files found in {NORMALIZED_DIR}")

    print("Summarizing normalized reference tables:")
    for f in files:
        print(f"  {f}")

    for path in files:
        print(f"\nProcessing: {path.name}")
        summarize_file(path)

    full = make_summary_dataframe(summary)

    detected = pd.DataFrame(detected_rows)
    detected.to_csv(OUTDIR / "detected_columns.tsv", sep="\t", index=False)

    full_out = OUTDIR / "normalized_reference_evidence_summary_full.tsv"
    full.to_csv(full_out, sep="\t", index=False)

    by_database = aggregate_summary(full, ["database"])
    by_database.to_csv(OUTDIR / "normalized_reference_summary_by_database.tsv", sep="\t", index=False)

    by_strength = aggregate_summary(full, ["database", "evidence_strength"])
    by_strength.to_csv(OUTDIR / "normalized_reference_summary_by_database_and_evidence_strength.tsv", sep="\t", index=False)

    by_class = aggregate_summary(full, ["database", "evidence_class"])
    by_class.to_csv(OUTDIR / "normalized_reference_summary_by_database_and_evidence_class.tsv", sep="\t", index=False)

    by_group_strength_class_hla = aggregate_summary(
        full,
        ["database", "evidence_group", "evidence_strength", "evidence_class", "hla_class"],
    )
    by_group_strength_class_hla.to_csv(
        OUTDIR / "normalized_reference_summary_by_database_group_strength_class_hla.tsv",
        sep="\t",
        index=False,
    )

    plot_stacked(
        by_strength,
        index_col="database",
        column_col="evidence_strength",
        value_col="n_rows",
        out_png=FIGDIR / "rows_by_database_and_evidence_strength.png",
        title="Rows by database and evidence strength",
        ylabel="Rows",
    )

    plot_stacked(
        by_class,
        index_col="database",
        column_col="evidence_class",
        value_col="n_rows",
        out_png=FIGDIR / "rows_by_database_and_evidence_class.png",
        title="Rows by database and evidence class",
        ylabel="Rows",
    )

    by_hla = aggregate_summary(full, ["database", "hla_class"])
    by_hla.to_csv(OUTDIR / "normalized_reference_summary_by_database_and_hla_class.tsv", sep="\t", index=False)

    plot_stacked(
        by_hla,
        index_col="database",
        column_col="hla_class",
        value_col="n_rows",
        out_png=FIGDIR / "rows_by_database_and_hla_class.png",
        title="Rows by database and HLA class",
        ylabel="Rows",
    )

    plot_top_classes(full)

    write_readme(files, full, by_database, by_strength)

    print("\nDone.")
    print(f"Output folder: {OUTDIR}")
    print("\nMain table:")
    print(full_out)

    print("\nQuick database summary:")
    print(by_database.to_string(index=False))


if __name__ == "__main__":
    main()
