#!/usr/bin/env python3

from pathlib import Path
from collections import Counter, defaultdict
import gzip
import csv
import re
import pandas as pd


NORMALIZED_DIR = Path("/data/fsluma/pipelines/Epitope_Generation_Gateway/Immunopeptidomics_Data/normalized_reference_tables")
CANDIDATE_DIR = Path("/data/fsluma/pipelines/Epitope_Generation_Gateway/TCGA_melanoma/epitopes_prioritisation/final_epitopes/HLA")

OUTDIR = Path("unresolved_hla_audit").resolve()
OUTDIR.mkdir(exist_ok=True)


def clean(x):
    if x is None:
        return ""
    x = str(x).strip()
    if x.upper() in {"", "NA", "NAN", "NONE", "NULL"}:
        return ""
    return x


def classify_hla_original(x):
    s = clean(x)
    u = s.upper().replace(" ", "")

    if not u:
        return "missing_original"

    # Non-human MHC examples
    if re.search(r"(^H2|^H-2|H2-|H-2|H2K|H2D|H2IA|H2-IA|H2-IE|RT1|MAMU|PATR|BOLA|DLA|SLA-|SLA\*)", u):
        return "nonhuman_mhc"

    # Explicitly broad/ambiguous
    if re.search(r"(NOTDETERMINED|UNKNOWN|UNDETERMINED|UNRESTRICTED|MHCCLASS|HLACLASS|CLASSI$|CLASSII$|MHC-I|MHC-II)", u):
        return "broad_or_ambiguous"

    # Multiple alternatives or complex entries
    if re.search(r"[/;,]|(\bOR\b)|(\bAND\b)", s, flags=re.I):
        return "multiple_or_complex"

    # Non-classical human class I
    if re.search(r"^(HLA-)?[EFG](\*|:|$)", u):
        return "human_nonclassical_class_I_EFG"

    # Classical human class I
    if re.search(r"^(HLA-)?A(\*|:|\d|$)", u):
        return "human_classical_class_I_A"
    if re.search(r"^(HLA-)?B(\*|:|\d|$)", u):
        return "human_classical_class_I_B"
    if re.search(r"^(HLA-)?C(W)?(\*|:|\d|$)", u):
        return "human_classical_class_I_C"

    # Classical human class II
    if re.search(r"^(HLA-)?DR", u):
        return "human_classical_class_II_DR"
    if re.search(r"^(HLA-)?DQ", u):
        return "human_classical_class_II_DQ"
    if re.search(r"^(HLA-)?DP", u):
        return "human_classical_class_II_DP"

    # Other HLA loci / MHC-related human molecules
    if re.search(r"^(HLA-)?D(MA|MB|OA|OB)", u):
        return "human_other_HLA_DM_DO"
    if re.search(r"^(HLA-)?H(\*|:|$)", u):
        return "human_other_HLA_H"
    if re.search(r"^(MICA|MICB|MR1|CD1A|CD1B|CD1C|CD1D|CD1E)", u):
        return "human_MHC_related_non_HLA_or_like"

    # Anything else containing HLA
    if "HLA" in u:
        return "other_HLA_like"

    return "other_unclassified"


def file_label(path):
    name = path.name.replace(".normalized.tsv.gz", "")
    return name


def audit_reference_tables():
    unresolved_counter = Counter()
    all_original_counter = Counter()
    normalized_counter = Counter()
    class_counter = Counter()
    per_file_category = Counter()
    examples = defaultdict(list)

    files = sorted(NORMALIZED_DIR.glob("*.normalized.tsv.gz"))

    for f in files:
        print(f"Auditing {f.name}")

        with gzip.open(f, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")

            required = {"hla_original", "hla_normalized", "hla_class"}
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise RuntimeError(f"{f} missing required columns: {missing}")

            for row in reader:
                orig = clean(row.get("hla_original", ""))
                norm = clean(row.get("hla_normalized", ""))
                cls = clean(row.get("hla_class", ""))

                if orig:
                    all_original_counter[(file_label(f), orig)] += 1

                if norm:
                    normalized_counter[(file_label(f), norm)] += 1

                if cls:
                    class_counter[(file_label(f), cls)] += 1

                if orig and not norm:
                    cat = classify_hla_original(orig)
                    unresolved_counter[(file_label(f), cat, orig)] += 1
                    per_file_category[(file_label(f), cat)] += 1

                    if len(examples[(file_label(f), cat)]) < 10:
                        examples[(file_label(f), cat)].append({
                            "file": file_label(f),
                            "category": cat,
                            "hla_original": orig,
                            "hla_normalized": norm,
                            "hla_class": cls,
                            "peptide": clean(row.get("peptide", "")),
                            "evidence_class": clean(row.get("evidence_class", "")),
                            "evidence_strength": clean(row.get("evidence_strength", "")),
                            "source_organism": clean(row.get("source_organism", "")),
                            "pmid": clean(row.get("pmid", "")),
                        })

    unresolved_rows = []
    for (f, cat, orig), n in unresolved_counter.items():
        unresolved_rows.append({
            "file": f,
            "category": cat,
            "hla_original": orig,
            "count": n,
        })

    unresolved_df = pd.DataFrame(unresolved_rows)
    if not unresolved_df.empty:
        unresolved_df = unresolved_df.sort_values(["count", "file"], ascending=[False, True])
    unresolved_df.to_csv(OUTDIR / "unresolved_hla_original_values.tsv", sep="\t", index=False)

    cat_rows = []
    for (f, cat), n in per_file_category.items():
        cat_rows.append({
            "file": f,
            "category": cat,
            "count": n,
        })

    cat_df = pd.DataFrame(cat_rows)
    if not cat_df.empty:
        cat_df = cat_df.sort_values(["file", "count"], ascending=[True, False])
    cat_df.to_csv(OUTDIR / "unresolved_hla_category_summary_by_file.tsv", sep="\t", index=False)

    example_rows = []
    for vals in examples.values():
        example_rows.extend(vals)

    example_df = pd.DataFrame(example_rows)
    example_df.to_csv(OUTDIR / "unresolved_hla_examples.tsv", sep="\t", index=False)

    # Human-like unresolved values that may be worth fixing
    patchable_categories = {
        "human_nonclassical_class_I_EFG",
        "human_classical_class_I_A",
        "human_classical_class_I_B",
        "human_classical_class_I_C",
        "human_classical_class_II_DR",
        "human_classical_class_II_DQ",
        "human_classical_class_II_DP",
        "human_other_HLA_DM_DO",
        "human_other_HLA_H",
        "other_HLA_like",
    }

    patchable_df = unresolved_df[unresolved_df["category"].isin(patchable_categories)].copy()
    patchable_df.to_csv(OUTDIR / "potentially_patchable_human_hla_like_values.tsv", sep="\t", index=False)

    return unresolved_df, cat_df, patchable_df


def audit_candidate_hlas():
    files = sorted(CANDIDATE_DIR.glob("Sample_*_epitopes_final.with_hla.csv"))
    rows = []

    if not files:
        return pd.DataFrame()

    for f in files:
        try:
            df = pd.read_csv(f, dtype=str)
        except Exception as e:
            rows.append({
                "sample_file": f.name,
                "hla_allele": "",
                "count": 0,
                "category": "read_error",
                "error": str(e),
            })
            continue

        if "HLA.Allele" not in df.columns:
            continue

        for hla, n in df["HLA.Allele"].fillna("").astype(str).str.strip().value_counts().items():
            if hla:
                rows.append({
                    "sample_file": f.name,
                    "hla_allele": hla,
                    "count": int(n),
                    "category": classify_hla_original(hla),
                    "error": "",
                })

    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["category", "count"], ascending=[True, False])

    out.to_csv(OUTDIR / "candidate_hla_alleles_by_file.tsv", sep="\t", index=False)

    if not out.empty:
        collapsed = (
            out.groupby(["hla_allele", "category"], dropna=False)
            .agg(
                total_candidate_rows=("count", "sum"),
                n_files=("sample_file", "nunique"),
            )
            .reset_index()
            .sort_values("total_candidate_rows", ascending=False)
        )
        collapsed.to_csv(OUTDIR / "candidate_hla_alleles_collapsed.tsv", sep="\t", index=False)

    return out


def main():
    unresolved_df, cat_df, patchable_df = audit_reference_tables()
    cand_df = audit_candidate_hlas()

    print()
    print(f"Wrote audit folder: {OUTDIR}")

    print()
    print("Top unresolved categories by file:")
    if cat_df.empty:
        print("No unresolved HLA original values found.")
    else:
        print(cat_df.head(50).to_string(index=False))

    print()
    print("Top potentially patchable human-like unresolved HLA values:")
    if patchable_df.empty:
        print("No potentially patchable human-like unresolved values found.")
    else:
        print(patchable_df.head(50).to_string(index=False))

    print()
    print("Candidate HLA categories:")
    if cand_df.empty:
        print("No candidate HLA files found or no HLA.Allele column.")
    else:
        c = (
            cand_df.groupby("category")["count"]
            .sum()
            .sort_values(ascending=False)
            .reset_index()
        )
        print(c.to_string(index=False))


if __name__ == "__main__":
    main()
