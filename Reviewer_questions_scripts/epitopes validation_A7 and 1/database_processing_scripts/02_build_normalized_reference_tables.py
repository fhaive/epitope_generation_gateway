#!/usr/bin/env python3

from pathlib import Path
import argparse
import csv
import gzip
import re
from collections import Counter, defaultdict


ROOT = Path("/data/fsluma/pipelines/Epitope_Generation_Gateway/Immunopeptidomics_Data")
OUTDIR = ROOT / "normalized_reference_tables"
OUTDIR.mkdir(exist_ok=True)

AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")

OUT_COLUMNS = [
    "database",
    "evidence_group",
    "evidence_class",
    "evidence_strength",
    "peptide",
    "peptide_length",
    "hla_original",
    "hla_normalized",
    "hla_class",
    "assay_outcome",
    "response_measured",
    "mhc_evidence_code",
    "source_molecule",
    "molecule_parent",
    "source_organism",
    "gene",
    "mutation",
    "cancer_type",
    "mutation_type",
    "pmid",
    "reference_title",
    "journal",
    "year",
    "score_mhcf_rank",
    "score_net4_rank",
    "source_file",
]


IEDB_CEDAR_CONFIGS = [
    {
        "database": "CEDAR",
        "evidence_group": "tcell",
        "path": ROOT / "cedar/tcell_full/tcell_full_v3.csv",
        "out": OUTDIR / "cedar_tcell.normalized.tsv.gz",
        "cols": {
            "pmid": 3,
            "reference_title": 8,
            "peptide": 11,
            "source_molecule": 19,
            "molecule_parent": 21,
            "source_organism": 23,
            "response_measured": 120,
            "assay_outcome": 123,
            "hla_original": 142,
            "mhc_evidence_code": 144,
            "hla_class": 146,
        },
    },
    {
        "database": "CEDAR",
        "evidence_group": "mhc_ligand",
        "path": ROOT / "cedar/mhc_ligand_single_file/mhc_ligand_full.csv",
        "out": OUTDIR / "cedar_mhc_ligand.normalized.tsv.gz",
        "cols": {
            "pmid": 3,
            "reference_title": 8,
            "peptide": 11,
            "source_molecule": 19,
            "molecule_parent": 21,
            "source_organism": 23,
            "response_measured": 92,
            "assay_outcome": 95,
            "hla_original": 108,
            "mhc_evidence_code": 110,
            "hla_class": 112,
        },
    },
    {
        "database": "IEDB",
        "evidence_group": "tcell",
        "path": ROOT / "iedb/tcell_full/tcell_full_v3.csv",
        "out": OUTDIR / "iedb_tcell.normalized.tsv.gz",
        "cols": {
            "pmid": 3,
            "reference_title": 8,
            "peptide": 11,
            "source_molecule": 19,
            "molecule_parent": 21,
            "source_organism": 23,
            "response_measured": 119,
            "assay_outcome": 122,
            "hla_original": 141,
            "mhc_evidence_code": 143,
            "hla_class": 145,
        },
    },
    {
        "database": "IEDB",
        "evidence_group": "mhc_ligand",
        "path": ROOT / "iedb/mhc_ligand_single_file/mhc_ligand_full.csv",
        "out": OUTDIR / "iedb_mhc_ligand.normalized.tsv.gz",
        "cols": {
            "pmid": 3,
            "reference_title": 8,
            "peptide": 11,
            "source_molecule": 19,
            "molecule_parent": 21,
            "source_organism": 23,
            "response_measured": 91,
            "assay_outcome": 94,
            "hla_original": 107,
            "mhc_evidence_code": 109,
            "hla_class": 111,
        },
    },
]


TSNADB_VALIDATED = {
    "database": "TSNAdb_v2",
    "evidence_group": "validated_collected",
    "path": ROOT / "tsnadb_v2/validated_tsnadb2_download/validated_tsnadb2_download.txt",
    "out": OUTDIR / "tsnadb_validated_collected.normalized.tsv.gz",
}


TSNADB_PREDICTED = [
    {
        "database": "TSNAdb_v2",
        "evidence_group": "predicted_snv",
        "path": ROOT / "tsnadb_v2/SNV-derived/SNV-derived.txt",
        "out": OUTDIR / "tsnadb_predicted_snv.normalized.tsv.gz",
    },
    {
        "database": "TSNAdb_v2",
        "evidence_group": "predicted_indel",
        "path": ROOT / "tsnadb_v2/INDEL-derived/INDEL-derived.txt",
        "out": OUTDIR / "tsnadb_predicted_indel.normalized.tsv.gz",
    },
    {
        "database": "TSNAdb_v2",
        "evidence_group": "predicted_fusion",
        "path": ROOT / "tsnadb_v2/Fusion-derived/Fusion-derived.txt",
        "out": OUTDIR / "tsnadb_predicted_fusion.normalized.tsv.gz",
    },
]


def clean_text(x):
    if x is None:
        return ""
    x = str(x).replace("\r", "").replace("\n", " ").strip()
    if x.upper() in {"NA", "N/A", "NONE", "NULL", "NAN"}:
        return ""
    return x


def normalize_peptide(x):
    x = clean_text(x).upper()
    x = re.sub(r"\s+", "", x)
    return x


def valid_peptide(x):
    if not x:
        return False
    if len(x) < 7 or len(x) > 35:
        return False
    return bool(AA_RE.match(x))


def normalize_hla(x):
    """
    Converts common class I and class II HLA strings to a consistent style.

    Examples:
      B44:03       -> HLA-B*44:03
      A02:01       -> HLA-A*02:01
      HLA-A02:01   -> HLA-A*02:01
      HLA-A*02:01  -> HLA-A*02:01
      HLA-A2       -> HLA-A*02
      DRB1*04:01   -> HLA-DRB1*04:01

    Non-human MHC such as H2-Kb is returned as blank for hla_normalized.
    """
    raw = clean_text(x).upper()
    if not raw:
        return ""

    raw = raw.replace(" ", "")
    raw = raw.replace("HLA_", "HLA-")
    raw = raw.replace("HLA", "HLA-", 1) if raw.startswith("HLA") and not raw.startswith("HLA-") else raw

    # Human class I: A, B, C
    m = re.search(r"(?:HLA-)?([ABC])\*?0?(\d{1,2})(?::?(\d{2}))?", raw)
    if m:
        locus = m.group(1)
        group = int(m.group(2))
        allele = m.group(3)
        if allele:
            return f"HLA-{locus}*{group:02d}:{allele}"
        return f"HLA-{locus}*{group:02d}"

    # Human class II
    m = re.search(r"(?:HLA-)?(DRB1|DRB3|DRB4|DRB5|DQB1|DQA1|DPB1|DPA1)\*?0?(\d{1,2})(?::?(\d{2}))?", raw)
    if m:
        locus = m.group(1)
        group = int(m.group(2))
        allele = m.group(3)
        if allele:
            return f"HLA-{locus}*{group:02d}:{allele}"
        return f"HLA-{locus}*{group:02d}"

    return ""


def infer_hla_class(hla_original, hla_normalized, provided_class=""):
    provided_class = clean_text(provided_class)
    if provided_class in {"I", "II"}:
        return provided_class

    h = hla_normalized.upper()
    raw = clean_text(hla_original).upper()

    if re.search(r"HLA-[ABC]\*", h):
        return "I"
    if re.search(r"HLA-(DR|DQ|DP)", h):
        return "II"
    if raw.startswith("H2-") or raw.startswith("H-2"):
        return "nonhuman"
    return ""


def get(row, idx):
    if idx is None:
        return ""
    if idx >= len(row):
        return ""
    return clean_text(row[idx])


def outcome_is_positive(outcome):
    x = clean_text(outcome).lower()
    return x.startswith("positive")


def outcome_is_negative(outcome):
    x = clean_text(outcome).lower()
    return x.startswith("negative")


def response_is_binding_like(response):
    x = clean_text(response).lower()
    return "binding" in x or "affinity" in x or "dissociation constant" in x


def is_ligand_presentation(response, evidence_code):
    r = clean_text(response).lower()
    e = clean_text(evidence_code).lower()

    if "ligand presentation" in r:
        return True
    if "elution" in e:
        return True
    if "mass spectrometry" in r or "mass spectrometry" in e:
        return True
    if re.search(r"\bms\b", r) or re.search(r"\bms\b", e):
        return True
    return False


def classify_iedb_cedar(evidence_group, assay_outcome, response_measured, mhc_evidence_code):
    positive = outcome_is_positive(assay_outcome)
    negative = outcome_is_negative(assay_outcome)
    binding_like = response_is_binding_like(response_measured)
    presentation = is_ligand_presentation(response_measured, mhc_evidence_code)

    if evidence_group == "tcell":
        if positive and not binding_like:
            return "tcell_positive", "strong_immunogenicity"
        if negative and not binding_like:
            return "tcell_negative", "negative_tcell_assay"
        if positive and binding_like:
            return "tcell_file_binding_positive", "binding_evidence"
        if negative and binding_like:
            return "tcell_file_binding_negative", "negative_binding_assay"
        return "tcell_other", "uncategorized_tcell"

    if evidence_group == "mhc_ligand":
        if positive and presentation:
            return "mhc_ligand_presentation_positive", "presentation"
        if negative and presentation:
            return "mhc_ligand_presentation_negative", "negative_presentation_assay"
        if positive and binding_like:
            return "mhc_binding_positive", "binding_evidence"
        if negative and binding_like:
            return "mhc_binding_negative", "negative_binding_assay"
        if presentation:
            return "mhc_ligand_presentation_other", "presentation"
        return "mhc_ligand_other", "uncategorized_mhc_ligand"

    return "other", "other"


def classify_tsnadb_level(level):
    level = clean_text(level).lower()

    if level == "tier1":
        return "tsnadb_validated_tier1", "strong_immunogenicity_and_presentation"
    if level == "tier2":
        return "tsnadb_validated_tier2", "strong_immunogenicity"
    if level == "tier3":
        return "tsnadb_validated_tier3", "presentation"

    return "tsnadb_validated_other", "curated_other"


def empty_record(database, evidence_group, source_file):
    return {
        "database": database,
        "evidence_group": evidence_group,
        "evidence_class": "",
        "evidence_strength": "",
        "peptide": "",
        "peptide_length": "",
        "hla_original": "",
        "hla_normalized": "",
        "hla_class": "",
        "assay_outcome": "",
        "response_measured": "",
        "mhc_evidence_code": "",
        "source_molecule": "",
        "molecule_parent": "",
        "source_organism": "",
        "gene": "",
        "mutation": "",
        "cancer_type": "",
        "mutation_type": "",
        "pmid": "",
        "reference_title": "",
        "journal": "",
        "year": "",
        "score_mhcf_rank": "",
        "score_net4_rank": "",
        "source_file": str(source_file),
    }


def write_header(writer):
    writer.writerow(OUT_COLUMNS)


def write_record(writer, rec):
    writer.writerow([rec.get(c, "") for c in OUT_COLUMNS])


def process_iedb_cedar(config, max_rows=None):
    database = config["database"]
    evidence_group = config["evidence_group"]
    path = config["path"]
    out = config["out"]
    cols = config["cols"]

    print(f"\nProcessing {database} {evidence_group}")
    print(f"Input : {path}")
    print(f"Output: {out}")

    counts = Counter()

    with path.open("r", encoding="utf-8", errors="replace", newline="") as fin, \
         gzip.open(out, "wt", encoding="utf-8", newline="") as fout:

        reader = csv.reader(fin)
        writer = csv.writer(fout, delimiter="\t")
        write_header(writer)

        # IEDB/CEDAR exports use two header rows.
        header1 = next(reader, None)
        header2 = next(reader, None)
        counts["header_rows_skipped"] = 2

        for i, row in enumerate(reader, start=1):
            if max_rows is not None and i > max_rows:
                break

            counts["rows_seen"] += 1

            peptide = normalize_peptide(get(row, cols["peptide"]))
            if not valid_peptide(peptide):
                counts["invalid_or_unusable_peptide"] += 1
                continue

            hla_original = get(row, cols["hla_original"])
            hla_normalized = normalize_hla(hla_original)
            hla_class = infer_hla_class(
                hla_original=hla_original,
                hla_normalized=hla_normalized,
                provided_class=get(row, cols["hla_class"]),
            )

            assay_outcome = get(row, cols["assay_outcome"])
            response_measured = get(row, cols["response_measured"])
            mhc_evidence_code = get(row, cols["mhc_evidence_code"])

            evidence_class, evidence_strength = classify_iedb_cedar(
                evidence_group=evidence_group,
                assay_outcome=assay_outcome,
                response_measured=response_measured,
                mhc_evidence_code=mhc_evidence_code,
            )

            rec = empty_record(database, evidence_group, path)
            rec.update({
                "evidence_class": evidence_class,
                "evidence_strength": evidence_strength,
                "peptide": peptide,
                "peptide_length": len(peptide),
                "hla_original": hla_original,
                "hla_normalized": hla_normalized,
                "hla_class": hla_class,
                "assay_outcome": assay_outcome,
                "response_measured": response_measured,
                "mhc_evidence_code": mhc_evidence_code,
                "source_molecule": get(row, cols["source_molecule"]),
                "molecule_parent": get(row, cols["molecule_parent"]),
                "source_organism": get(row, cols["source_organism"]),
                "pmid": get(row, cols["pmid"]),
                "reference_title": get(row, cols["reference_title"]),
            })

            write_record(writer, rec)

            counts["rows_written"] += 1
            counts[f"evidence_class::{evidence_class}"] += 1
            counts[f"evidence_strength::{evidence_strength}"] += 1
            if hla_normalized:
                counts["rows_with_normalized_hla"] += 1
            if hla_class:
                counts[f"hla_class::{hla_class}"] += 1

    return {
        "database": database,
        "evidence_group": evidence_group,
        "input_file": str(path),
        "output_file": str(out),
        **dict(counts),
    }


def process_tsnadb_validated(config, max_rows=None):
    database = config["database"]
    evidence_group = config["evidence_group"]
    path = config["path"]
    out = config["out"]

    print(f"\nProcessing {database} {evidence_group}")
    print(f"Input : {path}")
    print(f"Output: {out}")

    # The TSNAdb validated file has a 14-column header but 17-column data rows.
    names = [
        "Level",
        "Patient ID",
        "Tumor Type",
        "Tumor Type detaile",
        "Mutation Type",
        "Gene",
        "Mutation",
        "Position",
        "Mutant Peptide",
        "HLA",
        "Pubmed ID",
        "Reference Name",
        "Journal",
        "Year",
        "DOI_or_URL",
        "Source Databases",
        "Source Count",
    ]

    counts = Counter()

    with path.open("r", encoding="utf-8", errors="replace", newline="") as fin, \
         gzip.open(out, "wt", encoding="utf-8", newline="") as fout:

        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t")
        write_header(writer)

        header = next(reader, None)
        counts["header_rows_skipped"] = 1

        for i, row in enumerate(reader, start=1):
            if max_rows is not None and i > max_rows:
                break

            counts["rows_seen"] += 1

            if len(row) < len(names):
                row = row + [""] * (len(names) - len(row))

            d = {names[j]: clean_text(row[j]) for j in range(len(names))}

            peptide = normalize_peptide(d["Mutant Peptide"])
            if not valid_peptide(peptide):
                counts["invalid_or_unusable_peptide"] += 1
                continue

            hla_original = d["HLA"]
            hla_normalized = normalize_hla(hla_original)
            hla_class = infer_hla_class(hla_original, hla_normalized)

            evidence_class, evidence_strength = classify_tsnadb_level(d["Level"])

            rec = empty_record(database, evidence_group, path)
            rec.update({
                "evidence_class": evidence_class,
                "evidence_strength": evidence_strength,
                "peptide": peptide,
                "peptide_length": len(peptide),
                "hla_original": hla_original,
                "hla_normalized": hla_normalized,
                "hla_class": hla_class,
                "gene": d["Gene"],
                "mutation": d["Mutation"],
                "cancer_type": d["Tumor Type detaile"],
                "mutation_type": d["Mutation Type"],
                "pmid": d["Pubmed ID"],
                "reference_title": d["Reference Name"],
                "journal": d["Journal"],
                "year": d["Year"],
            })

            write_record(writer, rec)

            counts["rows_written"] += 1
            counts[f"evidence_class::{evidence_class}"] += 1
            counts[f"evidence_strength::{evidence_strength}"] += 1
            if hla_normalized:
                counts["rows_with_normalized_hla"] += 1
            if hla_class:
                counts[f"hla_class::{hla_class}"] += 1

    return {
        "database": database,
        "evidence_group": evidence_group,
        "input_file": str(path),
        "output_file": str(out),
        **dict(counts),
    }


def find_column(fieldnames, candidates):
    lower_to_original = {f.lower(): f for f in fieldnames}
    for c in candidates:
        if c.lower() in lower_to_original:
            return lower_to_original[c.lower()]
    return None


def process_tsnadb_predicted(config, max_rows=None):
    database = config["database"]
    evidence_group = config["evidence_group"]
    path = config["path"]
    out = config["out"]

    print(f"\nProcessing {database} {evidence_group}")
    print(f"Input : {path}")
    print(f"Output: {out}")

    counts = Counter()

    with path.open("r", encoding="utf-8", errors="replace", newline="") as fin, \
         gzip.open(out, "wt", encoding="utf-8", newline="") as fout:

        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t")
        write_header(writer)

        fieldnames = reader.fieldnames or []

        peptide_col = find_column(fieldnames, ["Peptide"])
        hla_col = find_column(fieldnames, ["HLA"])
        mhcf_col = find_column(fieldnames, ["MHCf_rank (%)", "MHCf_rank"])
        net4_col = find_column(fieldnames, ["Net4_rank (%)", "Net4_rank"])
        gene_col = find_column(fieldnames, ["Gene"])
        mutation_col = find_column(fieldnames, ["Mutation"])
        tumor_col = find_column(fieldnames, ["Tumor Type", "Cancer Type"])
        mutation_type_col = find_column(fieldnames, ["Mutation Type"])

        evidence_class = f"tsnadb_{evidence_group}"
        evidence_strength = "prediction_concordance"

        for i, row in enumerate(reader, start=1):
            if max_rows is not None and i > max_rows:
                break

            counts["rows_seen"] += 1

            peptide = normalize_peptide(row.get(peptide_col, ""))
            if not valid_peptide(peptide):
                counts["invalid_or_unusable_peptide"] += 1
                continue

            hla_original = clean_text(row.get(hla_col, ""))
            hla_normalized = normalize_hla(hla_original)
            hla_class = infer_hla_class(hla_original, hla_normalized)

            rec = empty_record(database, evidence_group, path)
            rec.update({
                "evidence_class": evidence_class,
                "evidence_strength": evidence_strength,
                "peptide": peptide,
                "peptide_length": len(peptide),
                "hla_original": hla_original,
                "hla_normalized": hla_normalized,
                "hla_class": hla_class,
                "gene": clean_text(row.get(gene_col, "")),
                "mutation": clean_text(row.get(mutation_col, "")),
                "cancer_type": clean_text(row.get(tumor_col, "")),
                "mutation_type": clean_text(row.get(mutation_type_col, "")),
                "score_mhcf_rank": clean_text(row.get(mhcf_col, "")),
                "score_net4_rank": clean_text(row.get(net4_col, "")),
            })

            write_record(writer, rec)

            counts["rows_written"] += 1
            counts[f"evidence_class::{evidence_class}"] += 1
            counts[f"evidence_strength::{evidence_strength}"] += 1
            if hla_normalized:
                counts["rows_with_normalized_hla"] += 1
            if hla_class:
                counts[f"hla_class::{hla_class}"] += 1

    return {
        "database": database,
        "evidence_group": evidence_group,
        "input_file": str(path),
        "output_file": str(out),
        **dict(counts),
    }


def write_summary(summary_rows):
    all_keys = set()
    for row in summary_rows:
        all_keys.update(row.keys())

    preferred = [
        "database",
        "evidence_group",
        "input_file",
        "output_file",
        "rows_seen",
        "rows_written",
        "invalid_or_unusable_peptide",
        "rows_with_normalized_hla",
    ]

    other_keys = sorted(k for k in all_keys if k not in preferred)
    columns = preferred + other_keys

    out = OUTDIR / "normalization_summary.tsv"
    with out.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=columns)
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    print(f"\nWrote summary: {out}")


def write_manifest():
    out = OUTDIR / "normalized_reference_manifest.tsv"
    rows = []

    for p in sorted(OUTDIR.glob("*.normalized.tsv.gz")):
        rows.append({
            "file": str(p),
            "size_mb": round(p.stat().st_size / (1024 ** 2), 3),
        })

    with out.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=["file", "size_mb"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote manifest: {out}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--max-rows",
        type=int,
        default=None,
        help="Optional test mode: process only this many data rows per file.",
    )
    args = parser.parse_args()

    summary_rows = []

    for config in IEDB_CEDAR_CONFIGS:
        summary_rows.append(process_iedb_cedar(config, max_rows=args.max_rows))

    summary_rows.append(process_tsnadb_validated(TSNADB_VALIDATED, max_rows=args.max_rows))

    for config in TSNADB_PREDICTED:
        summary_rows.append(process_tsnadb_predicted(config, max_rows=args.max_rows))

    write_summary(summary_rows)
    write_manifest()

    print("\nDone.")


if __name__ == "__main__":
    main()
