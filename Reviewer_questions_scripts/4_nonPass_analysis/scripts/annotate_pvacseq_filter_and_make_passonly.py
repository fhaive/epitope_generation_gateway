#!/usr/bin/env python3

import argparse
import gzip
import shutil
from pathlib import Path

import pandas as pd


PASS_EQUIVALENT = {"PASS", "."}
DEFAULT_CLASSES = ["MHC_Class_I", "MHC_Class_II"]


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def norm_chrom(chrom):
    chrom = str(chrom)
    return chrom[3:] if chrom.startswith("chr") else chrom


def is_pass_filter(value):
    if pd.isna(value):
        return False

    value = str(value).strip()

    if value == "":
        return False

    return value in PASS_EQUIVALENT


def read_vcf_filters(vcf_path):
    records = []

    with open_text(vcf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom, pos, variant_id, ref, alts, qual, filt, info = fields[:8]
            pos = int(pos)

            for alt in alts.split(","):
                records.append({
                    "chrom_norm": norm_chrom(chrom),
                    "vcf_pos": pos,
                    "Reference": ref,
                    "Variant": alt,
                    "VCF_FILTER": filt,
                })

    if not records:
        return pd.DataFrame(columns=["chrom_norm", "vcf_pos", "Reference", "Variant", "VCF_FILTER"])

    vcf_df = pd.DataFrame(records)

    key_cols = ["chrom_norm", "vcf_pos", "Reference", "Variant"]

    vcf_df = (
        vcf_df.groupby(key_cols, as_index=False)
        .agg(VCF_FILTER=("VCF_FILTER", lambda x: ";".join(sorted(set(map(str, x))))))
    )

    return vcf_df


def infer_sample_from_pvac_path(path):
    name = path.name

    if "_CancerDNA" in name:
        return name.split("_CancerDNA")[0]

    return path.parents[1].name


def infer_class_from_pvac_path(path):
    return path.parent.name


def find_vcf_for_sample(vcf_dir, sample):
    candidates = [
        vcf_dir / f"{sample}_somatic_quantification.full_rescue.vcf.gz",
        vcf_dir / f"{sample}_somatic_quantification.full_rescue.vcf",
        vcf_dir / f"{sample}_somatic_quantification.vcf.gz",
        vcf_dir / f"{sample}_somatic_quantification.vcf",
        vcf_dir / f"{sample}_somatic_VEP.full_rescue.vcf.gz",
        vcf_dir / f"{sample}_somatic_VEP.full_rescue.vcf",
        vcf_dir / f"{sample}_somatic_VEP.vcf.gz",
        vcf_dir / f"{sample}_somatic_VEP.vcf",
    ]

    for candidate in candidates:
        if candidate.exists():
            return candidate

    raise FileNotFoundError(
        f"No VCF found for sample {sample}. Tried: "
        + ", ".join(str(x) for x in candidates)
    )


def try_merge_with_coordinate_modes(pvac, vcf_filters):
    key_cols = ["chrom_norm", "vcf_pos", "Reference", "Variant"]

    pvac_start = pvac.copy()
    pvac_start["vcf_pos"] = pd.to_numeric(pvac_start["Start"], errors="coerce").astype("Int64")

    merged_start = pvac_start.merge(vcf_filters, how="left", on=key_cols)
    n_start_matches = int(merged_start["VCF_FILTER"].notna().sum())

    pvac_start_plus_one = pvac.copy()
    pvac_start_plus_one["vcf_pos"] = (
        pd.to_numeric(pvac_start_plus_one["Start"], errors="coerce") + 1
    ).astype("Int64")

    merged_start_plus_one = pvac_start_plus_one.merge(vcf_filters, how="left", on=key_cols)
    n_start_plus_one_matches = int(merged_start_plus_one["VCF_FILTER"].notna().sum())

    if n_start_plus_one_matches > n_start_matches:
        return merged_start_plus_one, "Start_plus_1", n_start_plus_one_matches

    return merged_start, "Start_as_VCF_POS", n_start_matches


def fill_unmatched_with_rowwise_alternative_exact_matches(merged, vcf_filters):
    """
    Recover remaining unmatched pVACseq rows by trying row-wise exact CHROM/POS/REF/ALT
    matches at nearby coordinate representations.

    This is used after selecting the best global coordinate mode. It prevents a small
    number of indel/VEP representation edge cases from being incorrectly classified
    as unmatched.
    """

    merged = merged.copy()

    if "rowwise_alternative_position_match" not in merged.columns:
        merged["rowwise_alternative_position_match"] = False

    # Fast exact lookup: chrom, pos, ref, alt -> filter
    lookup = {}
    for _, r in vcf_filters.iterrows():
        key = (
            str(r["chrom_norm"]),
            int(r["vcf_pos"]),
            str(r["Reference"]),
            str(r["Variant"]),
        )
        lookup[key] = str(r["VCF_FILTER"])

    def as_int(x):
        try:
            if pd.isna(x) or str(x).strip() == "":
                return None
            return int(float(x))
        except Exception:
            return None

    alternative_modes = [
        ("Start_minus_2", "Start", -2),
        ("Start_minus_1", "Start", -1),
        ("Start", "Start", 0),
        ("Start_plus_1", "Start", 1),
        ("Start_plus_2", "Start", 2),
        ("Stop_minus_1", "Stop", -1),
        ("Stop", "Stop", 0),
        ("Stop_plus_1", "Stop", 1),
    ]

    unmatched_idx = merged.index[merged["VCF_FILTER"].isna()].tolist()

    for idx in unmatched_idx:
        row = merged.loc[idx]

        chrom = str(row["chrom_norm"])
        ref = str(row["Reference"])
        alt = str(row["Variant"])

        for mode_name, base_col, offset in alternative_modes:
            base_pos = as_int(row.get(base_col, None))
            if base_pos is None:
                continue

            pos = base_pos + offset
            key = (chrom, pos, ref, alt)

            if key in lookup:
                merged.at[idx, "vcf_pos"] = pos
                merged.at[idx, "VCF_FILTER"] = lookup[key]
                merged.at[idx, "coordinate_match_mode"] = f"rowwise_{mode_name}"
                merged.at[idx, "rowwise_alternative_position_match"] = True
                break

    return merged


def find_pvac_files(pvac_dir, mhc_classes, report_name):
    pvac_files = []

    for mhc_class in mhc_classes:
        pattern = f"Sample_*/{mhc_class}/*{report_name}"
        pvac_files.extend(sorted(pvac_dir.glob(pattern)))

    return sorted(pvac_files)


def find_peptide_column(df):
    candidates = [
        "MT Epitope Seq",
        "MT.Epitope.Seq",
        "Epitope Seq",
        "Epitope.Seq",
        "MT Epitope",
        "Peptide",
        "Mutation",
    ]

    for col in candidates:
        if col in df.columns:
            return col

    return None


def summarize_one(merged, sample, mhc_class, pvac_tsv, vcf_path, coordinate_mode):
    total_rows = len(merged)
    matched_rows = int(merged["VCF_FILTER"].notna().sum())
    unmatched_rows = int(merged["VCF_FILTER"].isna().sum())
    pass_rows = int(merged["would_pass_pass_only"].astype(int).sum())
    nonpass_rows = int(merged["is_nonpass"].astype(int).sum())

    variant_key_cols = ["sample", "chrom_norm", "vcf_pos", "Reference", "Variant"]

    matched = merged[merged["VCF_FILTER"].notna()].copy()
    pass_df = matched[matched["would_pass_pass_only"]].copy()
    nonpass_df = matched[~matched["would_pass_pass_only"]].copy()

    peptide_col = find_peptide_column(merged)

    row = {
        "sample": sample,
        "mhc_class": mhc_class,
        "pvac_tsv": str(pvac_tsv),
        "vcf": str(vcf_path),
        "coordinate_match_mode": coordinate_mode,
        "pvacseq_rows_total": total_rows,
        "pvacseq_rows_matched_to_vcf": matched_rows,
        "pvacseq_rows_unmatched_to_vcf": unmatched_rows,
        "pvacseq_rows_pass": pass_rows,
        "pvacseq_rows_nonpass": nonpass_rows,
        "pvacseq_rows_pass_pct_of_matched": 100 * pass_rows / matched_rows if matched_rows else pd.NA,
        "pvacseq_rows_nonpass_pct_of_matched": 100 * nonpass_rows / matched_rows if matched_rows else pd.NA,
        "pvacseq_rows_unmatched_pct_of_total": 100 * unmatched_rows / total_rows if total_rows else pd.NA,
        "unique_variant_keys_total_in_pvacseq": int(merged[variant_key_cols].drop_duplicates().shape[0]),
        "unique_variant_keys_matched_to_vcf": int(matched[variant_key_cols].drop_duplicates().shape[0]),
        "unique_variant_keys_pass": int(pass_df[variant_key_cols].drop_duplicates().shape[0]),
        "unique_variant_keys_nonpass": int(nonpass_df[variant_key_cols].drop_duplicates().shape[0]),
        "peptide_column_used": peptide_col or "",
    }

    if peptide_col:
        peptide_key_cols = ["sample", peptide_col]
        row.update({
            "unique_peptides_total": int(merged[peptide_key_cols].drop_duplicates().shape[0]),
            "unique_peptides_matched_to_vcf": int(matched[peptide_key_cols].drop_duplicates().shape[0]),
            "unique_peptides_pass": int(pass_df[peptide_key_cols].drop_duplicates().shape[0]),
            "unique_peptides_nonpass": int(nonpass_df[peptide_key_cols].drop_duplicates().shape[0]),
        })
    else:
        row.update({
            "unique_peptides_total": pd.NA,
            "unique_peptides_matched_to_vcf": pd.NA,
            "unique_peptides_pass": pd.NA,
            "unique_peptides_nonpass": pd.NA,
        })

    return row


def process_one_pvac_file(pvac_tsv, vcf_dir):
    sample = infer_sample_from_pvac_path(pvac_tsv)
    mhc_class = infer_class_from_pvac_path(pvac_tsv)
    vcf_path = find_vcf_for_sample(vcf_dir, sample)

    pvac = pd.read_csv(pvac_tsv, sep="\t", dtype=str, low_memory=False).fillna("")
    original_cols = list(pvac.columns)

    pvac["sample"] = sample
    pvac["mhc_class"] = mhc_class
    pvac["pvac_tsv"] = str(pvac_tsv)
    pvac["source_vcf"] = str(vcf_path)

    required_cols = ["Chromosome", "Start", "Reference", "Variant"]
    missing = [c for c in required_cols if c not in pvac.columns]
    if missing:
        raise ValueError(f"{pvac_tsv} is missing required columns: {missing}")

    pvac["chrom_norm"] = pvac["Chromosome"].map(norm_chrom)

    vcf_filters = read_vcf_filters(vcf_path)

    merged, coordinate_mode, n_exact_matches = try_merge_with_coordinate_modes(
        pvac=pvac,
        vcf_filters=vcf_filters,
    )

    merged["coordinate_match_mode"] = coordinate_mode

    # Recover remaining exact CHROM/POS/REF/ALT matches using row-wise alternative
    # coordinate modes, mainly for indel/VEP representation edge cases.
    merged = fill_unmatched_with_rowwise_alternative_exact_matches(
        merged=merged,
        vcf_filters=vcf_filters,
    )

    merged["would_pass_pass_only"] = [bool(is_pass_filter(x)) for x in merged["VCF_FILTER"]]
    merged["would_pass_pass_only"] = merged["would_pass_pass_only"].astype(bool)

    merged["is_nonpass"] = (
        merged["VCF_FILTER"].notna().astype(bool)
        & (~merged["would_pass_pass_only"]).astype(bool)
    )
    merged["is_nonpass"] = merged["is_nonpass"].astype(bool)

    summary = summarize_one(
        merged=merged,
        sample=sample,
        mhc_class=mhc_class,
        pvac_tsv=pvac_tsv,
        vcf_path=vcf_path,
        coordinate_mode=coordinate_mode,
    )

    filter_counts = (
        merged["VCF_FILTER"]
        .fillna("UNMATCHED_TO_VCF")
        .astype(str)
        .value_counts(dropna=False)
        .rename_axis("VCF_FILTER")
        .reset_index(name="pvacseq_rows")
    )
    filter_counts.insert(0, "sample", sample)
    filter_counts.insert(1, "mhc_class", mhc_class)

    pass_only_original_schema = merged.loc[merged["would_pass_pass_only"], original_cols].copy()

    return merged, pass_only_original_schema, summary, filter_counts


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--vcf-dir", required=True, type=Path)
    parser.add_argument("--pvac-dir", required=True, type=Path)
    parser.add_argument("--passonly-pvac-dir", required=True, type=Path)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--mhc-classes", nargs="+", default=DEFAULT_CLASSES)
    parser.add_argument("--report-name", default=".filtered.tsv")

    args = parser.parse_args()

    pvac_files = find_pvac_files(
        pvac_dir=args.pvac_dir,
        mhc_classes=args.mhc_classes,
        report_name=args.report_name,
    )

    if not pvac_files:
        raise SystemExit(f"No pVACseq filtered.tsv files found under {args.pvac_dir}")

    args.passonly_pvac_dir.mkdir(parents=True, exist_ok=True)

    all_annotated = []
    all_summaries = []
    all_filter_counts = []
    passonly_rows = []

    for pvac_tsv in pvac_files:
        print(f"[INFO] Processing {pvac_tsv}")

        annotated, pass_only_df, summary, filter_counts = process_one_pvac_file(
            pvac_tsv=pvac_tsv,
            vcf_dir=args.vcf_dir,
        )

        all_annotated.append(annotated)
        all_summaries.append(summary)
        all_filter_counts.append(filter_counts)

        rel = pvac_tsv.relative_to(args.pvac_dir)
        out_tsv = args.passonly_pvac_dir / rel
        out_tsv.parent.mkdir(parents=True, exist_ok=True)

        pass_only_df.to_csv(out_tsv, sep="\t", index=False)

        passonly_rows.append({
            "source_pvac_tsv": str(pvac_tsv),
            "passonly_pvac_tsv": str(out_tsv),
            "source_rows": int(summary["pvacseq_rows_total"]),
            "passonly_rows": int(pass_only_df.shape[0]),
            "removed_nonpass_or_unmatched_rows": int(summary["pvacseq_rows_total"] - pass_only_df.shape[0]),
        })

        matched = summary["pvacseq_rows_matched_to_vcf"]
        total = summary["pvacseq_rows_total"]
        if total and matched / total < 0.95:
            print(f"[WARN] Low VCF match rate for {summary['sample']} {summary['mhc_class']}: {matched}/{total}")

    annotated_df = pd.concat(all_annotated, ignore_index=True)
    summary_df = pd.DataFrame(all_summaries)
    filter_counts_df = pd.concat(all_filter_counts, ignore_index=True)
    passonly_df = pd.DataFrame(passonly_rows)

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    annotated_df.to_csv(f"{out_prefix}.annotated_pvacseq_rows.tsv.gz", sep="\t", index=False)
    summary_df.to_csv(f"{out_prefix}.per_sample_class_summary.tsv", sep="\t", index=False)
    filter_counts_df.to_csv(f"{out_prefix}.per_sample_class_filter_counts.tsv", sep="\t", index=False)
    passonly_df.to_csv(f"{out_prefix}.passonly_file_build_summary.tsv", sep="\t", index=False)

    cohort_by_class = (
        summary_df.groupby("mhc_class", as_index=False)[
            [
                "pvacseq_rows_total",
                "pvacseq_rows_matched_to_vcf",
                "pvacseq_rows_unmatched_to_vcf",
                "pvacseq_rows_pass",
                "pvacseq_rows_nonpass",
                "unique_variant_keys_total_in_pvacseq",
                "unique_variant_keys_matched_to_vcf",
                "unique_variant_keys_pass",
                "unique_variant_keys_nonpass",
                "unique_peptides_total",
                "unique_peptides_matched_to_vcf",
                "unique_peptides_pass",
                "unique_peptides_nonpass",
            ]
        ]
        .sum(numeric_only=True)
    )

    cohort_by_class["pvacseq_rows_pass_pct_of_matched"] = (
        100 * cohort_by_class["pvacseq_rows_pass"] / cohort_by_class["pvacseq_rows_matched_to_vcf"]
    )
    cohort_by_class["pvacseq_rows_nonpass_pct_of_matched"] = (
        100 * cohort_by_class["pvacseq_rows_nonpass"] / cohort_by_class["pvacseq_rows_matched_to_vcf"]
    )

    cohort_by_class.to_csv(f"{out_prefix}.cohort_by_class_summary.tsv", sep="\t", index=False)

    overall = pd.DataFrame([{
        "n_samples_with_any_report": int(summary_df["sample"].nunique()),
        "n_sample_class_reports": int(summary_df.shape[0]),
        "pvacseq_rows_total": int(summary_df["pvacseq_rows_total"].sum()),
        "pvacseq_rows_matched_to_vcf": int(summary_df["pvacseq_rows_matched_to_vcf"].sum()),
        "pvacseq_rows_unmatched_to_vcf": int(summary_df["pvacseq_rows_unmatched_to_vcf"].sum()),
        "pvacseq_rows_pass": int(summary_df["pvacseq_rows_pass"].sum()),
        "pvacseq_rows_nonpass": int(summary_df["pvacseq_rows_nonpass"].sum()),
        "pvacseq_rows_pass_pct_of_matched": 100 * summary_df["pvacseq_rows_pass"].sum() / summary_df["pvacseq_rows_matched_to_vcf"].sum(),
        "pvacseq_rows_nonpass_pct_of_matched": 100 * summary_df["pvacseq_rows_nonpass"].sum() / summary_df["pvacseq_rows_matched_to_vcf"].sum(),
        "unique_peptides_total_sum_by_sample_class": pd.to_numeric(summary_df["unique_peptides_total"], errors="coerce").sum(),
        "unique_peptides_pass_sum_by_sample_class": pd.to_numeric(summary_df["unique_peptides_pass"], errors="coerce").sum(),
        "unique_peptides_nonpass_sum_by_sample_class": pd.to_numeric(summary_df["unique_peptides_nonpass"], errors="coerce").sum(),
    }])

    overall.to_csv(f"{out_prefix}.overall_summary.tsv", sep="\t", index=False)

    annotated_df[
        annotated_df["VCF_FILTER"].isna()
        | (annotated_df["VCF_FILTER"].astype(str).str.strip() == "")
    ].to_csv(f"{out_prefix}.unmatched_to_vcf.tsv", sep="\t", index=False)

    annotated_df[
        annotated_df["VCF_FILTER"].notna()
        & (~annotated_df["would_pass_pass_only"])
    ].to_csv(f"{out_prefix}.nonpass_rows.tsv", sep="\t", index=False)

    print("[DONE]")
    print(f"Annotated rows: {out_prefix}.annotated_pvacseq_rows.tsv.gz")
    print(f"Overall summary: {out_prefix}.overall_summary.tsv")
    print(f"PASS-only pVACseq dir: {args.passonly_pvac_dir}")


if __name__ == "__main__":
    main()
