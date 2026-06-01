#!/usr/bin/env python3

import gzip
from pathlib import Path
import pandas as pd


BASE = Path("TCGA_melanoma")
RR4 = BASE / "rescue_final_analysis/reviewer_updates/reviewer_4_nonpass_passonly"

ANNOTATED = RR4 / "complete_rescue_pvacseq_filter_status.annotated_pvacseq_rows.tsv.gz"
VCF_DIR = BASE / "2A_somatic_mutation_epitopes_complete_rescue/kallisto_somatic_VCF_full"

OUT_PREFIX = RR4 / "complete_rescue_pvacseq_unmatched_diagnosis"


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def norm_chrom(chrom):
    chrom = str(chrom)
    return chrom[3:] if chrom.startswith("chr") else chrom


def find_vcf_for_sample(sample):
    candidates = [
        VCF_DIR / f"{sample}_somatic_quantification.full_rescue.vcf.gz",
        VCF_DIR / f"{sample}_somatic_quantification.full_rescue.vcf",
        VCF_DIR / f"{sample}_somatic_quantification.vcf.gz",
        VCF_DIR / f"{sample}_somatic_quantification.vcf",
    ]

    for c in candidates:
        if c.exists():
            return c

    raise FileNotFoundError(f"No VCF found for {sample}")


def read_vcf_records(vcf_path):
    rows = []

    with open_text(vcf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom, pos, vid, ref, alts, qual, filt, info = fields[:8]
            pos = int(pos)

            for alt in alts.split(","):
                rows.append({
                    "chrom_norm": norm_chrom(chrom),
                    "vcf_pos": pos,
                    "Reference": ref,
                    "Variant": alt,
                    "VCF_FILTER": filt,
                    "vcf_key": f"{norm_chrom(chrom)}:{pos}:{ref}>{alt}",
                })

    return pd.DataFrame(rows)


def as_int(x):
    try:
        return int(float(x))
    except Exception:
        return None


def diagnose_row(row, vcf):
    chrom = norm_chrom(row.get("Chromosome", ""))
    start = as_int(row.get("Start", ""))
    stop = as_int(row.get("Stop", ""))
    ref = str(row.get("Reference", ""))
    alt = str(row.get("Variant", ""))

    tests = []

    positions = []
    if start is not None:
        positions.extend([
            ("Start_minus_2", start - 2),
            ("Start_minus_1", start - 1),
            ("Start", start),
            ("Start_plus_1", start + 1),
            ("Start_plus_2", start + 2),
        ])

    if stop is not None:
        positions.extend([
            ("Stop_minus_1", stop - 1),
            ("Stop", stop),
            ("Stop_plus_1", stop + 1),
        ])

    seen = set()
    unique_positions = []

    for mode, pos in positions:
        key = (mode, pos)
        if key not in seen:
            seen.add(key)
            unique_positions.append((mode, pos))

    # Exact CHROM/POS/REF/ALT tests
    for mode, pos in unique_positions:
        exact = vcf[
            (vcf["chrom_norm"] == chrom)
            & (vcf["vcf_pos"] == pos)
            & (vcf["Reference"] == ref)
            & (vcf["Variant"] == alt)
        ]

        if len(exact):
            return {
                "diagnosis": "EXACT_MATCH_WITH_ALTERNATIVE_POSITION",
                "matched_mode": mode,
                "matched_pos": pos,
                "matched_filter_values": ";".join(sorted(set(exact["VCF_FILTER"].astype(str)))),
                "matched_vcf_keys": ";".join(sorted(set(exact["vcf_key"].astype(str)))),
            }

    # Same CHROM/POS but different REF/ALT
    chrom_pos_hits = []

    for mode, pos in unique_positions:
        hits = vcf[
            (vcf["chrom_norm"] == chrom)
            & (vcf["vcf_pos"] == pos)
        ]

        if len(hits):
            chrom_pos_hits.append((mode, pos, hits))

    if chrom_pos_hits:
        mode, pos, hits = chrom_pos_hits[0]
        return {
            "diagnosis": "POSITION_MATCH_BUT_REF_ALT_DIFFER",
            "matched_mode": mode,
            "matched_pos": pos,
            "matched_filter_values": ";".join(sorted(set(hits["VCF_FILTER"].astype(str)))),
            "matched_vcf_keys": ";".join(sorted(set(hits["vcf_key"].astype(str)))[:10]),
        }

    # Same CHROM/REF/ALT anywhere nearby ±20 bp
    if start is not None:
        nearby = vcf[
            (vcf["chrom_norm"] == chrom)
            & (vcf["Reference"] == ref)
            & (vcf["Variant"] == alt)
            & (vcf["vcf_pos"] >= start - 20)
            & (vcf["vcf_pos"] <= start + 20)
        ]

        if len(nearby):
            return {
                "diagnosis": "SAME_REF_ALT_NEARBY_WITHIN_20BP",
                "matched_mode": "nearby_20bp",
                "matched_pos": ";".join(map(str, sorted(set(nearby["vcf_pos"])))),
                "matched_filter_values": ";".join(sorted(set(nearby["VCF_FILTER"].astype(str)))),
                "matched_vcf_keys": ";".join(sorted(set(nearby["vcf_key"].astype(str)))[:10]),
            }

    # Same chromosome only near position
    if start is not None:
        nearby_any = vcf[
            (vcf["chrom_norm"] == chrom)
            & (vcf["vcf_pos"] >= start - 20)
            & (vcf["vcf_pos"] <= start + 20)
        ]

        if len(nearby_any):
            return {
                "diagnosis": "NEARBY_VARIANT_WITHIN_20BP_BUT_NOT_SAME_REF_ALT",
                "matched_mode": "nearby_20bp",
                "matched_pos": ";".join(map(str, sorted(set(nearby_any["vcf_pos"])))),
                "matched_filter_values": ";".join(sorted(set(nearby_any["VCF_FILTER"].astype(str)))),
                "matched_vcf_keys": ";".join(sorted(set(nearby_any["vcf_key"].astype(str)))[:10]),
            }

    return {
        "diagnosis": "NO_MATCH_FOUND_IN_TESTED_MODES",
        "matched_mode": "",
        "matched_pos": "",
        "matched_filter_values": "",
        "matched_vcf_keys": "",
    }


def main():
    df = pd.read_csv(ANNOTATED, sep="\t", dtype=str, low_memory=False).fillna("")

    unmatched = df[
        (df["VCF_FILTER"].astype(str).str.strip() == "")
        | (df["VCF_FILTER"].isna())
    ].copy()

    print(f"Unmatched rows: {len(unmatched)}")

    out_rows = []
    vcf_cache = {}

    for idx, row in unmatched.iterrows():
        sample = row["sample"]

        if sample not in vcf_cache:
            vcf_path = find_vcf_for_sample(sample)
            vcf_cache[sample] = (vcf_path, read_vcf_records(vcf_path))

        vcf_path, vcf = vcf_cache[sample]

        diag = diagnose_row(row, vcf)

        out = row.to_dict()
        out["diagnosis"] = diag["diagnosis"]
        out["matched_mode"] = diag["matched_mode"]
        out["matched_pos"] = diag["matched_pos"]
        out["matched_filter_values"] = diag["matched_filter_values"]
        out["matched_vcf_keys"] = diag["matched_vcf_keys"]
        out["diagnostic_vcf"] = str(vcf_path)

        out_rows.append(out)

    out_df = pd.DataFrame(out_rows)

    out_tsv = Path(str(OUT_PREFIX) + ".tsv")
    out_summary = Path(str(OUT_PREFIX) + ".summary.tsv")

    out_df.to_csv(out_tsv, sep="\t", index=False)

    summary = (
        out_df
        .groupby("diagnosis", as_index=False)
        .size()
        .rename(columns={"size": "n_rows"})
        .sort_values("n_rows", ascending=False)
    )

    summary.to_csv(out_summary, sep="\t", index=False)

    print("Diagnosis summary:")
    print(summary.to_string(index=False))
    print()
    print("Wrote:")
    print(out_tsv)
    print(out_summary)


if __name__ == "__main__":
    main()
