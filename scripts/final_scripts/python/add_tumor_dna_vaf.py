#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys

import pandas as pd


def read_table(path: str, sep: str):
    if path in {"", "/dev/null"}:
        return None
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return None
    return pd.read_csv(p, sep=sep)


def norm_value(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def norm_class(x):
    s = norm_value(x).lower().replace("_", " ").replace("-", " ")
    if s in {"i", "1", "mhc class i", "class i", "mhci", "mhc i"}:
        return "I"
    if s in {"ii", "2", "mhc class ii", "class ii", "mhcii", "mhc ii"}:
        return "II"
    if "class i" in s and "ii" not in s:
        return "I"
    if "class ii" in s or "mhcii" in s:
        return "II"
    return norm_value(x)


def collapse_unique(series):
    vals = []
    seen = set()
    for x in series:
        if pd.isna(x):
            continue
        s = str(x).strip()
        if s == "" or s.lower() == "nan":
            continue
        if s not in seen:
            vals.append(s)
            seen.add(s)
    if not vals:
        return pd.NA
    return ";".join(vals)


def prepare_pvac(path: str, hla_class: str):
    df = read_table(path, sep="\t")
    if df is None:
        return pd.DataFrame()

    required = [
        "Gene Name",
        "HLA Allele",
        "MT Epitope Seq",
        "WT Epitope Seq",
        "Tumor DNA Depth",
        "Tumor DNA VAF",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")

    out = df[required].copy()
    out = out.rename(columns={
        "Gene Name": "Gene.Name",
        "HLA Allele": "HLA.Allele",
        "MT Epitope Seq": "MT.Epitope.Seq",
        "WT Epitope Seq": "WT.Epitope.Seq",
        "Tumor DNA Depth": "Tumor.DNA.Depth",
        "Tumor DNA VAF": "Tumor.DNA.VAF",
    })
    out["HLA.Class.norm"] = hla_class

    key_cols = [
        "Gene.Name",
        "HLA.Allele",
        "MT.Epitope.Seq",
        "WT.Epitope.Seq",
        "HLA.Class.norm",
    ]

    for c in key_cols:
        out[c] = out[c].map(norm_value)

    # Collapse duplicate pVACseq rows with the same final-table key.
    out = (
        out.groupby(key_cols, dropna=False, as_index=False)
        .agg({
            "Tumor.DNA.Depth": collapse_unique,
            "Tumor.DNA.VAF": collapse_unique,
        })
    )

    return out


def main():
    ap = argparse.ArgumentParser(
        description="Add Tumor DNA VAF from somatic pVACseq MHC-I/II outputs to final EGG epitope tables."
    )
    ap.add_argument("--sample", required=True)
    ap.add_argument("--final-in", required=True)
    ap.add_argument("--somatic-mhci", required=True)
    ap.add_argument("--somatic-mhcii", required=True)
    ap.add_argument("--final-out", required=True)
    args = ap.parse_args()

    final = read_table(args.final_in, sep=",")
    if final is None:
        raise FileNotFoundError(f"Could not read final table: {args.final_in}")

    required_final = [
        "Gene.Name",
        "HLA.Allele",
        "MT.Epitope.Seq",
        "WT.Epitope.Seq",
        "HLA.Class",
    ]
    missing_final = [c for c in required_final if c not in final.columns]
    if missing_final:
        raise ValueError(f"Final table is missing required columns: {missing_final}")

    pvac_i = prepare_pvac(args.somatic_mhci, "I")
    pvac_ii = prepare_pvac(args.somatic_mhcii, "II")
    pvac = pd.concat([pvac_i, pvac_ii], ignore_index=True)

    # Start with explicit NA columns so non-somatic/fusion/splicing rows remain present.
    final["Tumor.DNA.Depth"] = pd.NA
    final["Tumor.DNA.VAF"] = pd.NA

    if not pvac.empty:
        final["_row_id"] = range(len(final))
        final["HLA.Class.norm"] = final["HLA.Class"].map(norm_class)

        key_cols = [
            "Gene.Name",
            "HLA.Allele",
            "MT.Epitope.Seq",
            "WT.Epitope.Seq",
            "HLA.Class.norm",
        ]
        for c in key_cols:
            final[c] = final[c].map(norm_value)

        merged = final.merge(
            pvac,
            on=key_cols,
            how="left",
            suffixes=("", ".pvac"),
        )

        # Use matched pVAC values where available.
        merged["Tumor.DNA.Depth"] = merged["Tumor.DNA.Depth.pvac"]
        merged["Tumor.DNA.VAF"] = merged["Tumor.DNA.VAF.pvac"]

        merged = merged.drop(columns=[
            c for c in ["Tumor.DNA.Depth.pvac", "Tumor.DNA.VAF.pvac", "HLA.Class.norm"]
            if c in merged.columns
        ])

        # Restore original row order and remove helper.
        merged = merged.sort_values("_row_id").drop(columns=["_row_id"])
        final = merged

    out = Path(args.final_out)
    out.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(out, index=False)

    n_total = len(final)
    n_vaf = final["Tumor.DNA.VAF"].notna().sum()
    n_somatic = (
        final["Mutation.Source"].astype(str).str.lower().str.contains("somatic", na=False).sum()
        if "Mutation.Source" in final.columns else "NA"
    )

    print(f"[{args.sample}] wrote: {out}")
    print(f"[{args.sample}] rows: {n_total}")
    print(f"[{args.sample}] rows with Tumor.DNA.VAF: {n_vaf}")
    print(f"[{args.sample}] rows labelled somatic: {n_somatic}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise
