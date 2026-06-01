#!/usr/bin/env python3

from pathlib import Path
import re

import pandas as pd


BASE = Path("TCGA_melanoma")

ANNOTATED_PVAC = BASE / "rescue_final_analysis/reviewer_updates/reviewer_4_nonpass_passonly/complete_rescue_pvacseq_filter_status.annotated_pvacseq_rows.tsv.gz"

FULL_FINAL_DIR = BASE / "epitopes_prioritisation_complete_rescue/final_epitopes"
PASSONLY_FINAL_DIR = BASE / "epitopes_prioritisation_complete_rescue_PASS_ONLY/final_epitopes"

OUT_DIR = BASE / "rescue_final_analysis/reviewer_updates/reviewer_4_nonpass_passonly"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TOP_NS = [10, 20, 50]


def clean_str(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def norm_sample(x):
    x = clean_str(x)
    if x.endswith("_epitopes_final.csv"):
        x = x.replace("_epitopes_final.csv", "")
    if x.startswith("Sample_"):
        return x
    m = re.search(r"(TCGA-[A-Z0-9]+-[A-Z0-9]+)", x)
    if m:
        return "Sample_" + m.group(1)
    return x


def norm_mhc(x):
    x = clean_str(x)
    if x in {"MHC_Class_I", "MHC I", "MHC Class I", "I"}:
        return "MHC_Class_I"
    if x in {"MHC_Class_II", "MHC II", "MHC Class II", "II"}:
        return "MHC_Class_II"
    return x


def norm_hla(x):
    x = clean_str(x)
    x = x.replace("HLA-", "HLA-")
    return x


def get_first_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


def build_pvac_status_lookup():
    pvac = pd.read_csv(ANNOTATED_PVAC, sep="\t", dtype=str, low_memory=False).fillna("")

    gene_col = get_first_col(pvac, ["Gene Name", "Gene", "Gene.Name", "Gene.Symbol"])
    peptide_col = get_first_col(pvac, ["MT Epitope Seq", "MT.Epitope.Seq", "Epitope Seq", "Epitope.Seq", "Peptide"])
    hla_col = get_first_col(pvac, ["HLA Allele", "HLA.Allele", "HLA"])

    required = {
        "gene_col": gene_col,
        "peptide_col": peptide_col,
        "hla_col": hla_col,
    }

    missing = [k for k, v in required.items() if v is None]
    if missing:
        raise ValueError(f"Could not identify required pVAC columns: {missing}. Columns: {pvac.columns.tolist()}")

    pvac["key_sample"] = pvac["sample"].map(norm_sample)
    pvac["key_mhc"] = pvac["mhc_class"].map(norm_mhc)
    pvac["key_gene"] = pvac[gene_col].map(clean_str)
    pvac["key_peptide"] = pvac[peptide_col].map(clean_str)
    pvac["key_hla"] = pvac[hla_col].map(norm_hla)

    pvac["filter_class"] = "UNMATCHED_TO_VCF"
    pvac.loc[pvac["VCF_FILTER"].astype(str).str.len() > 0, "filter_class"] = "NONPASS"
    pvac.loc[pvac["would_pass_pass_only"].astype(str).str.lower().isin(["true", "1"]), "filter_class"] = "PASS"

    key_cols = ["key_sample", "key_mhc", "key_gene", "key_peptide", "key_hla"]

    rows = []

    for key, group in pvac.groupby(key_cols, dropna=False):
        statuses = set(group["filter_class"].astype(str))

        if statuses == {"PASS"}:
            status = "PASS"
        elif statuses == {"NONPASS"}:
            status = "NONPASS"
        elif statuses == {"UNMATCHED_TO_VCF"}:
            status = "UNMATCHED_TO_VCF"
        elif "PASS" in statuses and "NONPASS" in statuses:
            status = "MIXED_PASS_NONPASS"
        elif "PASS" in statuses:
            status = "PASS_WITH_UNMATCHED"
        elif "NONPASS" in statuses:
            status = "NONPASS_WITH_UNMATCHED"
        else:
            status = ";".join(sorted(statuses))

        row = dict(zip(key_cols, key))
        row.update({
            "somatic_filter_status": status,
            "n_pvac_rows_for_key": int(group.shape[0]),
            "VCF_FILTER_values": ";".join(sorted(set(group["VCF_FILTER"].astype(str)))),
            "coordinate_match_modes": ";".join(sorted(set(group["coordinate_match_mode"].astype(str)))),
        })
        rows.append(row)

    lookup = pd.DataFrame(rows)
    lookup.to_csv(OUT_DIR / "somatic_pvacseq_filter_status_lookup.tsv", sep="\t", index=False)

    return lookup


def read_final_tables(final_dir):
    rows = []

    for f in sorted(final_dir.glob("*_epitopes_final.csv")):
        sample = norm_sample(f.name)
        df = pd.read_csv(f, dtype=str).fillna("")
        df["sample"] = sample
        df["final_table"] = str(f)
        rows.append(df)

    if not rows:
        raise FileNotFoundError(f"No final tables found in {final_dir}")

    return pd.concat(rows, ignore_index=True)


def driver_flag_value(x):
    x = clean_str(x)
    if x == "":
        return False
    if x.lower() in {"nan", "na", "none", ".", "false"}:
        return False
    return True


def annotate_final_with_pvac_status(final_df, lookup):
    final = final_df.copy()

    source_col = get_first_col(final, ["Mutation.Source", "Mutation Source"])
    gene_col = get_first_col(final, ["Gene.Name", "Gene Name", "Gene"])
    peptide_col = get_first_col(final, ["MT.Epitope.Seq", "MT Epitope Seq", "Epitope Seq", "Peptide"])
    hla_col = get_first_col(final, ["HLA.Allele", "HLA Allele", "HLA"])
    mhc_col = get_first_col(final, ["HLA.Class", "MHC_Class", "MHC Class"])

    required = {
        "source_col": source_col,
        "gene_col": gene_col,
        "peptide_col": peptide_col,
        "hla_col": hla_col,
        "mhc_col": mhc_col,
    }

    missing = [k for k, v in required.items() if v is None]
    if missing:
        raise ValueError(f"Could not identify required final-table columns: {missing}. Columns: {final.columns.tolist()}")

    final["key_sample"] = final["sample"].map(norm_sample)
    final["key_mhc"] = final[mhc_col].map(norm_mhc)
    final["key_gene"] = final[gene_col].map(clean_str)
    final["key_peptide"] = final[peptide_col].map(clean_str)
    final["key_hla"] = final[hla_col].map(norm_hla)

    key_cols = ["key_sample", "key_mhc", "key_gene", "key_peptide", "key_hla"]

    annotated = final.merge(
        lookup,
        how="left",
        on=key_cols,
    )

    annotated["somatic_filter_status"] = annotated["somatic_filter_status"].fillna("NOT_SOMATIC_OR_NO_PVAC_MATCH")

    is_somatic = annotated[source_col].astype(str).str.lower().eq("somatic")
    annotated.loc[~is_somatic, "somatic_filter_status"] = "NOT_SOMATIC"

    annotated["is_somatic"] = is_somatic
    annotated["is_somatic_pass"] = is_somatic & annotated["somatic_filter_status"].isin(["PASS", "PASS_WITH_UNMATCHED"])
    annotated["is_somatic_nonpass"] = is_somatic & annotated["somatic_filter_status"].isin(["NONPASS", "NONPASS_WITH_UNMATCHED"])
    annotated["is_somatic_mixed"] = is_somatic & annotated["somatic_filter_status"].eq("MIXED_PASS_NONPASS")
    annotated["is_somatic_unmatched"] = is_somatic & annotated["somatic_filter_status"].isin(["UNMATCHED_TO_VCF", "NOT_SOMATIC_OR_NO_PVAC_MATCH"])

    return annotated


def topn_nonpass_summary(annotated):
    rank_col = get_first_col(annotated, ["Borda_Rank", "Borda Rank"])
    source_col = get_first_col(annotated, ["Mutation.Source", "Mutation Source"])

    annotated = annotated.copy()
    annotated["rank_numeric"] = pd.to_numeric(annotated[rank_col], errors="coerce")

    rows = []

    for sample, group in annotated.groupby("sample", dropna=False):
        group = group.sort_values("rank_numeric", na_position="last").copy()

        for n in TOP_NS + ["all"]:
            if n == "all":
                sub = group.copy()
                n_label = "all"
            else:
                sub = group.head(n).copy()
                n_label = str(n)

            total = int(sub.shape[0])
            somatic = int(sub["is_somatic"].sum())
            nonpass = int(sub["is_somatic_nonpass"].sum())
            pass_n = int(sub["is_somatic_pass"].sum())
            mixed = int(sub["is_somatic_mixed"].sum())
            unmatched = int(sub["is_somatic_unmatched"].sum())

            fusion = int(sub[source_col].astype(str).str.lower().eq("fusion").sum())
            splice = int(sub[source_col].astype(str).str.lower().isin(["splice", "splicing"]).sum())

            rows.append({
                "sample": sample,
                "top_n": n_label,
                "total_rows": total,
                "somatic_rows": somatic,
                "fusion_rows": fusion,
                "splice_rows": splice,
                "somatic_PASS_rows": pass_n,
                "somatic_nonPASS_rows": nonpass,
                "somatic_mixed_PASS_nonPASS_rows": mixed,
                "somatic_unmatched_rows": unmatched,
                "nonPASS_fraction_of_all_top_rows": nonpass / total if total else pd.NA,
                "nonPASS_fraction_of_somatic_top_rows": nonpass / somatic if somatic else pd.NA,
            })

    return pd.DataFrame(rows)


def driver_enrichment_summary(final_df, label):
    rank_col = get_first_col(final_df, ["Borda_Rank", "Borda Rank"])
    driver_col = get_first_col(final_df, ["Intogen_Driver_role", "Intogen Driver role", "Intogen.Driver.role"])
    gene_col = get_first_col(final_df, ["Gene.Name", "Gene Name", "Gene"])

    if driver_col is None:
        raise ValueError(f"No IntOGen driver column found. Columns: {final_df.columns.tolist()}")

    df = final_df.copy()
    df["rank_numeric"] = pd.to_numeric(df[rank_col], errors="coerce")
    df["is_driver"] = df[driver_col].map(driver_flag_value)

    rows = []

    for sample, group in df.groupby("sample", dropna=False):
        group = group.sort_values("rank_numeric", na_position="last").copy()
        top10 = group.head(10).copy()
        rest = group.iloc[10:].copy()

        top10_driver_genes = sorted(set(top10.loc[top10["is_driver"], gene_col].astype(str))) if gene_col else []

        rows.append({
            "analysis": label,
            "sample": sample,
            "total_rows": int(group.shape[0]),
            "total_driver_rows": int(group["is_driver"].sum()),
            "top10_rows": int(top10.shape[0]),
            "top10_driver_rows": int(top10["is_driver"].sum()),
            "top10_driver_fraction": top10["is_driver"].mean() if len(top10) else pd.NA,
            "rest_rows": int(rest.shape[0]),
            "rest_driver_rows": int(rest["is_driver"].sum()) if len(rest) else 0,
            "rest_driver_fraction": rest["is_driver"].mean() if len(rest) else pd.NA,
            "top10_unique_driver_genes": ";".join(top10_driver_genes),
        })

    return pd.DataFrame(rows)


def cohort_driver_summary(per_sample):
    rows = []

    for label, group in per_sample.groupby("analysis", dropna=False):
        rows.append({
            "analysis": label,
            "n_samples": int(group["sample"].nunique()),
            "top10_rows_total": int(group["top10_rows"].sum()),
            "top10_driver_rows_total": int(group["top10_driver_rows"].sum()),
            "top10_driver_fraction_cohort": group["top10_driver_rows"].sum() / group["top10_rows"].sum(),
            "all_rows_total": int(group["total_rows"].sum()),
            "all_driver_rows_total": int(group["total_driver_rows"].sum()),
            "all_driver_fraction_cohort": group["total_driver_rows"].sum() / group["total_rows"].sum(),
        })

    return pd.DataFrame(rows)


def main():
    lookup = build_pvac_status_lookup()

    full = read_final_tables(FULL_FINAL_DIR)
    full_annotated = annotate_final_with_pvac_status(full, lookup)

    full_annotated.to_csv(
        OUT_DIR / "complete_rescue_final_epitopes_with_somatic_filter_status.tsv.gz",
        sep="\t",
        index=False,
    )

    topn = topn_nonpass_summary(full_annotated)
    topn.to_csv(OUT_DIR / "complete_rescue_topN_nonPASS_contribution.tsv", sep="\t", index=False)

    rank_col = get_first_col(full_annotated, ["Borda_Rank", "Borda Rank"])
    top10_rows = (
        full_annotated.assign(rank_numeric=pd.to_numeric(full_annotated[rank_col], errors="coerce"))
        .sort_values(["sample", "rank_numeric"])
        .groupby("sample", as_index=False)
        .head(10)
    )
    top10_rows.to_csv(OUT_DIR / "complete_rescue_top10_epitopes_with_somatic_filter_status.tsv", sep="\t", index=False)

    driver_tables = []

    driver_tables.append(driver_enrichment_summary(full, "complete_rescue_all_somatic"))

    if PASSONLY_FINAL_DIR.exists():
        passonly = read_final_tables(PASSONLY_FINAL_DIR)
        driver_tables.append(driver_enrichment_summary(passonly, "complete_rescue_PASS_only_somatic"))

    driver_per_sample = pd.concat(driver_tables, ignore_index=True)
    driver_per_sample.to_csv(OUT_DIR / "driver_enrichment_top10_per_sample_full_vs_PASS_only.tsv", sep="\t", index=False)

    driver_cohort = cohort_driver_summary(driver_per_sample)
    driver_cohort.to_csv(OUT_DIR / "driver_enrichment_top10_cohort_summary_full_vs_PASS_only.tsv", sep="\t", index=False)

    # Compact text summary for reviewer response.
    overall_filter = pd.read_csv(
        OUT_DIR / "complete_rescue_pvacseq_filter_status.overall_summary.tsv",
        sep="\t",
    ) if (OUT_DIR / "complete_rescue_pvacseq_filter_status.overall_summary.tsv").exists() else None

    top10 = topn[topn["top_n"].astype(str) == "10"].copy()

    text_lines = []

    text_lines.append("Reviewer 4 PASS/non-PASS sensitivity summary")
    text_lines.append("")
    text_lines.append("Top-10 non-PASS contribution from final complete-rescue prioritisation:")
    text_lines.append(f"  Samples evaluated: {top10['sample'].nunique()}")
    text_lines.append(f"  Total top-10 rows: {int(top10['total_rows'].sum())}")
    text_lines.append(f"  Top-10 somatic rows: {int(top10['somatic_rows'].sum())}")
    text_lines.append(f"  Top-10 somatic PASS rows: {int(top10['somatic_PASS_rows'].sum())}")
    text_lines.append(f"  Top-10 somatic non-PASS rows: {int(top10['somatic_nonPASS_rows'].sum())}")
    text_lines.append(f"  Top-10 somatic mixed PASS/non-PASS rows: {int(top10['somatic_mixed_PASS_nonPASS_rows'].sum())}")
    text_lines.append(f"  Top-10 somatic unmatched rows: {int(top10['somatic_unmatched_rows'].sum())}")

    total_top10 = top10["total_rows"].sum()
    total_somatic_top10 = top10["somatic_rows"].sum()
    total_nonpass_top10 = top10["somatic_nonPASS_rows"].sum()

    text_lines.append(
        f"  Non-PASS fraction of all top-10 rows: {total_nonpass_top10 / total_top10:.4f}"
        if total_top10 else "  Non-PASS fraction of all top-10 rows: NA"
    )
    text_lines.append(
        f"  Non-PASS fraction of somatic top-10 rows: {total_nonpass_top10 / total_somatic_top10:.4f}"
        if total_somatic_top10 else "  Non-PASS fraction of somatic top-10 rows: NA"
    )

    text_lines.append("")
    text_lines.append("Driver enrichment / driver representation in top-10:")
    for _, row in driver_cohort.iterrows():
        text_lines.append(
            f"  {row['analysis']}: top10_driver_rows={int(row['top10_driver_rows_total'])}/"
            f"{int(row['top10_rows_total'])} "
            f"({row['top10_driver_fraction_cohort']:.4f}); "
            f"all_driver_rows={int(row['all_driver_rows_total'])}/"
            f"{int(row['all_rows_total'])} "
            f"({row['all_driver_fraction_cohort']:.4f})"
        )

    out_txt = OUT_DIR / "reviewer_4_key_numbers_for_response.txt"
    out_txt.write_text("\n".join(text_lines) + "\n")

    print("[DONE] Wrote reviewer response outputs to:")
    print(OUT_DIR)
    print()
    print(out_txt.read_text())


if __name__ == "__main__":
    main()
