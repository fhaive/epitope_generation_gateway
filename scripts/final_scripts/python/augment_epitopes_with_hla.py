#!/usr/bin/env python3
"""
augment_epitopes_with_hla.py

Add HLA_TPM (and a detail column) to a final epitope CSV using arcasHLA outputs:
  <SAMPLE>_CancerRNA_LRNA_sorted.genotype.json
  <SAMPLE>_CancerRNA_LRNA_sorted.genes.json

Interpretation of genes.json values: [TPM, READ_COUNT, REL_ABUNDANCE]
We use TPM for HLA_TPM. For class-II paired alleles (e.g., DQA1-DQB1, DPA1-DPB1),
HLA_TPM = min(TPM_chain1, TPM_chain2). The detail column preserves both chains'
TPM, read count, and relative abundance for reference.

The column HLA_TPM is inserted immediately after HLA.Allele.
"""
import argparse
import json
import os
import re
import sys
import pandas as pd

def gene_keys_from_allele(allele: str):
    if not allele:
        return []
    a = str(allele).strip()
    keys = []
    for part in a.split("-"):
        m = re.match(r"^(?:HLA-)?([A-Z0-9]+)\*", part)
        if m:
            keys.append(m.group(1))
    return keys

def load_arcas_genes(genes_json_path: str):
    with open(genes_json_path, "r") as f:
        d = json.load(f)
    out = {}
    for gene, vals in d.items():
        tpm = cnt = rel = None
        if isinstance(vals, (list, tuple)):
            if len(vals) >= 3:
                tpm, cnt, rel = vals[0], vals[1], vals[2]
            elif len(vals) == 2:
                cnt, rel = vals
        out[gene] = (tpm, cnt, rel)
    return out

def insert_after(df: pd.DataFrame, newcol: str, series, after_col: str) -> pd.DataFrame:
    df = df.copy()
    if after_col in df.columns:
        idx = list(df.columns).index(after_col) + 1
        left = df.iloc[:, :idx]
        right = df.iloc[:, idx:]
        return pd.concat([left, pd.Series(series, name=newcol), right], axis=1)
    else:
        df[newcol] = series
        return df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--final-in", required=True)
    ap.add_argument("--arcas-dir", required=True)
    ap.add_argument("--final-out", required=True)
    args = ap.parse_args()

    genes_json = os.path.join(args.arcas_dir, f"{args.sample}_CancerRNA_LRNA_sorted.genes.json")
    if not os.path.isfile(genes_json):
        sys.exit(f"genes.json not found: {genes_json}")

    df = pd.read_csv(args.final_in, dtype=str).fillna("")
    genes = load_arcas_genes(genes_json)

    hla_col = "HLA.Allele" if "HLA.Allele" in df.columns else ("HLA_Allele" if "HLA_Allele" in df.columns else None)
    if not hla_col:
        sys.exit("Could not find HLA allele column in the table.")

    tpm_vals = []
    details = []
    for _, r in df.iterrows():
        allele = r.get(hla_col, "")
        keys = gene_keys_from_allele(allele)
        per_chain = []
        for k in keys:
            if k in genes:
                tpm, cnt, rel = genes[k]
                per_chain.append((k, tpm, cnt, rel))
        if not per_chain:
            tpm_vals.append("")
            details.append("")
            continue
        if len(per_chain) == 1:
            k, tpm, cnt, rel = per_chain[0]
            tpm_vals.append("" if tpm is None else float(tpm))
            details.append(f"{k}: TPM={tpm}; count={cnt}; rel={rel}")
        else:
            chain_tpms = [p[1] for p in per_chain if p[1] is not None]
            limiting = min(chain_tpms) if chain_tpms else ""
            tpm_vals.append(limiting if limiting != "" else "")
            parts = [f"{k}: TPM={tpm}; count={cnt}; rel={rel}" for (k,tpm,cnt,rel) in per_chain]
            details.append("; ".join(parts))

    df = insert_after(df, "HLA_TPM", tpm_vals, after_col=hla_col)
    df["HLA_TPM_Detail"] = details

    os.makedirs(os.path.dirname(args.final_out), exist_ok=True)
    df.to_csv(args.final_out, index=False)

if __name__ == "__main__":
    main()
