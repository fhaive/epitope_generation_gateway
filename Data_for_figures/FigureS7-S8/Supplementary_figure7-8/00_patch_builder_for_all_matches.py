#!/usr/bin/env python3

from pathlib import Path
import re

ROOT = Path.cwd()
SRC = ROOT / "01_build_exact_relaxed_evidence_flags_complete_rescue_vaf.py"
WORK = ROOT / "reviewer2_coverage_sensitivity"
SCRIPTS = WORK / "scripts"
SCRIPTS.mkdir(parents=True, exist_ok=True)

DST = SCRIPTS / "01_build_exact_relaxed_evidence_flags_complete_rescue_vaf_ALL_MATCHES.py"

if not SRC.exists():
    raise FileNotFoundError(f"Could not find original script: {SRC}")

text = SRC.read_text()

# ---------------------------------------------------------------------
# 1. Redirect output files from the copied script into a new results folder
# ---------------------------------------------------------------------
lines = text.splitlines()
new_lines = []
inserted_coverage_out = False

for line in lines:
    if (not inserted_coverage_out) and re.match(r'^[A-Z0-9_]+_OUT\s*=\s*OUT\s*/', line):
        new_lines.append("# Reviewer round 2 coverage-sensitivity outputs")
        new_lines.append('COVERAGE_OUT = OUT / "reviewer2_coverage_sensitivity" / "patched_builder_outputs"')
        new_lines.append("COVERAGE_OUT.mkdir(parents=True, exist_ok=True)")
        inserted_coverage_out = True

    if re.match(r'^[A-Z0-9_]+_OUT\s*=\s*OUT\s*/', line):
        line = re.sub(r'=\s*OUT\s*/', '= COVERAGE_OUT /', line)

    new_lines.append(line)

    if line.startswith("MATCH_EXAMPLES_OUT"):
        new_lines.append('ALL_MATCHES_OUT = COVERAGE_OUT / "all_matches.exact_relaxed.with_coverage.tsv.gz"')

text = "\n".join(new_lines) + "\n"

# ---------------------------------------------------------------------
# 2. Add a full all_matches list, while keeping the original examples file
# ---------------------------------------------------------------------
needle = """    examples = []
    seen_examples = set()
    scan_rows = []"""

replacement = """    examples = []
    seen_examples = set()
    all_matches = []
    scan_rows = []"""

if needle not in text:
    raise RuntimeError("Could not find examples/seen_examples/scan_rows block to patch.")

text = text.replace(needle, replacement)

# ---------------------------------------------------------------------
# 3. Add full pairwise match recording inside mark_candidate()
# ---------------------------------------------------------------------
needle = """        for col in base_cols:
            ensure_flag(col)[cid] = 1

        ex_key = (cid, scope, database, evidence_strength, evidence_class)"""

replacement = """        for col in base_cols:
            ensure_flag(col)[cid] = 1

        # Full pairwise match record for reviewer coverage-sensitivity analysis.
        # This is NOT deduplicated and should therefore be used for percent-coverage thresholds.
        cpep = id_to_peptide[cid]
        ref_pep_clean = ref_pep or ""

        candidate_contains_reference = bool(ref_pep_clean and ref_pep_clean in cpep)
        reference_contains_candidate = bool(ref_pep_clean and cpep in ref_pep_clean)
        is_containment = candidate_contains_reference or reference_contains_candidate

        if is_containment:
            overlap_len = min(len(cpep), len(ref_pep_clean))
        else:
            overlap_len = 0

        longer_peptide_length = max(len(cpep), len(ref_pep_clean)) if ref_pep_clean else 0
        longer_peptide_coverage = (
            overlap_len / longer_peptide_length
            if longer_peptide_length > 0
            else np.nan
        )

        all_matches.append({
            "candidate_id": cand.loc[cid, "candidate_id"],
            "sample_id": cand.loc[cid, "sample_id"],
            "candidate_gene": cand.loc[cid, "Gene.Name"],
            "candidate_hla_original": cand.loc[cid, "HLA.Allele"],
            "candidate_hla_normalized_joined": cand.loc[cid, "candidate_hla_normalized_joined"],
            "candidate_peptide": cpep,
            "candidate_peptide_length": len(cpep),
            "match_scope": scope,
            "database": database,
            "evidence_group": row_meta.get("evidence_group", ""),
            "evidence_class": evidence_class,
            "evidence_strength": evidence_strength,
            "reference_peptide": ref_pep_clean,
            "reference_peptide_length": len(ref_pep_clean),
            "overlap_len": overlap_len,
            "longer_peptide_length": longer_peptide_length,
            "longer_peptide_coverage": longer_peptide_coverage,
            "length_difference": abs(len(cpep) - len(ref_pep_clean)) if ref_pep_clean else "",
            "candidate_contains_reference": int(candidate_contains_reference),
            "reference_contains_candidate": int(reference_contains_candidate),
            "reference_hla_original": row_meta.get("hla_original", ""),
            "reference_hla_normalized": ref_hla_norm,
            "hla_match": int(hla_match),
            "reference_gene": ref_gene,
            "gene_match": int(gene_match),
            "pmid": row_meta.get("pmid", ""),
            "reference_title": row_meta.get("reference_title", ""),
            "journal": row_meta.get("journal", ""),
            "year": row_meta.get("year", ""),
            "source_file": row_meta.get("source_file", ""),
        })

        ex_key = (cid, scope, database, evidence_strength, evidence_class)"""

if needle not in text:
    raise RuntimeError("Could not find mark_candidate flag-setting block to patch.")

text = text.replace(needle, replacement)

# ---------------------------------------------------------------------
# 4. Write the full all_matches table near the end of the script
# ---------------------------------------------------------------------
needle = """    print("Writing match examples...")
    pd.DataFrame(examples).to_csv(MATCH_EXAMPLES_OUT, sep="\\t", index=False, compression="gzip")"""

replacement = """    print("Writing match examples...")
    pd.DataFrame(examples).to_csv(MATCH_EXAMPLES_OUT, sep="\\t", index=False, compression="gzip")

    print("Writing all exact/relaxed matches with coverage...")
    pd.DataFrame(all_matches).to_csv(
        ALL_MATCHES_OUT,
        sep="\\t",
        index=False,
        compression="gzip"
    )"""

if needle not in text:
    raise RuntimeError("Could not find match_examples writing block to patch.")

text = text.replace(needle, replacement)

DST.write_text(text)

print(f"Patched script written to: {DST}")
print("Original script was not modified.")
