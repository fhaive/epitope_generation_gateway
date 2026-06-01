#!/usr/bin/env python3

from pathlib import Path
import gzip
import csv
import re
from collections import Counter

IN_DIR = Path("/data/fsluma/pipelines/Epitope_Generation_Gateway/Immunopeptidomics_Data/normalized_reference_tables")
OUT_DIR = Path("normalized_reference_tables_hla_patched").resolve()
OUT_DIR.mkdir(exist_ok=True)

SUMMARY_FILE = OUT_DIR / "hla_patch_summary.tsv"
DETAIL_FILE = OUT_DIR / "hla_patch_details.tsv"


def clean(x):
    if x is None:
        return ""
    x = str(x).strip()
    if x.upper() in {"", "NA", "NAN", "NONE", "NULL"}:
        return ""
    return x


def pad2(x):
    return str(int(x)).zfill(2)


def classify_hla(norm):
    """
    Assign broad human HLA class from a normalized HLA string.
    """
    if not norm:
        return ""

    u = norm.upper()

    if re.match(r"^HLA-[ABCEFG]\*", u) or re.match(r"^HLA-[ABCEFG]$", u):
        return "I"

    if re.match(r"^HLA-(DRB[1-9]?|DQA1|DQB1|DPA1|DPB1)\*", u):
        return "II"

    if re.match(r"^HLA-(DR|DQ|DP)$", u):
        return "II"

    return ""


def is_nonhuman_mhc(x):
    u = clean(x).upper().replace(" ", "")
    return bool(re.search(r"(^H2|^H-2|H2-|H-2|H2K|H2D|H2IA|H2-IA|H2-IE|RT1|MAMU|PATR|BOLA|DLA|SLA-|SLA\*)", u))


def normalize_exact_or_lowres_hla(x):
    """
    Return:
      normalized_hla, hla_class, match_level, patch_reason

    match_level:
      exact_allele      = allele-level or already valid specific allele, e.g. HLA-E*01:03
      low_resolution   = family/serotype-like, e.g. HLA-Cw6 -> HLA-C*06
      broad_locus      = HLA-DR, HLA-DQ, HLA-DP, HLA-E with no allele
      nonhuman_mhc     = H2-Kb etc.; no normalized human HLA
      unresolved       = do not patch
    """
    raw = clean(x)
    if not raw:
        return "", "", "missing", "missing_original"

    s = raw.strip()
    u = s.upper().replace(" ", "")

    if is_nonhuman_mhc(s):
        return "", "nonhuman", "nonhuman_mhc", "nonhuman_mhc"

    if re.search(r"(NOTDETERMINED|UNKNOWN|UNDETERMINED|UNRESTRICTED|MHCCLASS|HLACLASS|CLASSI$|CLASSII$|MHC-I|MHC-II)", u):
        if "II" in u:
            return "", "II", "broad_or_ambiguous", "ambiguous_class_II"
        if "I" in u:
            return "", "I", "broad_or_ambiguous", "ambiguous_class_I"
        return "", "", "broad_or_ambiguous", "ambiguous_unknown"

    # Do not try to normalize complex alternatives here.
    if re.search(r"[/;,]", s) or re.search(r"\bOR\b|\bAND\b", s, flags=re.I):
        return "", "", "multiple_or_complex", "multiple_or_complex"

    # Normalize common separators
    t = u.replace("_", "-")

    # ------------------------------------------------------------
    # Exact or allele-like class I: A/B/C/E/F/G
    # Examples:
    #   HLA-E*01:03 -> HLA-E*01:03
    #   HLA-F:02:02 -> HLA-F*02:02
    #   A02:01 -> HLA-A*02:01
    #   HLA-A0201 -> HLA-A*02:01
    # ------------------------------------------------------------
    m = re.fullmatch(r"(?:HLA-)?([ABCEFG])\*?(\d{2})(?::?(\d{2}))?(?::?(\d{2}))?", t)
    if m:
        locus, f1, f2, f3 = m.group(1), m.group(2), m.group(3), m.group(4)
        if f1 and f2 and f3:
            return f"HLA-{locus}*{f1}:{f2}:{f3}", "I", "exact_allele", "class_I_exact_or_allele_like"
        if f1 and f2:
            return f"HLA-{locus}*{f1}:{f2}", "I", "exact_allele", "class_I_exact_or_allele_like"
        if f1:
            return f"HLA-{locus}*{f1}", "I", "low_resolution", "class_I_low_resolution"

    # HLA-Cw6, HLA-Cw4, HLA-Cw12
    m = re.fullmatch(r"HLA-CW(\d{1,2})", t)
    if m:
        return f"HLA-C*{pad2(m.group(1))}", "I", "low_resolution", "HLA_Cw_serotype"

    # HLA-A2, HLA-B7 style
    m = re.fullmatch(r"HLA-([AB])(\d{1,2})", t)
    if m:
        return f"HLA-{m.group(1)}*{pad2(m.group(2))}", "I", "low_resolution", "class_I_serotype"

    # Broad class I locus only: HLA-A, HLA-B, HLA-C, HLA-E, HLA-F, HLA-G
    m = re.fullmatch(r"(?:HLA-)?([ABCEFG])", t)
    if m:
        return "", "I", "broad_locus", f"broad_HLA_{m.group(1)}"

    # ------------------------------------------------------------
    # Exact or allele-like class II single-chain values
    # Examples:
    #   DRB1*01:01 -> HLA-DRB1*01:01
    #   HLA-DRB1*01:01 -> HLA-DRB1*01:01
    #   DQB1*03:01 -> HLA-DQB1*03:01
    # ------------------------------------------------------------
    m = re.fullmatch(r"(?:HLA-)?(DRB[1-9]?|DQA1|DQB1|DPA1|DPB1)\*?(\d{2})(?::?(\d{2}))?(?::?(\d{2}))?", t)
    if m:
        locus, f1, f2, f3 = m.group(1), m.group(2), m.group(3), m.group(4)
        if f1 and f2 and f3:
            return f"HLA-{locus}*{f1}:{f2}:{f3}", "II", "exact_allele", "class_II_exact_or_allele_like"
        if f1 and f2:
            return f"HLA-{locus}*{f1}:{f2}", "II", "exact_allele", "class_II_exact_or_allele_like"
        if f1:
            return f"HLA-{locus}*{f1}", "II", "low_resolution", "class_II_low_resolution"

    # HLA-DR1, HLA-DR3, HLA-DR15 etc.
    # Treat as DRB1 family-level only.
    m = re.fullmatch(r"HLA-DR(\d{1,2})", t)
    if m:
        return f"HLA-DRB1*{pad2(m.group(1))}", "II", "low_resolution", "HLA_DR_serotype"

    # HLA-DQ2, HLA-DQ8 etc.
    # Do not assign exact DQA/DQB pair; keep as broad/serotype class II.
    m = re.fullmatch(r"HLA-DQ(\d{1,2})", t)
    if m:
        return f"HLA-DQ*{pad2(m.group(1))}", "II", "low_resolution", "HLA_DQ_serotype"

    # HLA-DPw4 etc.
    # Keep as low-resolution DPw label, not an exact DPB1 allele.
    m = re.fullmatch(r"HLA-DPW(\d{1,2})", t)
    if m:
        return f"HLA-DP*w{m.group(1)}", "II", "low_resolution", "HLA_DPw_serotype"

    # Broad class II locus only.
    if t in {"HLA-DR", "DR"}:
        return "", "II", "broad_locus", "broad_HLA_DR"
    if t in {"HLA-DQ", "DQ"}:
        return "", "II", "broad_locus", "broad_HLA_DQ"
    if t in {"HLA-DP", "DP"}:
        return "", "II", "broad_locus", "broad_HLA_DP"

    return "", "", "unresolved", "unresolved"


def patch_file(infile, outfile):
    stats = Counter()
    details = Counter()

    with gzip.open(infile, "rt", newline="") as fin, gzip.open(outfile, "wt", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        old_fields = reader.fieldnames

        if old_fields is None:
            raise RuntimeError(f"No header in {infile}")

        required = {"hla_original", "hla_normalized", "hla_class"}
        missing = required - set(old_fields)
        if missing:
            raise RuntimeError(f"{infile} missing required columns: {missing}")

        new_fields = list(old_fields)
        for extra in ["hla_match_level", "hla_patch_reason", "hla_patch_changed"]:
            if extra not in new_fields:
                new_fields.append(extra)

        writer = csv.DictWriter(fout, delimiter="\t", fieldnames=new_fields, lineterminator="\n")
        writer.writeheader()

        for row in reader:
            stats["rows"] += 1

            orig = clean(row.get("hla_original", ""))
            old_norm = clean(row.get("hla_normalized", ""))
            old_class = clean(row.get("hla_class", ""))

            row["hla_match_level"] = ""
            row["hla_patch_reason"] = ""
            row["hla_patch_changed"] = "0"

            if old_norm:
                # Existing normalized HLA remains untouched.
                row["hla_match_level"] = "preexisting_normalized"
                row["hla_patch_reason"] = "already_normalized"
                if not old_class:
                    inferred_class = classify_hla(old_norm)
                    if inferred_class:
                        row["hla_class"] = inferred_class
                        row["hla_patch_changed"] = "1"
                        stats["class_filled_from_existing_norm"] += 1
                stats["already_normalized"] += 1
                writer.writerow(row)
                continue

            new_norm, new_class, match_level, reason = normalize_exact_or_lowres_hla(orig)

            row["hla_match_level"] = match_level
            row["hla_patch_reason"] = reason

            changed = False

            if new_norm and not old_norm:
                row["hla_normalized"] = new_norm
                changed = True
                stats["newly_normalized"] += 1
                stats[f"newly_normalized::{match_level}"] += 1
                details[(reason, orig, new_norm, new_class, match_level)] += 1

            if new_class and not old_class:
                row["hla_class"] = new_class
                changed = True
                stats["newly_classified"] += 1
                stats[f"newly_classified::{new_class}"] += 1

            if not new_norm and not new_class:
                stats["still_unresolved"] += 1
                if orig:
                    details[(reason, orig, "", "", match_level)] += 1

            if match_level == "nonhuman_mhc":
                stats["nonhuman_mhc"] += 1
                details[(reason, orig, "", "nonhuman", match_level)] += 1

            if match_level == "broad_or_ambiguous":
                stats["broad_or_ambiguous"] += 1

            if match_level == "multiple_or_complex":
                stats["multiple_or_complex"] += 1

            row["hla_patch_changed"] = "1" if changed else "0"
            writer.writerow(row)

    return stats, details


def main():
    files = sorted(IN_DIR.glob("*.normalized.tsv.gz"))

    if not files:
        raise FileNotFoundError(f"No .normalized.tsv.gz files found in {IN_DIR}")

    summary_rows = []
    all_details = Counter()

    for infile in files:
        outfile = OUT_DIR / infile.name
        print(f"Patching {infile.name}")
        stats, details = patch_file(infile, outfile)

        row = {"file": infile.name}
        row.update(dict(stats))
        summary_rows.append(row)

        for key, n in details.items():
            reason, orig, norm, cls, level = key
            all_details[(infile.name, reason, orig, norm, cls, level)] += n

    # Write summary
    all_keys = sorted({k for row in summary_rows for k in row.keys() if k != "file"})
    with open(SUMMARY_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=["file"] + all_keys, lineterminator="\n")
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    # Write patch details
    with open(DETAIL_FILE, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            delimiter="\t",
            fieldnames=[
                "file",
                "patch_reason",
                "hla_original",
                "hla_normalized",
                "hla_class",
                "hla_match_level",
                "count",
            ],
            lineterminator="\n",
        )
        writer.writeheader()

        for (file, reason, orig, norm, cls, level), n in sorted(
            all_details.items(),
            key=lambda x: x[1],
            reverse=True,
        ):
            writer.writerow({
                "file": file,
                "patch_reason": reason,
                "hla_original": orig,
                "hla_normalized": norm,
                "hla_class": cls,
                "hla_match_level": level,
                "count": n,
            })

    print()
    print(f"Wrote patched tables to: {OUT_DIR}")
    print(f"Wrote summary: {SUMMARY_FILE}")
    print(f"Wrote details: {DETAIL_FILE}")


if __name__ == "__main__":
    main()
