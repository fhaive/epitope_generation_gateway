#!/usr/bin/env python3
"""
build_epitope_table.py

Generate standalone interactive HTML table(s) from epitope CSVs.

Features
- Sortable/searchable table via DataTables (CDN assets).
- `Intogen_Driver_role` rendered with hover tooltips:
    Act  -> Activating mutation (oncogene-like)
    LoF  -> Loss of Function (tumor suppressor-like)
    Act;LoF -> shows both explanations
- Compact hover badges for:
    GO_MF (MF), GO_BP (BP), GO_CC (CC), Cancer_Hallmark (CH),
    Intogen_Cohort (Coh), Intogen_CancerType (CT)
- If Mutation.Source == "Somatic", display it as "Genomic".
- Single-file mode or batch mode over a directory.

Usage
------
Single file:
    python build_epitope_table.py \
      -i /path/to/Sample_TCGA-XYZ_epitopes_final.csv \
      -o /path/to/Sample_TCGA-XYZ_epitopes.html

Batch (current directory, default pattern "*_epitopes_final.csv"):
    python build_epitope_table.py -i .

Batch with custom folder and/or glob:
    python build_epitope_table.py -i /data/.../final_epitopes --glob "*_epitopes_final.csv" --outdir /data/.../final_epitopes
"""
import argparse
import glob
import html
import os
import re
from typing import Iterable, List, Optional

import pandas as pd


# ---- Tooltip helpers ---------------------------------------------------------

DRIVER_EXPLAIN = {
    "Act": "Activating mutation (oncogene-like): mutation confers gain-of-function or increased activity.",
    "LoF": "Loss of Function (tumor suppressor-like): mutation likely inactivates or reduces the gene's function.",
}


def explain_driver_role(value: str) -> str:
    """Return an HTML <span> with title tooltip for Intogen_Driver_role."""
    if value is None:
        return ""
    raw = str(value).strip()
    if raw == "" or raw.lower() == "nan":
        return ""
    # Split on semicolons/commas/whitespace
    toks = [t for t in re.split(r"[;,\s]+", raw) if t]
    seen = []
    for t in toks:
        key = "Act" if t.lower().startswith("act") else ("LoF" if t.lower().startswith("lof") else t)
        if key not in seen:
            seen.append(key)
    # Build explanation lines for any recognized tokens; keep unknowns as-is
    lines: List[str] = []
    for t in seen:
        if t in DRIVER_EXPLAIN:
            lines.append(f"{t}: {DRIVER_EXPLAIN[t]}")
        else:
            lines.append(t)
    title = " \n".join(lines)
    return f'<span class="driver" title="{html.escape(title, quote=True)}">{html.escape(raw)}</span>'


def badge(short: str, full_text: str) -> str:
    """Return a compact badge with tooltip text, or empty if no value."""
    if full_text is None:
        return ""
    raw = str(full_text).strip()
    if raw == "" or raw.lower() == "nan":
        return ""
    return f'<span class="pill" title="{html.escape(raw, quote=True)}">{html.escape(short)}</span>'


# ---- Core renderer -----------------------------------------------------------

def render_html_table(df: pd.DataFrame, title: str = "Epitope Prioritisation Table") -> str:
    """Return a complete, standalone HTML document string for the given DataFrame."""
    df = df.fillna("")

    # Column resolver (be tolerant to minor variations)
    col = lambda name: name if name in df.columns else None

    driver_col = col("Intogen_Driver_role")

    go_mf_col = col("GO_MF") or col("GO.MF") or col("GO-MF")
    go_bp_col = col("GO_BP") or col("GO.BP") or col("GO-BP")
    go_cc_col = col("GO_CC") or col("GO.CC") or col("GO-CC")
    hallmark_col = col("Cancer_Hallmark") or col("Cancer.Hallmark") or col("Cancer-Hallmark")

    cohort_col = col("Intogen_Cohort") or col("Intogen.Cohort") or col("Intogen-Cohort")
    cancertype_col = col("Intogen_CancerType") or col("Intogen.CancerType") or col("Intogen-CancerType")

    mutation_source_col = (
        col("Mutation.Source") or col("Mutation_Source") or col("Mutation Source")
    )

    render_df = df.copy()

    # Normalize Mutation Source values
    if mutation_source_col:
        render_df[mutation_source_col] = render_df[mutation_source_col].replace("Somatic", "Genomic")

    # Driver role tooltip
    if driver_col:
        render_df[driver_col] = render_df[driver_col].apply(explain_driver_role)

    # Compact badges for GO & Hallmark columns
    if go_mf_col:
        render_df[go_mf_col] = render_df[go_mf_col].apply(lambda x: badge("MF", x))
    if go_bp_col:
        render_df[go_bp_col] = render_df[go_bp_col].apply(lambda x: badge("BP", x))
    if go_cc_col:
        render_df[go_cc_col] = render_df[go_cc_col].apply(lambda x: badge("CC", x))
    if hallmark_col:
        render_df[hallmark_col] = render_df[hallmark_col].apply(lambda x: badge("CH", x))

    # Compact badges for IntOGen cohort & cancer type
    if cohort_col:
        render_df[cohort_col] = render_df[cohort_col].apply(lambda x: badge("Coh", x))
    if cancertype_col:
        render_df[cancertype_col] = render_df[cancertype_col].apply(lambda x: badge("CT", x))

    # Build HTML table
    headers: List[str] = list(render_df.columns)

    # Columns containing intentional HTML fragments (do not escape again)
    html_cols = {
        c for c in [driver_col, go_mf_col, go_bp_col, go_cc_col, hallmark_col, cohort_col, cancertype_col] if c
    }

    def render_cell(colname: str, value: str) -> str:
        if colname in html_cols:
            return value  # already escaped appropriately inside helper functions
        return html.escape("" if value is None else str(value))

    # Create table rows
    rows_html: List[str] = []
    for _, row in render_df.iterrows():
        tds = "".join(f"<td>{render_cell(col, row[col])}</td>" for col in headers)
        rows_html.append(f"<tr>{tds}</tr>")

    # Legend pieces (show only what exists)
    badges_legend = []
    if go_mf_col: badges_legend.append('<span class="pill" title="Gene Ontology: Molecular Function terms">MF</span> Molecular Function')
    if go_bp_col: badges_legend.append('<span class="pill" title="Gene Ontology: Biological Process terms">BP</span> Biological Process')
    if go_cc_col: badges_legend.append('<span class="pill" title="Gene Ontology: Cellular Component terms">CC</span> Cellular Component')
    if hallmark_col: badges_legend.append('<span class="pill" title="MSigDB Hallmark pathway set(s)">CH</span> Cancer Hallmark')
    if cohort_col: badges_legend.append('<span class="pill" title="IntOGen cohort(s)">Coh</span> IntOGen Cohort')
    if cancertype_col: badges_legend.append('<span class="pill" title="IntOGen cancer type(s)">CT</span> IntOGen Cancer Type')

    legend_html = " &nbsp;|&nbsp; ".join(badges_legend)

    html_doc = f"""
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>{html.escape(title)}</title>

<!-- DataTables (CDN) -->
<link rel="stylesheet" href="https://cdn.datatables.net/v/dt/dt-1.13.8/datatables.min.css"/>
<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
<script src="https://cdn.datatables.net/v/dt/dt-1.13.8/datatables.min.js"></script>

<style>
  body {{ font-family: system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif; margin: 16px; }}
  h1 {{ font-size: 20px; margin: 0 0 10px 0; }}
  .legend {{ font-size: 12px; color: #555; margin-bottom: 12px; }}
  .legend .pill {{ margin-right: 6px; }}
  .pill {{
    display: inline-block;
    padding: 2px 6px;
    border-radius: 9999px;
    border: 1px solid #cfe1f3;
    background: #eaf3ff;
    font-size: 11px;
    font-weight: 600;
    line-height: 1.4;
    cursor: help;
    white-space: nowrap;
  }}
  .driver {{ cursor: help; border-bottom: 1px dotted #888; }}
  table.dataTable thead th {{ white-space: nowrap; }}
  table.dataTable tbody td {{ vertical-align: top; }}
  .table-wrap {{ overflow-x: auto; }}
</style>
</head>
<body>
  <h1>{html.escape(title)}</h1>
  <div class="legend">
    <strong>Hover tips:</strong>
    {legend_html}
    &nbsp;|&nbsp;
    <span class="driver" title="Act: Activating mutation (oncogene-like)
LoF: Loss of Function (tumor suppressor-like)">Intogen_Driver_role</span>
  </div>

  <div class="table-wrap">
    <table id="epitopes" class="display compact" style="width:100%">
      <thead>
        <tr>
          {"".join(f"<th>{html.escape(h)}</th>" for h in headers)}
        </tr>
      </thead>
      <tbody>
        {"".join(rows_html)}
      </tbody>
    </table>
  </div>

<script>
$(document).ready(function() {{
  $('#epitopes').DataTable({{
    pageLength: 25,
    lengthMenu: [[10, 25, 50, 100, -1], [10, 25, 50, 100, "All"]],
    deferRender: true,
    autoWidth: false,
    scrollX: true,
    order: []
  }});
}});
</script>
</body>
</html>
"""
    return html_doc


# ---- Batch wrapper -----------------------------------------------------------

def build_one(input_csv: str, output_html: Optional[str] = None, title: Optional[str] = None) -> str:
    """Build one HTML file from a CSV. Returns output path."""
    if title is None:
        title = "Epitope Prioritisation Table"

    df = pd.read_csv(input_csv, dtype=str)

    html_doc = render_html_table(df, title=title)

    if output_html is None:
        base = os.path.basename(input_csv)
        # default: replace trailing "_epitopes_final.csv" with "_epitopes.html"
        if base.endswith("_epitopes_final.csv"):
            base = base.replace("_epitopes_final.csv", "_epitopes.html")
        else:
            base = os.path.splitext(base)[0] + ".html"
        output_html = os.path.join(os.path.dirname(os.path.abspath(input_csv)), base)

    os.makedirs(os.path.dirname(os.path.abspath(output_html)), exist_ok=True)
    with open(output_html, "w", encoding="utf-8") as f:
        f.write(html_doc)
    return output_html


def iter_csvs(path: str, pattern: str) -> List[str]:
    """Yield CSV paths for batch mode."""
    if os.path.isdir(path):
        return sorted(glob.glob(os.path.join(path, pattern)))
    elif os.path.isfile(path):
        return [path]
    else:
        # Treat as glob relative to CWD
        return sorted(glob.glob(path))


def main():
    ap = argparse.ArgumentParser(description="Build standalone interactive HTML table(s) from epitope CSVs.")
    ap.add_argument("-i", "--input", required=True,
                    help="Path to a CSV file OR a directory OR a glob pattern.")
    ap.add_argument("-o", "--output", default=None,
                    help="Output HTML path (single-file mode only).")
    ap.add_argument("--title", default=None, help="Page title (single-file mode or applied to all).")
    ap.add_argument("--glob", default="*_epitopes_final.csv",
                    help="When INPUT is a directory, CSV filename pattern to match (default: *_epitopes_final.csv).")
    ap.add_argument("--outdir", default=None,
                    help="When batching, directory to write HTML files (defaults to same dir as CSVs).")
    args = ap.parse_args()

    csvs = iter_csvs(args.input, args.glob)

    if not csvs:
        raise SystemExit("No CSVs found for the given input.")

    # Single-file mode
    if len(csvs) == 1 and os.path.isfile(csvs[0]) and (args.output is not None or not os.path.isdir(args.input)):
        out = args.output
        if out is None:
            out = None  # build_one computes default
        final = build_one(csvs[0], output_html=out, title=args.title)
        print(f"✔ Wrote: {final}")
        return

    # Batch mode
    for csv_path in csvs:
        if args.outdir:
            base = os.path.basename(csv_path)
            if base.endswith("_epitopes_final.csv"):
                base = base.replace("_epitopes_final.csv", "_epitopes.html")
            else:
                base = os.path.splitext(base)[0] + ".html"
            out_path = os.path.join(args.outdir, base)
        else:
            out_path = None  # build_one computes default next to CSV
        final = build_one(csv_path, output_html=out_path, title=args.title)
        print(f"✔ Wrote: {final}")


if __name__ == "__main__":
    main()
