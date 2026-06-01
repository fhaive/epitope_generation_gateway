#!/usr/bin/env python3

import json
import os
import sys
import time
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

BASE_URL = "https://api.gdc.cancer.gov"
PROJECT = "TCGA-SKCM"

OUT_DIR = Path("gdc_tcga_skcm_star_counts")
FILES_DIR = OUT_DIR / "files"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FILES_DIR.mkdir(parents=True, exist_ok=True)

MANIFEST_TSV = OUT_DIR / "gdc_tcga_skcm_star_counts_manifest.tsv"

FILTERS = {
    "op": "and",
    "content": [
        {
            "op": "in",
            "content": {
                "field": "cases.project.project_id",
                "value": [PROJECT],
            },
        },
        {
            "op": "in",
            "content": {
                "field": "files.data_category",
                "value": ["Transcriptome Profiling"],
            },
        },
        {
            "op": "in",
            "content": {
                "field": "files.data_type",
                "value": ["Gene Expression Quantification"],
            },
        },
        {
            "op": "in",
            "content": {
                "field": "files.analysis.workflow_type",
                "value": ["STAR - Counts"],
            },
        },
        {
            "op": "in",
            "content": {
                "field": "files.access",
                "value": ["open"],
            },
        },
    ],
}

FIELDS = [
    "file_id",
    "file_name",
    "data_type",
    "data_category",
    "experimental_strategy",
    "analysis.workflow_type",
    "cases.submitter_id",
    "cases.case_id",
    "cases.samples.submitter_id",
    "cases.samples.sample_id",
    "cases.samples.sample_type",
    "cases.samples.tissue_type",
    "cases.samples.portions.analytes.aliquots.submitter_id",
]


def post_files(size: int = 2000) -> dict:
    payload = {
        "filters": FILTERS,
        "fields": ",".join(FIELDS),
        "format": "JSON",
        "size": str(size),
    }
    r = requests.post(f"{BASE_URL}/files", json=payload, timeout=120)
    r.raise_for_status()
    return r.json()


def parse_hits(resp: dict) -> pd.DataFrame:
    rows = []
    hits = resp["data"]["hits"]

    for h in hits:
        cases = h.get("cases", [])
        case_submitter_id = None
        case_id = None
        sample_submitter_id = None
        sample_type = None
        tissue_type = None
        aliquot_submitter_id = None

        if cases:
            c = cases[0]
            case_submitter_id = c.get("submitter_id")
            case_id = c.get("case_id")
            samples = c.get("samples", [])
            if samples:
                s = samples[0]
                sample_submitter_id = s.get("submitter_id")
                sample_type = s.get("sample_type")
                tissue_type = s.get("tissue_type")

                portions = s.get("portions", [])
                if portions:
                    analytes = portions[0].get("analytes", [])
                    if analytes:
                        aliquots = analytes[0].get("aliquots", [])
                        if aliquots:
                            aliquot_submitter_id = aliquots[0].get("submitter_id")

        analysis = h.get("analysis") or {}

        rows.append(
            {
                "file_id": h.get("file_id"),
                "file_name": h.get("file_name"),
                "case_submitter_id": case_submitter_id,
                "case_id": case_id,
                "sample_submitter_id": sample_submitter_id,
                "sample_type": sample_type,
                "tissue_type": tissue_type,
                "aliquot_submitter_id": aliquot_submitter_id,
                "data_category": h.get("data_category"),
                "data_type": h.get("data_type"),
                "experimental_strategy": h.get("experimental_strategy"),
                "workflow_type": analysis.get("workflow_type"),
            }
        )

    return pd.DataFrame(rows)


def download_one(file_id: str, file_name: str) -> Path:
    out = FILES_DIR / file_name
    if out.exists() and out.stat().st_size > 0:
        return out

    url = f"{BASE_URL}/data/{file_id}"

    # Retry because large public GDC downloads sometimes transiently fail.
    for attempt in range(1, 6):
        try:
            with requests.get(url, stream=True, timeout=300) as r:
                r.raise_for_status()
                total = int(r.headers.get("content-length", 0))
                tmp = out.with_suffix(out.suffix + ".tmp")

                with open(tmp, "wb") as fh, tqdm(
                    total=total,
                    unit="B",
                    unit_scale=True,
                    unit_divisor=1024,
                    desc=file_name[:60],
                    leave=False,
                ) as pbar:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            fh.write(chunk)
                            pbar.update(len(chunk))

                tmp.rename(out)
                return out

        except Exception as e:
            print(f"[WARN] download failed attempt {attempt}/5 for {file_id}: {e}", file=sys.stderr)
            time.sleep(10 * attempt)

    raise RuntimeError(f"Failed to download {file_id} after 5 attempts")


def main():
    print(f"Querying GDC for {PROJECT} STAR-count files...")
    resp = post_files(size=2000)
    total = resp["data"]["pagination"]["total"]
    print(f"GDC query returned {total} files")

    if total > 2000:
        raise RuntimeError("More than 2000 files returned; implement pagination")

    manifest = parse_hits(resp)
    manifest = manifest.sort_values(["case_submitter_id", "sample_submitter_id", "file_name"])
    manifest.to_csv(MANIFEST_TSV, sep="\t", index=False)

    print("\nSample type counts:")
    print(manifest["sample_type"].value_counts(dropna=False).to_string())

    print(f"\nManifest written: {MANIFEST_TSV}")
    print("\nDownloading files...")

    local_paths = []
    for _, row in tqdm(manifest.iterrows(), total=len(manifest), desc="files"):
        path = download_one(row["file_id"], row["file_name"])
        local_paths.append(str(path))

    manifest["local_path"] = local_paths
    manifest.to_csv(MANIFEST_TSV, sep="\t", index=False)

    print("\nDone.")
    print(f"Manifest: {MANIFEST_TSV}")
    print(f"Files dir: {FILES_DIR}")
    print(f"Downloaded files: {len(local_paths)}")


if __name__ == "__main__":
    main()
