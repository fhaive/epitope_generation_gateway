#!/usr/bin/env bash
set -euo pipefail

ROOT="/data/fsluma/pipelines/Epitope_Generation_Gateway"
WORK="$ROOT/TCGA_melanoma/sample_specific_networks/reviewer_C_external_benchmark"

IN_DIR="$WORK/tcga_primary_benchmark_lioness_wci/Sample_Specific_Networks"

OUT_BASE="$WORK/tcga_primary_benchmark_filtered_metrics"
OUT_RDS="$OUT_BASE/Sample_Specific_Networks_PPI_filtered/filtered_networks_rds"
OUT_MTX="$OUT_BASE/Sample_Specific_Networks_PPI_filtered/filtered_networks_matrix"
LOG_DIR="$OUT_BASE/logs_filtering"

CONFIG="$ROOT/scripts/final_scripts/config/filter_config.yaml"
HUMANNET="$ROOT/TCGA_melanoma/sample_specific_networks/gene_lists/HumanNet_ensembl.rds"
MEMBRANE="$ROOT/TCGA_melanoma/sample_specific_networks/gene_lists/membrane_ensembl.rds"
FILTER_SCRIPT="$ROOT/scripts/final_scripts/R_scripts/Filter_Sample_Networks.R"

mkdir -p "$OUT_RDS" "$OUT_MTX" "$LOG_DIR"

run_one() {
  set -euo pipefail

  sample="$1"

  in_rds="$IN_DIR/${sample}.rds"
  out_rds="$OUT_RDS/${sample}_filtered.rds"
  out_mtx="$OUT_MTX/${sample}.mtx"
  out_genes="$OUT_MTX/${sample}_genes.txt"
  log="$LOG_DIR/${sample}.log"

  if [[ -s "$out_rds" && -s "$out_mtx" && -s "$out_genes" ]]; then
    echo "[SKIP] $sample already filtered"
    exit 0
  fi

  echo "[RUN] $sample"

  if [[ ! -s "$in_rds" ]]; then
    echo "[FAIL] Missing input RDS: $in_rds" | tee "$log"
    exit 1
  fi

  if ! Rscript "$FILTER_SCRIPT" \
    "$in_rds" \
    "$out_rds" \
    "$out_mtx" \
    "$out_genes" \
    "$CONFIG" \
    "$HUMANNET" \
    "$MEMBRANE" \
    > "$log" 2>&1; then

    echo "[FAIL] $sample filtering failed. See: $log"
    tail -20 "$log"
    exit 1
  fi

  if [[ ! -s "$out_rds" ]]; then
    echo "[FAIL] $sample finished but output RDS missing: $out_rds"
    tail -20 "$log"
    exit 1
  fi

  echo "[DONE] $sample"
}

export -f run_one
export ROOT WORK IN_DIR OUT_RDS OUT_MTX LOG_DIR CONFIG HUMANNET MEMBRANE FILTER_SCRIPT

find "$IN_DIR" -maxdepth 1 -name "Sample_*.rds" \
  | sed 's|.*/||; s|\.rds$||' \
  | sort \
  | xargs -I{} -P "${N_JOBS:-4}" bash -c 'set -euo pipefail; run_one "$@"' _ {}
