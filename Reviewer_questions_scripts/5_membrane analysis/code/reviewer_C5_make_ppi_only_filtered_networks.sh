#!/usr/bin/env bash
set -euo pipefail

ROOT="/data/fsluma/pipelines/Epitope_Generation_Gateway"

IN_DIR="$ROOT/TCGA_melanoma/sample_specific_networks/Sample_Specific_Networks"

OUT_BASE="$ROOT/TCGA_melanoma/rescue_final_analysis/reviewer_C5_membrane_impact/ppi_only_filtered_networks"
OUT_RDS="$OUT_BASE/filtered_networks_rds"
OUT_MTX="$OUT_BASE/filtered_networks_matrix"
LOG_DIR="$OUT_BASE/logs_filtering"

CONFIG="$ROOT/scripts/final_scripts/config/filter_config.yaml"
HUMANNET="$ROOT/TCGA_melanoma/sample_specific_networks/gene_lists/HumanNet_ensembl.rds"
EMPTY_MEMBRANE="$ROOT/TCGA_melanoma/rescue_final_analysis/reviewer_C5_membrane_impact/empty_membrane_ensembl.rds"
FILTER_SCRIPT="$ROOT/scripts/final_scripts/R_scripts/Filter_Sample_Networks.R"

mkdir -p "$OUT_RDS" "$OUT_MTX" "$LOG_DIR"

run_one() {
  set -euo pipefail

  sample="$1"

  in_rds="$IN_DIR/${sample}.rds"
  out_rds="$OUT_RDS/${sample}_ppi_only_filtered.rds"
  out_mtx="$OUT_MTX/${sample}_ppi_only.mtx"
  out_genes="$OUT_MTX/${sample}_ppi_only_genes.txt"
  log="$LOG_DIR/${sample}.log"

  if [[ -s "$out_rds" && -s "$out_mtx" && -s "$out_genes" ]]; then
    echo "[SKIP] $sample"
    exit 0
  fi

  echo "[RUN] $sample"

  if [[ ! -s "$in_rds" ]]; then
    echo "[FAIL] Missing input network: $in_rds" | tee "$log"
    exit 1
  fi

  if ! Rscript "$FILTER_SCRIPT" \
    "$in_rds" \
    "$out_rds" \
    "$out_mtx" \
    "$out_genes" \
    "$CONFIG" \
    "$HUMANNET" \
    "$EMPTY_MEMBRANE" \
    > "$log" 2>&1; then

    echo "[FAIL] $sample filtering failed. See: $log"
    tail -30 "$log"
    exit 1
  fi

  if [[ ! -s "$out_rds" ]]; then
    echo "[FAIL] $sample finished but output RDS missing: $out_rds"
    tail -30 "$log"
    exit 1
  fi

  echo "[DONE] $sample"
}

export -f run_one
export ROOT IN_DIR OUT_RDS OUT_MTX LOG_DIR CONFIG HUMANNET EMPTY_MEMBRANE FILTER_SCRIPT

find "$IN_DIR" -maxdepth 1 -name "Sample_*.rds" \
  | sed 's|.*/||; s|\.rds$||' \
  | sort \
  | xargs -I{} -P "${N_JOBS:-8}" bash -c 'set -euo pipefail; run_one "$@"' _ {}
