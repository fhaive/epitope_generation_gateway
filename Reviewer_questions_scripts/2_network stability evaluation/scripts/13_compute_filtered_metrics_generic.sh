#!/usr/bin/env bash
set -euo pipefail

ROOT="/data/fsluma/pipelines/Epitope_Generation_Gateway"

if [[ -z "${OUT_BASE:-}" ]]; then
  echo "ERROR: Set OUT_BASE to the benchmark filtered output folder"
  exit 1
fi

IN_RDS="$OUT_BASE/Sample_Specific_Networks_PPI_filtered/filtered_networks_rds"

DEG_DIR="$OUT_BASE/Network_Metrics_Degree"
STR_DIR="$OUT_BASE/Network_Metrics_Strength"
BET_DIR="$OUT_BASE/Network_Metrics_Betweenness"
LCI_DIR="$OUT_BASE/Network_Metrics_LargestComponentImpact"
LOG_DIR="$OUT_BASE/logs_metrics"

METRIC_SCRIPT="$ROOT/scripts/final_scripts/R_scripts/compute_network_metrics.R"
LCI_SCRIPT="$ROOT/scripts/final_scripts/R_scripts/compute_largest_component_impact.R"

mkdir -p "$DEG_DIR" "$STR_DIR" "$BET_DIR" "$LCI_DIR" "$LOG_DIR"

run_one() {
  set -euo pipefail

  sample="$1"

  in_rds="$IN_RDS/${sample}_filtered.rds"

  degree="$DEG_DIR/${sample}_degree.tsv"
  strength="$STR_DIR/${sample}_strength.tsv"
  betweenness="$BET_DIR/${sample}_betweenness.tsv"
  impact="$LCI_DIR/${sample}_impact.tsv"

  log="$LOG_DIR/${sample}"

  if [[ -s "$degree" && -s "$strength" && -s "$betweenness" && -s "$impact" ]]; then
    echo "[SKIP] $sample metrics already exist"
    exit 0
  fi

  if [[ ! -s "$in_rds" ]]; then
    echo "[FAIL] Missing filtered RDS: $in_rds"
    exit 1
  fi

  echo "[RUN] $sample"

  if [[ ! -s "$degree" || ! -s "$strength" || ! -s "$betweenness" ]]; then
    if ! Rscript "$METRIC_SCRIPT" \
      "$in_rds" \
      "$betweenness" \
      "$degree" \
      "$strength" \
      > "${log}.compute_network_metrics.log" 2>&1; then
      echo "[FAIL] compute_network_metrics failed for $sample"
      tail -20 "${log}.compute_network_metrics.log"
      exit 1
    fi
  fi

  if [[ ! -s "$impact" ]]; then
    if ! Rscript "$LCI_SCRIPT" \
      "$in_rds" \
      "$impact" \
      > "${log}.lci.log" 2>&1; then
      echo "[FAIL] LCI failed for $sample"
      tail -20 "${log}.lci.log"
      exit 1
    fi
  fi

  echo "[DONE] $sample"
}

export -f run_one
export ROOT OUT_BASE IN_RDS DEG_DIR STR_DIR BET_DIR LCI_DIR LOG_DIR METRIC_SCRIPT LCI_SCRIPT

find "$IN_RDS" -maxdepth 1 -name "Sample_*_filtered.rds" \
  | sed 's|.*/||; s|_filtered\.rds$||' \
  | sort \
  | xargs -I{} -P "${N_JOBS:-8}" bash -c 'set -euo pipefail; run_one "$@"' _ {}
