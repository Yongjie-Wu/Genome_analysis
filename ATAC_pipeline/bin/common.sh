#!/usr/bin/env bash
set -euo pipefail

usage_common() {
  cat >&2 <<EOF
Usage: bash $0 <config.sh>

The config.sh file stores all input paths and parameters.
See config.template.sh for required variables.
EOF
}

load_config() {
  if [ "$#" -lt 1 ]; then usage_common; exit 1; fi
  CONFIG="$1"
  [ -s "$CONFIG" ] || { echo "ERROR: missing config: $CONFIG" >&2; exit 1; }
  # shellcheck source=/dev/null
  source "$CONFIG"
  : "${PROJECT_NAME:?PROJECT_NAME is required}"
  : "${PROJECT_DIR:?PROJECT_DIR is required}"
  : "${SAMPLE_TABLE:?SAMPLE_TABLE is required}"
  : "${HAP1_FA:?HAP1_FA is required}"
  : "${HAP2_FA:?HAP2_FA is required}"
  : "${HAP1_GFF:?HAP1_GFF is required}"
  : "${HAP2_GFF:?HAP2_GFF is required}"
  THREADS="${THREADS:-24}"
  SORT_THREADS="${SORT_THREADS:-12}"
  FASTP_THREADS="${FASTP_THREADS:-8}"
  BIGWIG_THREADS="${BIGWIG_THREADS:-16}"
  SUMMARY_BIN_SIZE="${SUMMARY_BIN_SIZE:-10000}"
  PEAK_WINDOW="${PEAK_WINDOW:-3000}"
  BIGWIG_BIN_SIZE="${BIGWIG_BIN_SIZE:-10}"
  WINDOW_SIZE="${WINDOW_SIZE:-50000}"
  GENOME_SIZE="${GENOME_SIZE:-5e8}"
  WORK_ENV="${WORK_ENV:-work}"
  MACS2_ENV="${MACS2_ENV:-macs2_clean}"
  DEEPTOOLS_ENV="${DEEPTOOLS_ENV:-deeptools_local}"
  PY_ENV="${PY_ENV:-base}"
}

init_project() {
  mkdir -p "$PROJECT_DIR"/{00.clean,01.index,02.bam,03.bam_filter,04.peak/hap1,04.peak/hap2,05.bigwig,06.qc,07.annotation,08.cen,09.atac_windows,logs}
}

log() { echo "[$(date '+%F %T')] $*"; }
need_file() { [ -s "$1" ] || { echo "ERROR: missing or empty file: $1" >&2; exit 1; }; }
gzip_head_ok(){ [ "$(od -An -tx1 -N2 "$1" 2>/dev/null | tr -d ' \n')" = "1f8b" ]; }

activate_env() {
  local env="$1"
  export PYTHONNOUSERSITE=1
  unset PYTHONPATH
  # shellcheck disable=SC1091
  source "${HOME}/miniconda3/etc/profile.d/conda.sh"
  conda activate "$env"
}

sample_names() {
  awk 'BEGIN{FS="\t"} NR>1 && $1!="" {print $1}' "$SAMPLE_TABLE"
}

sample_rows() {
  awk 'BEGIN{FS="\t"} NR>1 && $1!="" {print $1"\t"$2"\t"$3}' "$SAMPLE_TABLE"
}
