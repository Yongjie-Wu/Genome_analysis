#!/usr/bin/env bash
set -euo pipefail

usage() {
cat <<'EOF'
Usage:
  bash 00.index_reference.sh <genome.fa> <ref_prefix> [project_dir]

Example:
  bash 00.index_reference.sh hap1/Pyrps_dianli_hap1.fa dianli_hap1 .

Output:
  00.index/<ref_prefix>.mmi
  <genome.fa>.fai
EOF
}

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
  usage
  exit 1
fi

GENOME="$1"
REF_PREFIX="$2"
PROJECT_DIR="${3:-.}"

if [ ! -s "$GENOME" ]; then
  echo "ERROR: genome fasta not found: $GENOME" >&2
  exit 1
fi

mkdir -p "${PROJECT_DIR}/00.index" "${PROJECT_DIR}/logs"

command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found" >&2; exit 1; }
command -v pbmm2 >/dev/null 2>&1 || { echo "ERROR: pbmm2 not found" >&2; exit 1; }

GENOME_ABS="$(readlink -f "$GENOME")"
MMI="${PROJECT_DIR}/00.index/${REF_PREFIX}.mmi"

echo "[00] samtools faidx: ${GENOME_ABS}"
if [ ! -s "${GENOME_ABS}.fai" ]; then
  samtools faidx "$GENOME_ABS"
else
  echo "[00] fasta index exists: ${GENOME_ABS}.fai"
fi

echo "[00] pbmm2 index: ${MMI}"
if [ ! -s "$MMI" ]; then
  pbmm2 index "$GENOME_ABS" "$MMI" > "${PROJECT_DIR}/logs/00.${REF_PREFIX}.pbmm2_index.log" 2>&1
else
  echo "[00] pbmm2 index exists: ${MMI}"
fi

echo "[00] done"
