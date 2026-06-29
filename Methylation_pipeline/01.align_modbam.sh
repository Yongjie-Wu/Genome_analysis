#!/usr/bin/env bash
set -euo pipefail

usage() {
cat <<'EOF'
Usage:
  bash 01.align_modbam.sh <sample_name> <mod.bam> <genome.fa> <ref.mmi> [project_dir] [threads]

Example:
  bash 01.align_modbam.sh dianli_hap1 dianli/mod.bam dianli/hap1/Pyrps_dianli_hap1.fa 00.index/dianli_hap1.mmi . 48

Output:
  01.align/<sample_name>.pbmm2.bam
  01.align/<sample_name>.pbmm2.bam.bai
EOF
}

if [ "$#" -lt 4 ] || [ "$#" -gt 6 ]; then
  usage
  exit 1
fi

SAMPLE="$1"
MOD_BAM="$2"
GENOME="$3"
MMI="$4"
PROJECT_DIR="${5:-.}"
THREADS="${6:-48}"

for f in "$MOD_BAM" "$GENOME" "$MMI"; do
  if [ ! -s "$f" ]; then
    echo "ERROR: missing file: $f" >&2
    exit 1
  fi
done

mkdir -p "${PROJECT_DIR}/01.align" "${PROJECT_DIR}/logs"
command -v pbmm2 >/dev/null 2>&1 || { echo "ERROR: pbmm2 not found" >&2; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found" >&2; exit 1; }

BAM="${PROJECT_DIR}/01.align/${SAMPLE}.pbmm2.bam"

echo "[01] pbmm2 align sample=${SAMPLE}"
if [ ! -s "$BAM" ]; then
  pbmm2 align \
    --preset CCS \
    --sort \
    -j "$THREADS" \
    "$MMI" \
    "$MOD_BAM" \
    "$BAM" \
    > "${PROJECT_DIR}/logs/01.${SAMPLE}.pbmm2_align.log" 2>&1
else
  echo "[01] aligned BAM exists: ${BAM}"
fi

if [ ! -s "${BAM}.bai" ]; then
  samtools index -@ "$THREADS" "$BAM"
fi

echo "[01] done: ${BAM}"
