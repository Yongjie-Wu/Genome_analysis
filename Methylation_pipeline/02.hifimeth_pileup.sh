#!/usr/bin/env bash
set -euo pipefail

usage() {
cat <<'EOF'
Usage:
  bash 02.hifimeth_pileup.sh <sample_name> <genome.fa> <aligned.pbmm2.bam> [project_dir]

Example:
  bash 02.hifimeth_pileup.sh dianli_hap1 dianli/hap1/Pyrps_dianli_hap1.fa 01.align/dianli_hap1.pbmm2.bam .

Output:
  02.methylation/<sample_name>.CpG.cov.bed
  02.methylation/<sample_name>.CHG.cov.bed
  02.methylation/<sample_name>.CHH.cov.bed
EOF
}

if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
  usage
  exit 1
fi

SAMPLE="$1"
GENOME="$2"
BAM="$3"
PROJECT_DIR="${4:-.}"

for f in "$GENOME" "$BAM"; do
  if [ ! -s "$f" ]; then
    echo "ERROR: missing file: $f" >&2
    exit 1
  fi
done

mkdir -p "${PROJECT_DIR}/02.methylation" "${PROJECT_DIR}/logs"

if command -v hifimeth >/dev/null 2>&1; then
  HIFIMETH="$(command -v hifimeth)"
elif [ -x /storage3/wuyongjie/software/hifimeth_1.1.10_Linux-amd64/bin/hifimeth ]; then
  HIFIMETH="/storage3/wuyongjie/software/hifimeth_1.1.10_Linux-amd64/bin/hifimeth"
else
  echo "ERROR: hifimeth not found" >&2
  exit 1
fi

OUT_PREFIX="${PROJECT_DIR}/02.methylation/${SAMPLE}"

echo "[02] hifimeth pileup sample=${SAMPLE}"
if [ ! -s "${OUT_PREFIX}.CpG.cov.bed" ] || [ ! -s "${OUT_PREFIX}.CHG.cov.bed" ] || [ ! -s "${OUT_PREFIX}.CHH.cov.bed" ]; then
  (
    cd "${PROJECT_DIR}/02.methylation"
    "$HIFIMETH" pileup "$(readlink -f "$GENOME")" "$(readlink -f "$BAM")" "$SAMPLE" \
      > "../logs/02.${SAMPLE}.hifimeth_pileup.log" 2>&1
  )
else
  echo "[02] methylation cov.bed files exist"
fi

echo "[02] done"
