#!/usr/bin/env bash
set -euo pipefail

usage() {
cat <<'EOF'
Usage:
  bash 06.summary.sh <sample_name> [project_dir]

Example:
  bash 06.summary.sh dianli_hap1 .

Output:
  06.summary/<sample_name>.methylation_pipeline_outputs.txt
EOF
}

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
  usage
  exit 1
fi

SAMPLE="$1"
PROJECT_DIR="${2:-.}"
OUTDIR="${PROJECT_DIR}/06.summary"
mkdir -p "$OUTDIR"
OUT="${OUTDIR}/${SAMPLE}.methylation_pipeline_outputs.txt"

{
  echo "sample: ${SAMPLE}"
  echo "project_dir: $(readlink -f "$PROJECT_DIR")"
  echo "date: $(date)"
  echo
  echo "[01.align]"
  find "${PROJECT_DIR}/01.align" -maxdepth 1 -type f -name "${SAMPLE}*" -printf "%p\t%s bytes\n" 2>/dev/null | sort || true
  echo
  echo "[02.methylation]"
  find "${PROJECT_DIR}/02.methylation" -maxdepth 1 -type f -name "${SAMPLE}*.cov.bed" -printf "%p\t%s bytes\n" 2>/dev/null | sort || true
  echo
  echo "[03.windows]"
  find "${PROJECT_DIR}/03.windows" -maxdepth 1 -type f -name "${SAMPLE}*" -printf "%p\t%s bytes\n" 2>/dev/null | sort || true
  echo
  echo "[04.features]"
  find "${PROJECT_DIR}/04.features" -maxdepth 1 -type f -name "${SAMPLE}*" -printf "%p\t%s bytes\n" 2>/dev/null | sort || true
  echo
  echo "[05.plots]"
  find "${PROJECT_DIR}/05.plots" -maxdepth 1 -type f -name "${SAMPLE}*" -printf "%p\t%s bytes\n" 2>/dev/null | sort || true
} > "$OUT"

cat "$OUT"
echo "[06] done: ${OUT}"
