#!/usr/bin/env bash
set -euo pipefail

usage() {
cat <<'EOF'
Usage:
  bash 03.window_methylation.sh <sample_name> <genome.fa> <CpG.cov.bed> <CHG.cov.bed> <CHH.cov.bed> [window_size] [project_dir]

Example:
  bash 03.window_methylation.sh dianli_hap1 dianli/hap1/Pyrps_dianli_hap1.fa 02.methylation/dianli_hap1.CpG.cov.bed 02.methylation/dianli_hap1.CHG.cov.bed 02.methylation/dianli_hap1.CHH.cov.bed 30000 .

Output:
  03.windows/<sample_name>.CG.<window_size>.bed
  03.windows/<sample_name>.CHG.<window_size>.bed
  03.windows/<sample_name>.CHH.<window_size>.bed
EOF
}

if [ "$#" -lt 5 ] || [ "$#" -gt 7 ]; then
  usage
  exit 1
fi

SAMPLE="$1"
GENOME="$2"
CPG="$3"
CHG="$4"
CHH="$5"
WIN="${6:-30000}"
PROJECT_DIR="${7:-.}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for f in "$GENOME" "$CPG" "$CHG" "$CHH"; do
  if [ ! -s "$f" ]; then
    echo "ERROR: missing file: $f" >&2
    exit 1
  fi
done

mkdir -p "${PROJECT_DIR}/03.windows" "${PROJECT_DIR}/logs"

python3 "${SCRIPT_DIR}/bin/window_methylation.py" \
  --sample "$SAMPLE" \
  --genome "$GENOME" \
  --cpg "$CPG" \
  --chg "$CHG" \
  --chh "$CHH" \
  --window "$WIN" \
  --outdir "${PROJECT_DIR}/03.windows" \
  > "${PROJECT_DIR}/logs/03.${SAMPLE}.window_methylation.log" 2>&1

echo "[03] done"
