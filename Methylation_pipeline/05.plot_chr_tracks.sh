#!/usr/bin/env bash
set -euo pipefail

usage() {
cat <<'EOF'
Usage:
  bash 05.plot_chr_tracks.sh <sample_name> <chr_id> [window_size] [project_dir]

Example:
  bash 05.plot_chr_tracks.sh dianli_hap1 Chr02A 30000 .

Input:
  03.windows/<sample_name>.CG.<window_size>.bed
  03.windows/<sample_name>.CHG.<window_size>.bed
  03.windows/<sample_name>.CHH.<window_size>.bed
  04.features/<sample_name>.feature_density.<window_size>.tsv

Output:
  05.plots/<sample_name>.<chr_id>.methylation_tracks.pdf
  05.plots/<sample_name>.<chr_id>.methylation_tracks.png
EOF
}

if [ "$#" -lt 2 ] || [ "$#" -gt 4 ]; then
  usage
  exit 1
fi

SAMPLE="$1"
CHR="$2"
WIN="${3:-30000}"
PROJECT_DIR="${4:-.}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "${PROJECT_DIR}/05.plots" "${PROJECT_DIR}/logs"

Rscript "${SCRIPT_DIR}/bin/plot_chr_tracks.R" \
  "$SAMPLE" "$CHR" "$WIN" "$PROJECT_DIR" \
  > "${PROJECT_DIR}/logs/05.${SAMPLE}.${CHR}.plot_chr_tracks.log" 2>&1

echo "[05] done"
