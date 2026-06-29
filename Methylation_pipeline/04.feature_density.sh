#!/usr/bin/env bash
set -euo pipefail

usage() {
cat <<'EOF'
Usage:
  bash 04.feature_density.sh <sample_name> <genome.fa> <gene.gff3> <centromere_range.txt> <SatDNA.bed> <Copia.bed|none> <Gypsy.bed|none> [window_size] [project_dir]

Example:
  bash 04.feature_density.sh dianli_hap1 dianli/hap1/Pyrps_dianli_hap1.fa dianli/hap1/Pyrps_dianli_hap1.gff3 dianli/hap1/0.Pyrps_dianli.hap1.FINAL.fixTelo_range.txt dianli/hap1/Pyrps_dianli.hap1.SatDNA.bed dianli/hap1/Copia.bed dianli/hap1/Gypsy.bed 30000 .

If Copia or Gypsy is not available, use:
  none

Output:
  04.features/<sample_name>.gene.bed
  04.features/<sample_name>.feature_density.<window_size>.tsv
EOF
}

if [ "$#" -lt 7 ] || [ "$#" -gt 9 ]; then
  usage
  exit 1
fi

SAMPLE="$1"
GENOME="$2"
GFF3="$3"
CEN="$4"
SAT="$5"
COPIA="$6"
GYPSY="$7"
WIN="${8:-30000}"
PROJECT_DIR="${9:-.}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for f in "$GENOME" "$GFF3" "$CEN" "$SAT"; do
  if [ ! -s "$f" ]; then
    echo "ERROR: missing file: $f" >&2
    exit 1
  fi
done

if [ "$COPIA" != "none" ] && [ ! -s "$COPIA" ]; then
  echo "ERROR: missing Copia file: $COPIA" >&2
  exit 1
fi
if [ "$GYPSY" != "none" ] && [ ! -s "$GYPSY" ]; then
  echo "ERROR: missing Gypsy file: $GYPSY" >&2
  exit 1
fi

mkdir -p "${PROJECT_DIR}/04.features" "${PROJECT_DIR}/logs"

python3 "${SCRIPT_DIR}/bin/feature_density.py" \
  --sample "$SAMPLE" \
  --genome "$GENOME" \
  --gff3 "$GFF3" \
  --centromere "$CEN" \
  --satdna "$SAT" \
  --copia "$COPIA" \
  --gypsy "$GYPSY" \
  --window "$WIN" \
  --outdir "${PROJECT_DIR}/04.features" \
  > "${PROJECT_DIR}/logs/04.${SAMPLE}.feature_density.log" 2>&1

echo "[04] done"
