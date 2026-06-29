#!/usr/bin/env bash
set -euo pipefail

usage() {
cat <<'EOF'
Usage:
  bash plot_hic_triangle_tracks.sh <cool_file> <chromosome> <centromere_range.txt> <Copia_track.txt> <Gypsy_track.txt> <SatDNA_track.txt> <Gene_track.txt> <LeafFruit_expression.tsv> <CG_track.txt> <CHG_track.txt> <CHH_track.txt> <output_prefix> [options]

Required input:
  cool_file
      Hi-C .cool file, or .mcool URI such as sample.mcool::/resolutions/10000.

  chromosome
      Chromosome to plot, for example Chr02A.

  centromere_range.txt
      Three-column centromere file:
        chr  start  end

  Copia/Gypsy/SatDNA/Gene/CG/CHG/CHH track files
      Four-column window files:
        chr  start  end  value

  LeafFruit_expression.tsv
      Expression window file with columns:
        chr  start  end  leaf  fruit
      If no expression track is available, use "none".

  output_prefix
      Output file prefix. The script writes .png and .pdf.

Options:
  --gypsy-ltr <Gypsy_LTR.bed>
      Optional tick track near the bottom.

  --gypsy-int <Gypsy_INT.bed>
      Optional tick track near the bottom.

  --dpi <N>
      Output PNG dpi. Default: 600

  --no-pdf
      Do not save PDF.

  --hic-bin-step <N>
      Downsample Hi-C matrix by keeping every Nth bin for plotting. Default: 1.
      Use 2-4 for quick previews.

  --font-file <Arial.ttf>
      Optional Arial font file. Use this when the server does not have Arial installed.

  --max-distance <bp>
      Maximum Hi-C triangular height. Default: chromosome_length/2.

  --hic-vmax-percentile <N>
      Percentile for Hi-C color max on log1p scale. Default: 99.5

Environment:
  HIC_TRACK_CONDA_ENV
      Conda environment containing cooler/matplotlib/pandas/numpy. Default: base

Example:
  bash plot_hic_triangle_tracks.sh \
    /storage3/wuyongjie/Project/centromere/hic/Pyrbe_duli/hap1/Pyrbe_duli_hap1.10k.cool \
    Chr02A \
    Pyrbe_duli.hap1_centromere_range.txt \
    hic_tracks_input/Pyrbe_duli_hap1.Chr02A.Copia.30000bp.txt \
    hic_tracks_input/Pyrbe_duli_hap1.Chr02A.Gypsy.30000bp.txt \
    hic_tracks_input/Pyrbe_duli_hap1.Chr02A.SatDNA.30000bp.txt \
    hic_tracks_input/Pyrbe_duli_hap1.Chr02A.Gene.30000bp.txt \
    hic_tracks_input/Pyrbe_duli_hap1.Chr02A.leaf_fruit_expression.30000bp.tsv \
    hic_tracks_input/Pyrbe_duli_hap1.Chr02A.CG.30000bp.txt \
    hic_tracks_input/Pyrbe_duli_hap1.Chr02A.CHG.30000bp.txt \
    hic_tracks_input/Pyrbe_duli_hap1.Chr02A.CHH.30000bp.txt \
    output/Chr02A_hic_triangle_plus_multi_tracks \
    --gypsy-ltr hic_tracks_input/Pyrbe_duli_hap1.Chr02A.Gypsy_LTR.bed \
    --gypsy-int hic_tracks_input/Pyrbe_duli_hap1.Chr02A.Gypsy_INT.bed
EOF
}

if [ "$#" -lt 12 ]; then
  usage
  exit 1
fi

COOL="$1"
CHR="$2"
CENTROMERE="$3"
COPIA="$4"
GYPSY="$5"
SATDNA="$6"
GENE="$7"
EXPRESSION="$8"
CG="$9"
CHG="${10}"
CHH="${11}"
OUT_PREFIX="${12}"
shift 12

for f in "$CENTROMERE" "$COPIA" "$GYPSY" "$SATDNA" "$GENE" "$CG" "$CHG" "$CHH"; do
  [ -s "$f" ] || { echo "ERROR: missing or empty file: $f" >&2; exit 1; }
done
[ "$EXPRESSION" = "none" ] || [ -s "$EXPRESSION" ] || { echo "ERROR: missing expression file: $EXPRESSION" >&2; exit 1; }

COOL_PATH="${COOL%%::*}"
[ -s "$COOL_PATH" ] || { echo "ERROR: missing cool file: $COOL_PATH" >&2; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$(dirname "$OUT_PREFIX")"

CONDA_ENV="${HIC_TRACK_CONDA_ENV:-base}"
if [ -s "${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
  set +u
  source "${HOME}/miniconda3/etc/profile.d/conda.sh"
  conda activate "${CONDA_ENV}"
  set -u
fi

PYTHON_BIN="$(command -v python || command -v python3)"
"$PYTHON_BIN" "${SCRIPT_DIR}/scripts/plot_hic_triangle_tracks.py" \
  --cool "$COOL" \
  --chromosome "$CHR" \
  --centromere "$CENTROMERE" \
  --copia "$COPIA" \
  --gypsy "$GYPSY" \
  --satdna "$SATDNA" \
  --gene "$GENE" \
  --expression "$EXPRESSION" \
  --cg "$CG" \
  --chg "$CHG" \
  --chh "$CHH" \
  --out-prefix "$OUT_PREFIX" \
  "$@"

echo "[plot] done"
