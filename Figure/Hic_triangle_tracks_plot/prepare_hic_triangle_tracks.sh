#!/usr/bin/env bash
set -euo pipefail

usage() {
cat <<'EOF'
Usage:
  bash prepare_hic_triangle_tracks.sh <sample_name> <chromosome> <genome.fa.fai> <centromere_range.txt> <Gypsy.bed> <Copia.bed> <gene.gff3|Gene.bed> <SatDNA.bed> <CpG.cov.bed> <CHG.cov.bed> <CHH.cov.bed> <output_dir> [options]

Required input:
  sample_name
      Prefix for output files, for example Pyrbe_duli_hap1.

  chromosome
      Chromosome to plot, for example Chr02A.

  genome.fa.fai
      FASTA index file. First two columns must be chromosome name and length.

  centromere_range.txt
      Three-column centromere file:
        chr  start  end

  Gypsy.bed / Copia.bed / SatDNA.bed
      BED files used to calculate density in each window.

  gene.gff3 or Gene.bed
      Gene annotation used to calculate gene count in each window.

  CpG.cov.bed / CHG.cov.bed / CHH.cov.bed
      hifimeth methylation files used to calculate CG/CHG/CHH methylation levels.

  output_dir
      Directory for generated chromosome track tables.

Options:
  --window <N>
      Window size for tracks. Default: 30000

  --rpkm <RPKM.tsv>
      Optional expression table. Required columns:
        ID, leaf, fruit

  --gypsy-ltr <Gypsy_LTR.bed>
      Optional Gypsy LTR position BED, drawn as vertical ticks.

  --gypsy-int <Gypsy_INT.bed>
      Optional Gypsy INT position BED, drawn as vertical ticks.

Output:
  <sample>.<chr>.windows.<window>bp.bed
  <sample>.<chr>.Gypsy.<window>bp.txt
  <sample>.<chr>.Copia.<window>bp.txt
  <sample>.<chr>.Gene.<window>bp.txt
  <sample>.<chr>.SatDNA.<window>bp.txt
  <sample>.<chr>.CG.<window>bp.txt
  <sample>.<chr>.CHG.<window>bp.txt
  <sample>.<chr>.CHH.<window>bp.txt
  <sample>.<chr>.leaf_fruit_expression.<window>bp.tsv
  <sample>.<chr>.Gypsy_LTR.bed
  <sample>.<chr>.Gypsy_INT.bed

Example:
  bash prepare_hic_triangle_tracks.sh Pyrbe_duli_hap1 Chr02A \
    Pyrbe_duli_hap1.fa.fai \
    Pyrbe_duli.hap1_centromere_range.txt \
    Gypsy.bed \
    Copia.bed \
    Pyrbe_duli_hap1.gff3 \
    Pyrbe_duli.hap1.SatDNA.bed \
    hap1.CpG.cov.bed \
    hap1.CHG.cov.bed \
    hap1.CHH.cov.bed \
    hic_tracks_input \
    --rpkm Pyrbe_duli_hap1.RPKM \
    --gypsy-ltr Chr02A/Chr02A_Gypsy_LTR.bed \
    --gypsy-int Chr02A/Chr02A_Gypsy_INT.bed
EOF
}

if [ "$#" -lt 12 ]; then
  usage
  exit 1
fi

SAMPLE="$1"
CHR="$2"
FAI="$3"
CENTROMERE="$4"
GYPSY="$5"
COPIA="$6"
GENE="$7"
SATDNA="$8"
CPG="$9"
CHG="${10}"
CHH="${11}"
OUTDIR="${12}"
shift 12

WINDOW=30000
RPKM=""
GYPSY_LTR=""
GYPSY_INT=""

while [ "$#" -gt 0 ]; do
  case "$1" in
    --window)
      WINDOW="$2"; shift 2 ;;
    --rpkm)
      RPKM="$2"; shift 2 ;;
    --gypsy-ltr)
      GYPSY_LTR="$2"; shift 2 ;;
    --gypsy-int)
      GYPSY_INT="$2"; shift 2 ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      echo "ERROR: unknown option: $1" >&2
      usage
      exit 1 ;;
  esac
done

for f in "$FAI" "$CENTROMERE" "$GYPSY" "$COPIA" "$GENE" "$SATDNA" "$CPG" "$CHG" "$CHH"; do
  [ -s "$f" ] || { echo "ERROR: missing or empty file: $f" >&2; exit 1; }
done
[ -z "$RPKM" ] || [ -s "$RPKM" ] || { echo "ERROR: missing RPKM file: $RPKM" >&2; exit 1; }
[ -z "$GYPSY_LTR" ] || [ -s "$GYPSY_LTR" ] || { echo "ERROR: missing Gypsy LTR file: $GYPSY_LTR" >&2; exit 1; }
[ -z "$GYPSY_INT" ] || [ -s "$GYPSY_INT" ] || { echo "ERROR: missing Gypsy INT file: $GYPSY_INT" >&2; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$OUTDIR"

CONDA_ENV="${HIC_TRACK_CONDA_ENV:-base}"
if [ -s "${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
  set +u
  source "${HOME}/miniconda3/etc/profile.d/conda.sh"
  conda activate "${CONDA_ENV}"
  set -u
fi

PYTHON_BIN="$(command -v python || command -v python3)"
PY_ARGS=(
  --sample "$SAMPLE"
  --chromosome "$CHR"
  --fai "$FAI"
  --centromere "$CENTROMERE"
  --gypsy "$GYPSY"
  --copia "$COPIA"
  --gene "$GENE"
  --satdna "$SATDNA"
  --cpg "$CPG"
  --chg "$CHG"
  --chh "$CHH"
  --outdir "$OUTDIR"
  --window "$WINDOW"
)
[ -z "$RPKM" ] || PY_ARGS+=(--rpkm "$RPKM")
[ -z "$GYPSY_LTR" ] || PY_ARGS+=(--gypsy-ltr "$GYPSY_LTR")
[ -z "$GYPSY_INT" ] || PY_ARGS+=(--gypsy-int "$GYPSY_INT")

"$PYTHON_BIN" "${SCRIPT_DIR}/scripts/prepare_hic_triangle_tracks.py" "${PY_ARGS[@]}"

echo "[prepare] done: $OUTDIR"
