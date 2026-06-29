#!/usr/bin/env bash
set -euo pipefail

usage() {
cat <<'EOF'
Usage:
  bash run_genome_feature_circos.sh <genome.fa.fai> <centromere_range.txt> <Gypsy_30kb.txt> <Copia_30kb.txt> <Gene_30kb.txt> <SatDNA_30kb.txt> <CG_30kb.txt> <CHG_30kb.txt> <CHH_30kb.txt> <output_prefix> [options]

Required input:
  genome.fa.fai
      FASTA index file. The first two columns must be chromosome name and length.

  centromere_range.txt
      Three-column file used to mark centromere positions on the chromosome ring:
        chr  start  end

  output_prefix
      Output file prefix. The script writes .png, .pdf and .svg.

  Gypsy_30kb.txt / Copia_30kb.txt / Gene_30kb.txt / SatDNA_30kb.txt
      Four-column genome-wide window files:
        chr  start  end  value

  CG_30kb.txt / CHG_30kb.txt / CHH_30kb.txt
      Four-column genome-wide methylation window files:
        chr  start  end  methylation_level
      CG and CHG are plotted from 0 to 1.
      CHH is plotted from 0 to 0.30.

Options:
  --expression <expression_30kb.tsv>
      Optional genome-wide expression file. Required columns:
        chr, start, end, leaf, fruit

  --atac <ATAC_window.tsv>
      Optional ATAC window file. Required columns:
        chr, start, end, CPM_per_chr_mapped

Tracks drawn:
  Outer chromosome ring
  Centromere ranges, drawn as black blocks on the chromosome ring
  Copia density
  Gypsy density
  SatDNA density
  ATAC signal, if --atac is provided
  Gene density
  Leaf expression, if expression files exist
  Fruit expression, if expression files exist
  CG methylation
  CHG methylation
  CHH methylation

  --chr-regex <regex>
      Chromosomes to plot. Default: ^Chr[0-9]+A$

  --chr-suffix <suffix>
      Remove this suffix in chromosome labels. Default: A

  --gap-bp <N>
      Gap size between chromosomes in pseudo genomic coordinates. Default: 1500000

  --plot-span-deg <N>
      Angular span used by chromosomes. Default: 350

Environment:
  CIRCOS_CONDA_ENV
      Conda environment used for Python packages. Default: base

Example:
  bash run_genome_feature_circos.sh \
    /storage3/wuyongjie/Project/centromere/jiajihua/Pyrbe_duli/zuohaodehap1/Pyrbe_duli_hap1.fa.fai \
    /storage3/wuyongjie/Project/centromere/jiajihua/Pyrbe_duli/zuohaodehap1/Pyrbe_duli.hap1_centromere_range.txt \
    circos_input/Pyrbe_duli_hap1.Gypsy.30000bp.txt \
    circos_input/Pyrbe_duli_hap1.Copia.30000bp.txt \
    circos_input/Pyrbe_duli_hap1.Gene.30000bp.txt \
    circos_input/Pyrbe_duli_hap1.SatDNA.30000bp.txt \
    circos_input/Pyrbe_duli_hap1.CG.30000bp.txt \
    circos_input/Pyrbe_duli_hap1.CHG.30000bp.txt \
    circos_input/Pyrbe_duli_hap1.CHH.30000bp.txt \
    /storage3/wuyongjie/Project/centromere/jiajihua/Pyrbe_duli/zuohaodehap1/circos_plot/output/Pyrbe_duli_hap1_17chr_multi_track_circos_fast \
    --atac circos_input/Pyrbe_duli_hap1.ATAC.30000bp.tsv \
    --expression circos_input/Pyrbe_duli_hap1.leaf_fruit_expression.30000bp.tsv
EOF
}

if [ "$#" -lt 10 ]; then
  usage
  exit 1
fi

FAI="$1"
CENTROMERE="$2"
GYPSY="$3"
COPIA="$4"
GENE="$5"
SATDNA="$6"
CG="$7"
CHG="$8"
CHH="$9"
OUT_PREFIX="${10}"
shift 10

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="${SCRIPT_DIR}/scripts/plot_genome_feature_circos.py"

if [ ! -s "$PY_SCRIPT" ]; then
  echo "ERROR: cannot find Python script: $PY_SCRIPT" >&2
  exit 1
fi

for f in "$FAI" "$CENTROMERE" "$GYPSY" "$COPIA" "$GENE" "$SATDNA" "$CG" "$CHG" "$CHH"; do
  if [ ! -s "$f" ]; then
    echo "ERROR: input file not found or empty: $f" >&2
    exit 1
  fi
done

mkdir -p "$(dirname "$OUT_PREFIX")"

CONDA_ENV="${CIRCOS_CONDA_ENV:-base}"
if [ -s "${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
  # The original project uses conda base for matplotlib/pandas/numpy.
  set +u
  source "${HOME}/miniconda3/etc/profile.d/conda.sh"
  conda activate "${CONDA_ENV}"
  set -u
fi

if command -v python >/dev/null 2>&1; then
  PYTHON_BIN="$(command -v python)"
elif command -v python3 >/dev/null 2>&1; then
  PYTHON_BIN="$(command -v python3)"
else
  echo "ERROR: python not found" >&2
  exit 1
fi

echo "[circos] python: ${PYTHON_BIN}"
"${PYTHON_BIN}" "${PY_SCRIPT}" \
  --fai "$FAI" \
  --centromere "$CENTROMERE" \
  --gypsy "$GYPSY" \
  --copia "$COPIA" \
  --gene "$GENE" \
  --satdna "$SATDNA" \
  --cg "$CG" \
  --chg "$CHG" \
  --chh "$CHH" \
  --out-prefix "$OUT_PREFIX" \
  "$@"

echo "[circos] done"
