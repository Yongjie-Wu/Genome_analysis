#!/usr/bin/env bash
set -euo pipefail

usage() {
cat <<'EOF'
Usage:
  bash prepare_circos_inputs.sh <sample_name> <genome.fa.fai> <centromere_range.txt> <Gypsy.bed> <Copia.bed> <gene.gff3|Gene.bed> <SatDNA.bed> <CpG.cov.bed> <CHG.cov.bed> <CHH.cov.bed> <output_dir> [options]

Required input:
  sample_name
      Prefix for output files, for example Pyrbe_duli_hap1.

  genome.fa.fai
      FASTA index file. First two columns must be chromosome name and length.

  centromere_range.txt
      Three-column centromere file:
        chr  start  end

  Gypsy.bed / Copia.bed / SatDNA.bed
      BED files used to calculate repeat density in each window.

  gene.gff3 or Gene.bed
      Gene annotation used to calculate gene density in each window.
      GFF3 gene records are parsed by feature type "gene".

  CpG.cov.bed / CHG.cov.bed / CHH.cov.bed
      hifimeth methylation files used to calculate CG/CHG/CHH methylation levels.

  output_dir
      Directory for generated circos input tables.

Options:
  --window <N>
      Window size. Default: 30000

  --chr-regex <regex>
      Chromosomes to keep. Default: ^Chr[0-9]+A$

  --rpkm <RPKM.tsv>
      Optional expression table. Required columns:
        ID, leaf, fruit
      Output:
        <sample_name>.leaf_fruit_expression.<window>bp.tsv

  --atac-bam <merged.rmdup.bam>
      Optional ATAC BAM file. The script runs bedtools coverage on windows.
      Output:
        <sample_name>.ATAC.<window>bp.tsv

  --threads <N>
      Threads for samtools. Default: 4

Output:
  <sample_name>.windows.<window>bp.bed
  <sample_name>.Gypsy.<window>bp.txt
  <sample_name>.Copia.<window>bp.txt
  <sample_name>.Gene.<window>bp.txt
  <sample_name>.SatDNA.<window>bp.txt
  <sample_name>.CG.<window>bp.txt
  <sample_name>.CHG.<window>bp.txt
  <sample_name>.CHH.<window>bp.txt
  <sample_name>.ATAC.<window>bp.tsv                         if --atac-bam is used
  <sample_name>.leaf_fruit_expression.<window>bp.tsv        if --rpkm is used

Example:
  bash prepare_circos_inputs.sh Pyrbe_duli_hap1 \
    Pyrbe_duli_hap1.fa.fai \
    Pyrbe_duli.hap1_centromere_range.txt \
    Gypsy.bed \
    Copia.bed \
    Pyrbe_duli_hap1.gff3 \
    Pyrbe_duli.hap1.SatDNA.bed \
    hap1.CpG.cov.bed \
    hap1.CHG.cov.bed \
    hap1.CHH.cov.bed \
    circos_input \
    --rpkm Pyrbe_duli_hap1.RPKM \
    --atac-bam /storage3/wuyongjie/Project/centromere/atct/duli/03.bam_filter/duli.hap1.merge.rmdup.bam
EOF
}

if [ "$#" -lt 11 ]; then
  usage
  exit 1
fi

SAMPLE="$1"
FAI="$2"
CENTROMERE="$3"
GYPSY="$4"
COPIA="$5"
GENE="$6"
SATDNA="$7"
CPG="$8"
CHG="$9"
CHH="${10}"
OUTDIR="${11}"
shift 11

WINDOW=30000
CHR_REGEX='^Chr[0-9]+A$'
RPKM=""
ATAC_BAM=""
THREADS=4

while [ "$#" -gt 0 ]; do
  case "$1" in
    --window)
      WINDOW="$2"; shift 2 ;;
    --chr-regex)
      CHR_REGEX="$2"; shift 2 ;;
    --rpkm)
      RPKM="$2"; shift 2 ;;
    --atac-bam)
      ATAC_BAM="$2"; shift 2 ;;
    --threads)
      THREADS="$2"; shift 2 ;;
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
[ -z "$ATAC_BAM" ] || [ -s "$ATAC_BAM" ] || { echo "ERROR: missing ATAC BAM: $ATAC_BAM" >&2; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$OUTDIR"

CONDA_ENV="${CIRCOS_CONDA_ENV:-base}"
if [ -s "${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
  set +u
  source "${HOME}/miniconda3/etc/profile.d/conda.sh"
  conda activate "${CONDA_ENV}"
  set -u
fi

PYTHON_BIN="$(command -v python || command -v python3)"

echo "[prepare] make genome feature and methylation window tracks"
PY_ARGS=(
  --sample "$SAMPLE"
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
  --chr-regex "$CHR_REGEX"
)
if [ -n "$RPKM" ]; then
  PY_ARGS+=(--rpkm "$RPKM")
fi

"$PYTHON_BIN" "${SCRIPT_DIR}/scripts/prepare_circos_inputs.py" \
  "${PY_ARGS[@]}"

if [ -n "$ATAC_BAM" ]; then
  command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found" >&2; exit 1; }
  command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found" >&2; exit 1; }

  WINDOWS="${OUTDIR}/${SAMPLE}.windows.${WINDOW}bp.bed"
  ATAC_OUT="${OUTDIR}/${SAMPLE}.ATAC.${WINDOW}bp.tsv"

  echo "[prepare] count ATAC reads by window: $ATAC_BAM"
  MAPPED=$(samtools view -@ "$THREADS" -c -F 4 "$ATAC_BAM")
  {
    echo -e "chr\tstart\tend\twindow_size\tread_count\treads_per_kb\tCPM_per_chr_mapped"
    bedtools coverage -a "$WINDOWS" -b "$ATAC_BAM" -counts \
      | awk -v OFS="\t" -v mapped="$MAPPED" '{w=$3-$2; count=$4; cpm=(mapped>0 ? count*1000000/mapped : 0); print $1,$2,$3,w,count,count/(w/1000),cpm}'
  } > "$ATAC_OUT"
  echo "[prepare] wrote $ATAC_OUT"
fi

echo "[prepare] done: $OUTDIR"
