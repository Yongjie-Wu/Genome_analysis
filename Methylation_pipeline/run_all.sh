#!/usr/bin/env bash
set -euo pipefail

usage() {
cat <<'EOF'
Usage:
  bash run_all.sh <sample_name> <mod.bam> <genome.fa> <gene.gff3> <centromere_range.txt> <SatDNA.bed> <Copia.bed|none> <Gypsy.bed|none> <genome_ref_prefix> [project_dir] [options]

Options:
  --threads <N>          default: 48
  --window <N>           default: 30000
  --plot-chr <chr_id>    draw one chromosome track after window analysis
  --skip-align           use existing 01.align/<sample_name>.pbmm2.bam
  --skip-plot            do not draw chromosome plot

Example, full single-haplotype methylation workflow:
  bash run_all.sh dianli_hap1 dianli/mod.bam dianli/hap1/Pyrps_dianli_hap1.fa dianli/hap1/Pyrps_dianli_hap1.gff3 dianli/hap1/0.Pyrps_dianli.hap1.FINAL.fixTelo_range.txt dianli/hap1/Pyrps_dianli.hap1.SatDNA.bed dianli/hap1/Copia.bed dianli/hap1/Gypsy.bed dianli_hap1 . --threads 48 --window 30000 --plot-chr Chr02A

Example, if Copia/Gypsy files are not available:
  bash run_all.sh liuyeli_hap1 liuyeli/mod.bam liuyeli/hap1/Pyrsa_liuyeli_hap1.fa liuyeli/hap1/Pyrsa_liuyeli_hap1.gff3 liuyeli/hap1/centromere_range.txt liuyeli/hap1/SatDNA.bed none none liuyeli_hap1 .
EOF
}

if [ "$#" -lt 9 ]; then
  usage
  exit 1
fi

SAMPLE="$1"
MOD_BAM="$2"
GENOME="$3"
GFF3="$4"
CEN="$5"
SAT="$6"
COPIA="$7"
GYPSY="$8"
REF_PREFIX="$9"
PROJECT_DIR="${10:-.}"
shift $(( $# >= 10 ? 10 : 9 ))

THREADS=48
WIN=30000
PLOT_CHR=""
SKIP_ALIGN=0
SKIP_PLOT=0

while [ "$#" -gt 0 ]; do
  case "$1" in
    --threads)
      THREADS="$2"; shift 2 ;;
    --window)
      WIN="$2"; shift 2 ;;
    --plot-chr)
      PLOT_CHR="$2"; shift 2 ;;
    --skip-align)
      SKIP_ALIGN=1; shift ;;
    --skip-plot)
      SKIP_PLOT=1; shift ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      echo "ERROR: unknown option: $1" >&2
      usage
      exit 1 ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "${PROJECT_DIR}/logs"

echo "[run_all] sample=${SAMPLE}"
echo "[run_all] project_dir=${PROJECT_DIR}"
echo "[run_all] window=${WIN}"

bash "${SCRIPT_DIR}/00.index_reference.sh" "$GENOME" "$REF_PREFIX" "$PROJECT_DIR"

MMI="${PROJECT_DIR}/00.index/${REF_PREFIX}.mmi"
BAM="${PROJECT_DIR}/01.align/${SAMPLE}.pbmm2.bam"

if [ "$SKIP_ALIGN" -eq 0 ]; then
  bash "${SCRIPT_DIR}/01.align_modbam.sh" "$SAMPLE" "$MOD_BAM" "$GENOME" "$MMI" "$PROJECT_DIR" "$THREADS"
else
  if [ ! -s "$BAM" ]; then
    echo "ERROR: --skip-align was used, but BAM does not exist: $BAM" >&2
    exit 1
  fi
fi

bash "${SCRIPT_DIR}/02.hifimeth_pileup.sh" "$SAMPLE" "$GENOME" "$BAM" "$PROJECT_DIR"

bash "${SCRIPT_DIR}/03.window_methylation.sh" \
  "$SAMPLE" "$GENOME" \
  "${PROJECT_DIR}/02.methylation/${SAMPLE}.CpG.cov.bed" \
  "${PROJECT_DIR}/02.methylation/${SAMPLE}.CHG.cov.bed" \
  "${PROJECT_DIR}/02.methylation/${SAMPLE}.CHH.cov.bed" \
  "$WIN" "$PROJECT_DIR"

bash "${SCRIPT_DIR}/04.feature_density.sh" \
  "$SAMPLE" "$GENOME" "$GFF3" "$CEN" "$SAT" "$COPIA" "$GYPSY" "$WIN" "$PROJECT_DIR"

if [ "$SKIP_PLOT" -eq 0 ] && [ -n "$PLOT_CHR" ]; then
  bash "${SCRIPT_DIR}/05.plot_chr_tracks.sh" "$SAMPLE" "$PLOT_CHR" "$WIN" "$PROJECT_DIR"
fi

bash "${SCRIPT_DIR}/06.summary.sh" "$SAMPLE" "$PROJECT_DIR"

echo "[run_all] done"
