#!/usr/bin/env bash
set -euo pipefail
usage(){ cat <<'USAGE'
Usage:
  bash run_all.sh <project_name> <samples.tsv> <genome.fa> <gene.gff3> <genome_size> [project_dir] [options]

samples.tsv format:
  sample<TAB>R1<TAB>R2
  CRRxxxx<TAB>CRRxxxx_r1.fq.gz<TAB>CRRxxxx_r2.fq.gz

Options:
  --skip-clean
      Use existing 00.clean/<sample>.clean.r1/r2.fq.gz and skip fastp.

  --skip-step08
      Skip ATAC window feature analysis.

  --step08-extra <centromere_range.txt> <SatDNA.bed> <RPKM.tsv>
      Run step 08 ATAC/SatDNA/expression window analysis.

Example, run 00-07 and 09 summary:
  bash run_all.sh duli_hap1 samples.tsv genome/Pyrbe_hap1.fa genome/Pyrbe_hap1.gff3 5e8 .

Example, run 00-09:
  bash run_all.sh duli_hap1 samples.tsv genome/Pyrbe_hap1.fa genome/Pyrbe_hap1.gff3 5e8 . \
    --step08-extra hap1_centromere_range.txt hap1.SatDNA.bed hap1.RPKM

Example, start from existing clean reads:
  bash run_all.sh duli_hap1 samples.clean.tsv genome/Pyrbe_hap1.fa genome/Pyrbe_hap1.gff3 5e8 . --skip-clean
USAGE
}
[ "$#" -ge 5 ] || { usage; exit 1; }
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
project="$1"; samples_tsv="$2"; genome_fa="$3"; gene_gff="$4"; genome_size="$5"; project_dir="${6:-.}"
shift 5
if [ "$#" -gt 0 ] && [[ "$1" != --* ]]; then project_dir="$1"; shift; fi
skip_clean=0; skip_step08=0; step08_extra=()
while [ "$#" -gt 0 ]; do
  case "$1" in
    --skip-clean) skip_clean=1; shift ;;
    --skip-step08) skip_step08=1; shift ;;
    --step08-extra) step08_extra=("$2" "$3" "$4"); shift 4 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage; exit 1 ;;
  esac
done
[ -s "$samples_tsv" ] || { echo "ERROR: missing samples.tsv: $samples_tsv" >&2; exit 1; }
cd "$project_dir"
mkdir -p logs
if [ "$skip_clean" -eq 0 ]; then
  awk 'BEGIN{FS="\t"} NR>1 && $1!="" {print $1"\t"$2"\t"$3}' "$samples_tsv" | while IFS=$'\t' read -r sample r1 r2; do
    bash "$SCRIPT_DIR/00.clean.sh" "$sample" "$r1" "$r2" "$PWD"
  done
fi
bash "$SCRIPT_DIR/01.index.sh" "$genome_fa" "$PWD"
awk 'BEGIN{FS="\t"} NR>1 && $1!="" {print $1}' "$samples_tsv" | while read -r sample; do
  bash "$SCRIPT_DIR/02.map.sh" "$sample" "00.clean/${sample}.clean.r1.fq.gz" "00.clean/${sample}.clean.r2.fq.gz" "$genome_fa" "$PWD"
done
mapfile -t sample_names < <(awk 'BEGIN{FS="\t"} NR>1 && $1!="" {print $1}' "$samples_tsv")
bash "$SCRIPT_DIR/03.filter_markdup_merge.sh" "$project" "$PWD" "${sample_names[@]}"
bash "$SCRIPT_DIR/04.macs2_callpeak.sh" "$project" "$genome_size" "$PWD"
bash "$SCRIPT_DIR/05.make_bigwig.sh" "$project" "$PWD"
bash "$SCRIPT_DIR/06.qc_peak_signal.sh" "$project" "$PWD" "${sample_names[@]}"
bash "$SCRIPT_DIR/07.peak_gene_annotation.sh" "$project" "$gene_gff" "$PWD"
if [ "$skip_step08" -eq 0 ] && [ "${#step08_extra[@]}" -eq 3 ]; then
  bash "$SCRIPT_DIR/08.atac_window_features.sh" "$project" "$genome_fa" "$gene_gff" "${step08_extra[0]}" "${step08_extra[1]}" "${step08_extra[2]}" "$PWD"
fi
bash "$SCRIPT_DIR/09.summary.sh" "$project" "$PWD"
