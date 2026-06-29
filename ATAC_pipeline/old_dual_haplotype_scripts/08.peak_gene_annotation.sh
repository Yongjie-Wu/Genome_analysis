#!/usr/bin/env bash
set -euo pipefail
usage(){ cat <<'USAGE'
Usage:
  bash 08.peak_gene_annotation.sh <project_name> <hap1.gff3> <hap2.gff3> [project_dir]

Example:
  bash 08.peak_gene_annotation.sh duli hap1/Pyrbe_duli_hap1.gff3 hap2/Pyrbe_duli_hap2.gff3 .

Input:
  04.peak/hap1/<project_name>.hap1.ATAC_peaks.narrowPeak
  04.peak/hap2/<project_name>.hap2.ATAC_peaks.narrowPeak
Output:
  07.annotation/<project_name>.hap*.ATAC_peaks.overlap_genes.tsv
  07.annotation/<project_name>.hap*.ATAC_peaks.closest_gene.tsv
USAGE
}
[ "$#" -ge 3 ] || { usage; exit 1; }
project="$1"; hap1_gff="$2"; hap2_gff="$3"; project_dir="${4:-.}"
cd "$project_dir"; mkdir -p 07.annotation logs
LOG="logs/08.peak_gene_annotation.$(date +%Y%m%d_%H%M%S).log"; exec > >(tee -a "$LOG") 2>&1
log(){ echo "[$(date '+%F %T')] $*"; }; need_file(){ [ -s "$1" ] || { echo "ERROR: missing or empty file: $1" >&2; exit 1; }; }
source "${HOME}/miniconda3/etc/profile.d/conda.sh"; conda activate "${WORK_ENV:-work}"
command -v bedtools
for hap in hap1 hap2; do
  gff="$hap1_gff"; [ "$hap" = hap2 ] && gff="$hap2_gff"
  peaks="04.peak/${hap}/${project}.${hap}.ATAC_peaks.narrowPeak"; gene_bed="07.annotation/${project}.${hap}.genes.bed"
  need_file "$gff"; need_file "$peaks"
  awk 'BEGIN{OFS="\t"} $0 !~ /^#/ && ($3=="gene" || $3=="mRNA") {id="."; if (match($9,/ID=([^;]+)/,a)) id=a[1]; else if (match($9,/Name=([^;]+)/,b)) id=b[1]; print $1,$4-1,$5,id,".",$7}' "$gff" | sort -k1,1 -k2,2n > "$gene_bed"
  bedtools intersect -a "$peaks" -b "$gene_bed" -wa -wb > "07.annotation/${project}.${hap}.ATAC_peaks.overlap_genes.tsv"
  bedtools closest -a "$peaks" -b "$gene_bed" -d > "07.annotation/${project}.${hap}.ATAC_peaks.closest_gene.tsv"
done
log "Finished 08.peak_gene_annotation"
