#!/usr/bin/env bash
set -euo pipefail
usage(){ cat <<'USAGE'
Usage:
  bash 07.peak_gene_annotation.sh <project_name> <gene.gff3> [project_dir]

Example:
  bash 07.peak_gene_annotation.sh duli_hap1 genome/Pyrbe_hap1.gff3 .

Input:
  04.peak/<project_name>.ATAC_peaks.narrowPeak
Output:
  07.annotation/<project_name>.ATAC_peaks.overlap_genes.tsv
  07.annotation/<project_name>.ATAC_peaks.closest_gene.tsv
USAGE
}
[ "$#" -ge 2 ] || { usage; exit 1; }
project="$1"; gff="$2"; project_dir="${3:-.}"
cd "$project_dir"
mkdir -p 07.annotation logs
LOG="logs/07.peak_gene_annotation.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1
need_file(){ [ -s "$1" ] || { echo "ERROR: missing or empty file: $1" >&2; exit 1; }; }
source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate "${WORK_ENV:-work}"
command -v bedtools
peaks="04.peak/${project}.ATAC_peaks.narrowPeak"
gene_bed="07.annotation/${project}.genes.bed"
need_file "$gff"; need_file "$peaks"
awk 'BEGIN{OFS="\t"} $0 !~ /^#/ && ($3=="gene" || $3=="mRNA") {id="."; if (match($9,/ID=([^;]+)/,a)) id=a[1]; else if (match($9,/Name=([^;]+)/,b)) id=b[1]; print $1,$4-1,$5,id,".",$7}' "$gff" | sort -k1,1 -k2,2n > "$gene_bed"
bedtools intersect -a "$peaks" -b "$gene_bed" -wa -wb > "07.annotation/${project}.ATAC_peaks.overlap_genes.tsv"
bedtools closest -a "$peaks" -b "$gene_bed" -d > "07.annotation/${project}.ATAC_peaks.closest_gene.tsv"
echo "[$(date '+%F %T')] Finished 07.peak_gene_annotation"
