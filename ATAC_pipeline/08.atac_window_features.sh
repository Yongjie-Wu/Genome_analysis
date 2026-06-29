#!/usr/bin/env bash
set -euo pipefail
usage(){ cat <<'USAGE'
Usage:
  bash 08.atac_window_features.sh <project_name> <genome.fa> <gene.gff3> <centromere_range.txt> <SatDNA.bed> <RPKM.tsv> [project_dir]

Example:
  bash 08.atac_window_features.sh duli_hap1 genome/Pyrbe_hap1.fa genome/Pyrbe_hap1.gff3 hap1_centromere_range.txt hap1.SatDNA.bed hap1.RPKM .

Input:
  03.bam_filter/<project_name>.merge.rmdup.bam
Output:
  08.atac_windows/features/*ATAC_SatDNA_gene_expression.tsv
  08.atac_windows/plots/*ATAC_correlations_scatter.png/pdf
USAGE
}
[ "$#" -ge 6 ] || { usage; exit 1; }
project="$1"; fa="$2"; gff="$3"; cent="$4"; sat="$5"; rpkm="$6"; project_dir="${7:-.}"
cd "$project_dir"
mkdir -p logs
LOG="logs/08.atac_window_features.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1
source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate "${PY_ENV:-base}"
tmpcfg="logs/08.${project}.$$.config.sh"
cat > "$tmpcfg" <<CFG
PROJECT_NAME="$project"
PROJECT_DIR="$PWD"
SAMPLE_TABLE="unused"
HAP1_FA="$fa"
HAP2_FA="$fa"
HAP1_GFF="$gff"
HAP2_GFF="$gff"
HAP1_CENTROMERE="$cent"
HAP1_SATDNA_BED="$sat"
HAP1_RPKM="$rpkm"
HAP2_CENTROMERE=""
HAP2_SATDNA_BED=""
HAP2_RPKM=""
WINDOW_SIZE="${WINDOW_SIZE:-50000}"
CFG
python /storage3/wuyongjie/Project/centromere/atct/ATAC_pipeline/bin/atac_window_features.py "$tmpcfg" hap1
rm -f "$tmpcfg"
echo "[$(date '+%F %T')] Finished 08.atac_window_features"
