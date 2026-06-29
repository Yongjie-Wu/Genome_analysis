#!/usr/bin/env bash
set -euo pipefail
usage(){ cat <<'USAGE'
Usage:
  bash 09.atac_window_features.sh <project_name> <hap1|hap2> <hap.fa> <hap.gff3> <centromere_range.txt> <SatDNA.bed> <RPKM.tsv> [project_dir]

Example:
  bash 09.atac_window_features.sh duli hap1 hap1/Pyrbe_duli_hap1.fa hap1/Pyrbe_duli_hap1.gff3 Pyrbe_duli.hap1_centromere_range.txt Pyrbe_duli.hap1.SatDNA.bed Pyrbe_duli_hap1.RPKM .

Input:
  03.bam_filter/<project_name>.<hap>.merge.rmdup.bam
Output:
  09.atac_windows/<hap>/features/*ATAC_SatDNA_gene_expression.tsv
  09.atac_windows/<hap>/plots/*ATAC_correlations_scatter.png/pdf
USAGE
}
[ "$#" -ge 7 ] || { usage; exit 1; }
project="$1"; hap="$2"; fa="$3"; gff="$4"; cent="$5"; sat="$6"; rpkm="$7"; project_dir="${8:-.}"
cd "$project_dir"; mkdir -p logs
LOG="logs/09.atac_window_features.${hap}.$(date +%Y%m%d_%H%M%S).log"; exec > >(tee -a "$LOG") 2>&1
source "${HOME}/miniconda3/etc/profile.d/conda.sh"; conda activate "${PY_ENV:-base}"
tmpcfg="logs/09.${project}.${hap}.$$.config.sh"
if [ "$hap" = hap1 ]; then
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
else
cat > "$tmpcfg" <<CFG
PROJECT_NAME="$project"
PROJECT_DIR="$PWD"
SAMPLE_TABLE="unused"
HAP1_FA="$fa"
HAP2_FA="$fa"
HAP1_GFF="$gff"
HAP2_GFF="$gff"
HAP1_CENTROMERE=""
HAP1_SATDNA_BED=""
HAP1_RPKM=""
HAP2_CENTROMERE="$cent"
HAP2_SATDNA_BED="$sat"
HAP2_RPKM="$rpkm"
WINDOW_SIZE="${WINDOW_SIZE:-50000}"
CFG
fi
python /storage3/wuyongjie/Project/centromere/atct/ATAC_pipeline/bin/atac_window_features.py "$tmpcfg" "$hap"
rm -f "$tmpcfg"
