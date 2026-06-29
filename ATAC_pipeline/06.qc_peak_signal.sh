#!/usr/bin/env bash
set -euo pipefail
usage(){ cat <<'USAGE'
Usage:
  bash 06.qc_peak_signal.sh <project_name> <project_dir> <sample1> <sample2> [sample3 ...]

Example:
  bash 06.qc_peak_signal.sh duli_hap1 . CRR2978833 CRR2978834

Input:
  03.bam_filter/<sample>.rmdup.bam
  04.peak/<project_name>.ATAC_peaks.narrowPeak
  05.bigwig/<project_name>.merge.rmdup.CPM.bw
Output:
  06.qc/*correlation*.pdf, *PCA*.pdf, *heatmap*.pdf, *profile*.pdf
USAGE
}
[ "$#" -ge 4 ] || { usage; exit 1; }
project="$1"; project_dir="$2"; shift 2; samples=("$@")
cd "$project_dir"
mkdir -p 06.qc logs
LOG="logs/06.qc_peak_signal.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1
need_file(){ [ -s "$1" ] || { echo "ERROR: missing or empty file: $1" >&2; exit 1; }; }
source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate "${DEEPTOOLS_ENV:-deeptools_local}"
command -v multiBamSummary; command -v plotCorrelation; command -v plotPCA; command -v computeMatrix; command -v plotHeatmap; command -v plotProfile
rep_bams=(); labels=()
for sample in "${samples[@]}"; do
  b="03.bam_filter/${sample}.rmdup.bam"
  need_file "$b"
  rep_bams+=("$b")
  labels+=("$sample")
done
npz="06.qc/${project}.replicates.multiBamSummary.${SUMMARY_BIN_SIZE:-10000}bp.npz"
raw="06.qc/${project}.replicates.multiBamSummary.${SUMMARY_BIN_SIZE:-10000}bp.rawCounts.tsv"
[ -s "$npz" ] || multiBamSummary bins -b "${rep_bams[@]}" --labels "${labels[@]}" --binSize "${SUMMARY_BIN_SIZE:-10000}" --numberOfProcessors "${BIGWIG_THREADS:-16}" --outFileName "$npz" --outRawCounts "$raw"
[ -s "06.qc/${project}.replicate_correlation.pearson.heatmap.pdf" ] || plotCorrelation -in "$npz" --corMethod pearson --skipZeros --whatToPlot heatmap --plotFile "06.qc/${project}.replicate_correlation.pearson.heatmap.pdf" --outFileCorMatrix "06.qc/${project}.replicate_correlation.pearson.tsv"
[ -s "06.qc/${project}.replicate_PCA.pdf" ] || plotPCA -in "$npz" --plotFile "06.qc/${project}.replicate_PCA.pdf" --outFileNameData "06.qc/${project}.replicate_PCA.tsv"
bw="05.bigwig/${project}.merge.rmdup.CPM.bw"
peaks="04.peak/${project}.ATAC_peaks.narrowPeak"
matrix="06.qc/${project}.ATAC_peak_center.matrix.gz"
need_file "$bw"; need_file "$peaks"
[ -s "$matrix" ] || computeMatrix reference-point --referencePoint center -b "${PEAK_WINDOW:-3000}" -a "${PEAK_WINDOW:-3000}" -R "$peaks" -S "$bw" --skipZeros --missingDataAsZero --numberOfProcessors "${BIGWIG_THREADS:-16}" --outFileName "$matrix" --outFileNameMatrix "06.qc/${project}.ATAC_peak_center.matrix.tsv"
[ -s "06.qc/${project}.ATAC_peak_center.heatmap.pdf" ] || plotHeatmap -m "$matrix" --outFileName "06.qc/${project}.ATAC_peak_center.heatmap.pdf" --refPointLabel "peak center" --samplesLabel "${project} CPM"
[ -s "06.qc/${project}.ATAC_peak_center.profile.pdf" ] || plotProfile -m "$matrix" --outFileName "06.qc/${project}.ATAC_peak_center.profile.pdf" --refPointLabel "peak center" --samplesLabel "${project} CPM"
echo "[$(date '+%F %T')] Finished 06.qc_peak_signal"
