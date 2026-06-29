#!/usr/bin/env bash
set -euo pipefail
usage(){ cat <<'USAGE'
Usage:
  bash 04.macs2_callpeak.sh <project_name> <genome_size> [project_dir]

Example:
  bash 04.macs2_callpeak.sh duli_hap1 5e8 .

Input:
  03.bam_filter/<project_name>.merge.rmdup.bam
Output:
  04.peak/<project_name>.ATAC_peaks.narrowPeak
USAGE
}
[ "$#" -ge 2 ] || { usage; exit 1; }
project="$1"; genome_size="$2"; project_dir="${3:-.}"
cd "$project_dir"
mkdir -p 04.peak 06.qc logs
LOG="logs/04.macs2_callpeak.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1
log(){ echo "[$(date '+%F %T')] $*"; }
need_file(){ [ -s "$1" ] || { echo "ERROR: missing or empty file: $1" >&2; exit 1; }; }
source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate "${MACS2_ENV:-macs2_clean}"
command -v macs2; macs2 --version
bam="03.bam_filter/${project}.merge.rmdup.bam"
peak="04.peak/${project}.ATAC_peaks.narrowPeak"
need_file "$bam"
if [ -s "$peak" ]; then
  log "Skip ${peak}"
else
  log "MACS2 callpeak"
  macs2 callpeak -t "$bam" -f BAMPE -g "$genome_size" -n "${project}.ATAC" \
    --outdir 04.peak --nomodel --shift -100 --extsize 200 -q 0.01 \
    > "06.qc/${project}.macs2.log" 2>&1
fi
{
  echo -e "Sample\tPeaks\tSummits\tPeak_file"
  echo -e "${project}\t$(wc -l < "04.peak/${project}.ATAC_peaks.narrowPeak")\t$(wc -l < "04.peak/${project}.ATAC_summits.bed")\t04.peak/${project}.ATAC_peaks.narrowPeak"
} > 06.qc/peak_summary.txt
log "Finished 04.macs2_callpeak"
