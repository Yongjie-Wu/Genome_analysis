#!/usr/bin/env bash
set -euo pipefail
usage(){ cat <<'USAGE'
Usage:
  bash 05.macs2_callpeak.sh <project_name> <genome_size> [project_dir]

Example:
  bash 05.macs2_callpeak.sh duli 5e8 .

Input:
  03.bam_filter/<project_name>.hap1.merge.rmdup.bam
  03.bam_filter/<project_name>.hap2.merge.rmdup.bam
Output:
  04.peak/hap1/<project_name>.hap1.ATAC_peaks.narrowPeak
  04.peak/hap2/<project_name>.hap2.ATAC_peaks.narrowPeak
USAGE
}
[ "$#" -ge 2 ] || { usage; exit 1; }
project="$1"; genome_size="$2"; project_dir="${3:-.}"
cd "$project_dir"; mkdir -p 04.peak/hap1 04.peak/hap2 06.qc logs
LOG="logs/05.macs2_callpeak.$(date +%Y%m%d_%H%M%S).log"; exec > >(tee -a "$LOG") 2>&1
log(){ echo "[$(date '+%F %T')] $*"; }; need_file(){ [ -s "$1" ] || { echo "ERROR: missing or empty file: $1" >&2; exit 1; }; }
source "${HOME}/miniconda3/etc/profile.d/conda.sh"; conda activate "${MACS2_ENV:-macs2_clean}"
command -v macs2; macs2 --version
for hap in hap1 hap2; do
  bam="03.bam_filter/${project}.${hap}.merge.rmdup.bam"; peak="04.peak/${hap}/${project}.${hap}.ATAC_peaks.narrowPeak"
  need_file "$bam"
  if [ -s "$peak" ]; then log "Skip $peak"; else
    log "MACS2 ${hap}"
    macs2 callpeak -t "$bam" -f BAMPE -g "$genome_size" -n "${project}.${hap}.ATAC" --outdir "04.peak/${hap}" --nomodel --shift -100 --extsize 200 -q 0.01 > "06.qc/${project}.${hap}.macs2.log" 2>&1
  fi
done
log "Finished 05.macs2_callpeak"
