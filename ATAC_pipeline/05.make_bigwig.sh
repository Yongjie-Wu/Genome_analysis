#!/usr/bin/env bash
set -euo pipefail
usage(){ cat <<'USAGE'
Usage:
  bash 05.make_bigwig.sh <project_name> [project_dir]

Example:
  bash 05.make_bigwig.sh duli_hap1 .

Input:
  03.bam_filter/<project_name>.merge.rmdup.bam
Output:
  05.bigwig/<project_name>.merge.rmdup.CPM.bw
USAGE
}
[ "$#" -ge 1 ] || { usage; exit 1; }
project="$1"; project_dir="${2:-.}"
cd "$project_dir"
mkdir -p 05.bigwig logs
LOG="logs/05.make_bigwig.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1
log(){ echo "[$(date '+%F %T')] $*"; }
need_file(){ [ -s "$1" ] || { echo "ERROR: missing or empty file: $1" >&2; exit 1; }; }
source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate "${DEEPTOOLS_ENV:-deeptools_local}"
command -v bamCoverage; bamCoverage --version
bam="03.bam_filter/${project}.merge.rmdup.bam"
bw="05.bigwig/${project}.merge.rmdup.CPM.bw"
need_file "$bam"
if [ -s "$bw" ]; then
  log "Skip ${bw}"
else
  bamCoverage -b "$bam" -o "$bw" --binSize "${BIGWIG_BIN_SIZE:-10}" --normalizeUsing CPM --numberOfProcessors "${BIGWIG_THREADS:-16}"
fi
log "Finished 05.make_bigwig"
