#!/usr/bin/env bash
set -euo pipefail
usage(){ cat <<'USAGE'
Usage:
  bash 01.index.sh <genome.fa> [project_dir]

Example:
  bash 01.index.sh genome/Pyrbe_hap1.fa .

Output:
  <genome.fa>.fai
  <genome.fa>.bwt/.pac/.ann/.amb/.sa
USAGE
}
[ "$#" -ge 1 ] || { usage; exit 1; }
genome_fa="$1"; project_dir="${2:-.}"
cd "$project_dir"
mkdir -p 06.qc logs
LOG="logs/01.index.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1
log(){ echo "[$(date '+%F %T')] $*"; }
need_file(){ [ -s "$1" ] || { echo "ERROR: missing or empty file: $1" >&2; exit 1; }; }
source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate "${WORK_ENV:-work}"
command -v bwa; command -v samtools
need_file "$genome_fa"
[ -s "${genome_fa}.fai" ] || samtools faidx "$genome_fa"
if [ ! -s "${genome_fa}.bwt" ]; then
  log "bwa index ${genome_fa}"
  bwa index "$genome_fa" > "06.qc/bwa_index.log" 2>&1
else
  log "Skip existing bwa index"
fi
log "Finished 01.index"
