#!/usr/bin/env bash
set -euo pipefail
usage(){ cat <<'USAGE'
Usage:
  bash 00.clean.sh <sample_name> <R1.fq.gz> <R2.fq.gz> [project_dir]

Example:
  bash 00.clean.sh CRR2978833 CRR2978833_r1.fq.gz CRR2978833_r2.fq.gz .

Output:
  00.clean/<sample_name>.clean.r1.fq.gz
  00.clean/<sample_name>.clean.r2.fq.gz
  06.qc/<sample_name>.fastp.html/json/log
USAGE
}
[ "$#" -ge 3 ] || { usage; exit 1; }
sample="$1"; r1="$2"; r2="$3"; project_dir="${4:-.}"
cd "$project_dir"
mkdir -p 00.clean 06.qc logs
LOG="logs/00.clean.${sample}.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1
log(){ echo "[$(date '+%F %T')] $*"; }
need_file(){ [ -s "$1" ] || { echo "ERROR: missing or empty file: $1" >&2; exit 1; }; }
gzip_head_ok(){ [ "$(od -An -tx1 -N2 "$1" 2>/dev/null | tr -d ' \n')" = "1f8b" ]; }
need_file "$r1"; need_file "$r2"
source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate "${WORK_ENV:-work}"
command -v fastp
out1="00.clean/${sample}.clean.r1.fq.gz"
out2="00.clean/${sample}.clean.r2.fq.gz"
if [ -s "$out1" ] && [ -s "$out2" ] && gzip_head_ok "$out1" && gzip_head_ok "$out2"; then
  log "Skip existing clean reads for ${sample}"
else
  log "fastp ${sample}"
  rm -f "${out1}.tmp" "${out2}.tmp"
  fastp -i "$r1" -I "$r2" -o "${out1}.tmp" -O "${out2}.tmp" \
    -h "06.qc/${sample}.fastp.html" -j "06.qc/${sample}.fastp.json" \
    -w "${FASTP_THREADS:-8}" > "06.qc/${sample}.fastp.log" 2>&1
  mv -f "${out1}.tmp" "$out1"
  mv -f "${out2}.tmp" "$out2"
fi
log "Finished 00.clean"
