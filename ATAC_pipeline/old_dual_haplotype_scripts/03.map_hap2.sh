#!/usr/bin/env bash
set -euo pipefail
usage(){ cat <<'USAGE'
Usage:
  bash 03.map_hap2.sh <sample_name> <clean_R1.fq.gz> <clean_R2.fq.gz> <hap2.fa> [project_dir]

Example:
  bash 03.map_hap2.sh CRR2978833 00.clean/CRR2978833.clean.r1.fq.gz 00.clean/CRR2978833.clean.r2.fq.gz hap2/Pyrbe_hap2.fa .

Output:
  02.bam/<sample_name>.hap2.sorted.bam
USAGE
}
[ "$#" -ge 4 ] || { usage; exit 1; }
sample="$1"; r1="$2"; r2="$3"; fa="$4"; project_dir="${5:-.}"
cd "$project_dir"
mkdir -p 02.bam 06.qc logs
LOG="logs/03.map_hap2.${sample}.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1
need_file(){ [ -s "$1" ] || { echo "ERROR: missing or empty file: $1" >&2; exit 1; }; }
log(){ echo "[$(date '+%F %T')] $*"; }
source "${HOME}/miniconda3/etc/profile.d/conda.sh"; conda activate "${WORK_ENV:-work}"
command -v bwa; command -v samtools
need_file "$r1"; need_file "$r2"; need_file "$fa"
out="02.bam/${sample}.hap2.sorted.bam"
if [ -s "$out" ] && [ -s "${out}.bai" ] && samtools quickcheck -v "$out"; then log "Skip $out"; else
  log "Map ${sample} to hap2"
  bwa mem -t "${THREADS:-24}" "$fa" "$r1" "$r2" 2> "06.qc/${sample}.hap2.bwa.log" | samtools sort -@ "${SORT_THREADS:-12}" -o "$out" -
  samtools index "$out"
fi
log "Finished 03.map_hap2"
