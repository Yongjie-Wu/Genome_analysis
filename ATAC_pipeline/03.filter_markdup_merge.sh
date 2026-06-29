#!/usr/bin/env bash
set -euo pipefail
usage(){ cat <<'USAGE'
Usage:
  bash 03.filter_markdup_merge.sh <project_name> <project_dir> <sample1> <sample2> [sample3 ...]

Example:
  bash 03.filter_markdup_merge.sh duli_hap1 . CRR2978833 CRR2978834

Input:
  02.bam/<sample>.sorted.bam
Output:
  03.bam_filter/<sample>.rmdup.bam
  03.bam_filter/<project_name>.merge.rmdup.bam
USAGE
}
[ "$#" -ge 4 ] || { usage; exit 1; }
project="$1"; project_dir="$2"; shift 2; samples=("$@")
cd "$project_dir"
mkdir -p 03.bam_filter 06.qc logs
LOG="logs/03.filter_markdup_merge.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1
log(){ echo "[$(date '+%F %T')] $*"; }
need_file(){ [ -s "$1" ] || { echo "ERROR: missing or empty file: $1" >&2; exit 1; }; }
source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate "${WORK_ENV:-work}"
command -v samtools
bam_list=()
for sample in "${samples[@]}"; do
  inbam="02.bam/${sample}.sorted.bam"
  out="03.bam_filter/${sample}.rmdup.bam"
  need_file "$inbam"
  if [ -s "$out" ] && [ -s "${out}.bai" ] && samtools quickcheck -v "$out"; then
    log "Skip ${out}"
  else
    log "Process ${sample}"
    rm -f "03.bam_filter/${sample}."{name,fixmate,fixmate.sorted,rmdup.raw}".bam" "$out" "${out}.bai"
    samtools sort -@ "${SORT_THREADS:-12}" -n -o "03.bam_filter/${sample}.name.bam" "$inbam"
    samtools fixmate -@ "${SORT_THREADS:-12}" -m "03.bam_filter/${sample}.name.bam" "03.bam_filter/${sample}.fixmate.bam"
    samtools sort -@ "${SORT_THREADS:-12}" -o "03.bam_filter/${sample}.fixmate.sorted.bam" "03.bam_filter/${sample}.fixmate.bam"
    samtools markdup -@ "${SORT_THREADS:-12}" -r "03.bam_filter/${sample}.fixmate.sorted.bam" "03.bam_filter/${sample}.rmdup.raw.bam"
    samtools view -@ "${SORT_THREADS:-12}" -b -f 2 -q 30 -F 1804 "03.bam_filter/${sample}.rmdup.raw.bam" -o "$out"
    samtools index "$out"
    samtools flagstat -@ "${SORT_THREADS:-12}" "$out" > "06.qc/${sample}.rmdup.flagstat.txt"
    rm -f "03.bam_filter/${sample}."{name,fixmate,fixmate.sorted,rmdup.raw}".bam"
  fi
  bam_list+=("$out")
done
merged="03.bam_filter/${project}.merge.rmdup.bam"
if [ -s "$merged" ] && [ -s "${merged}.bai" ] && samtools quickcheck -v "$merged"; then
  log "Skip ${merged}"
else
  log "Merge replicates"
  samtools merge -@ "${SORT_THREADS:-12}" -f "$merged" "${bam_list[@]}"
  samtools index "$merged"
fi
samtools flagstat -@ "${SORT_THREADS:-12}" "$merged" > "06.qc/${project}.merge.rmdup.flagstat.txt"
log "Finished 03.filter_markdup_merge"
