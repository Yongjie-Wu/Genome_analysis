#!/usr/bin/env bash
set -euo pipefail
usage(){ cat <<'USAGE'
Usage:
  bash 04.filter_markdup_merge.sh <project_name> <project_dir> <sample1> <sample2> [sample3 ...]

Example:
  bash 04.filter_markdup_merge.sh duli . CRR2978833 CRR2978834

Input:
  02.bam/<sample>.hap1.sorted.bam
  02.bam/<sample>.hap2.sorted.bam
Output:
  03.bam_filter/<sample>.hap1.rmdup.bam
  03.bam_filter/<sample>.hap2.rmdup.bam
  03.bam_filter/<project_name>.hap1.merge.rmdup.bam
  03.bam_filter/<project_name>.hap2.merge.rmdup.bam
USAGE
}
[ "$#" -ge 4 ] || { usage; exit 1; }
project="$1"; project_dir="$2"; shift 2; samples=("$@")
cd "$project_dir"
mkdir -p 03.bam_filter 06.qc logs
LOG="logs/04.filter_markdup_merge.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1
log(){ echo "[$(date '+%F %T')] $*"; }
need_file(){ [ -s "$1" ] || { echo "ERROR: missing or empty file: $1" >&2; exit 1; }; }
source "${HOME}/miniconda3/etc/profile.d/conda.sh"; conda activate "${WORK_ENV:-work}"
command -v samtools
for hap in hap1 hap2; do
  bam_list=()
  for sample in "${samples[@]}"; do
    inbam="02.bam/${sample}.${hap}.sorted.bam"
    out="03.bam_filter/${sample}.${hap}.rmdup.bam"
    need_file "$inbam"
    if [ -s "$out" ] && [ -s "${out}.bai" ] && samtools quickcheck -v "$out"; then
      log "Skip $out"
    else
      log "Process ${sample}.${hap}"
      rm -f "03.bam_filter/${sample}.${hap}."{name,fixmate,fixmate.sorted,rmdup.raw}".bam" "$out" "${out}.bai"
      samtools sort -@ "${SORT_THREADS:-12}" -n -o "03.bam_filter/${sample}.${hap}.name.bam" "$inbam"
      samtools fixmate -@ "${SORT_THREADS:-12}" -m "03.bam_filter/${sample}.${hap}.name.bam" "03.bam_filter/${sample}.${hap}.fixmate.bam"
      samtools sort -@ "${SORT_THREADS:-12}" -o "03.bam_filter/${sample}.${hap}.fixmate.sorted.bam" "03.bam_filter/${sample}.${hap}.fixmate.bam"
      samtools markdup -@ "${SORT_THREADS:-12}" -r "03.bam_filter/${sample}.${hap}.fixmate.sorted.bam" "03.bam_filter/${sample}.${hap}.rmdup.raw.bam"
      samtools view -@ "${SORT_THREADS:-12}" -b -f 2 -q 30 -F 1804 "03.bam_filter/${sample}.${hap}.rmdup.raw.bam" -o "$out"
      samtools index "$out"
      samtools flagstat -@ "${SORT_THREADS:-12}" "$out" > "06.qc/${sample}.${hap}.rmdup.flagstat.txt"
      rm -f "03.bam_filter/${sample}.${hap}."{name,fixmate,fixmate.sorted,rmdup.raw}".bam"
    fi
    bam_list+=("$out")
  done
  merged="03.bam_filter/${project}.${hap}.merge.rmdup.bam"
  if [ -s "$merged" ] && [ -s "${merged}.bai" ] && samtools quickcheck -v "$merged"; then log "Skip $merged"; else
    log "Merge ${hap}"
    samtools merge -@ "${SORT_THREADS:-12}" -f "$merged" "${bam_list[@]}"
    samtools index "$merged"
  fi
  samtools flagstat -@ "${SORT_THREADS:-12}" "$merged" > "06.qc/${project}.${hap}.merge.rmdup.flagstat.txt"
done
log "Finished 04.filter_markdup_merge"
