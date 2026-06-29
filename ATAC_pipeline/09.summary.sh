#!/usr/bin/env bash
set -euo pipefail
usage(){ cat <<'USAGE'
Usage:
  bash 09.summary.sh <project_name> [project_dir]

Example:
  bash 09.summary.sh duli_hap1 .

Output:
  09.summary/<project_name>.pipeline_outputs.txt
USAGE
}
[ "$#" -ge 1 ] || { usage; exit 1; }
project="$1"; project_dir="${2:-.}"
cd "$project_dir"
mkdir -p 09.summary logs
out="09.summary/${project}.pipeline_outputs.txt"
{
  echo "Project: $project"
  echo "Date: $(date '+%F %T')"
  echo
  echo "Clean FASTQ:"; find 00.clean -type f -name '*.fq.gz' 2>/dev/null | sort
  echo
  echo "BAM:"; find 02.bam 03.bam_filter -type f \( -name '*.bam' -o -name '*.bai' \) 2>/dev/null | sort
  echo
  echo "Peaks:"; find 04.peak -type f 2>/dev/null | sort
  echo
  echo "bigWig:"; find 05.bigwig -type f 2>/dev/null | sort
  echo
  echo "QC:"; find 06.qc -type f 2>/dev/null | sort
  echo
  echo "Annotation:"; find 07.annotation -type f 2>/dev/null | sort
  echo
  echo "Window features:"; find 08.atac_windows -type f 2>/dev/null | sort
} > "$out"
echo "$out"
