#!/bin/bash
set -euo pipefail
PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project
source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"
OUT="$PROJECT/09_metrics/simulation1"
mkdir -p "$OUT"
JFILE="$SIM1_MAN/submitted_jobs.tsv"
{
  echo 'component,JobID,JobName,State,ExitCode,ElapsedRaw,TotalCPU,MaxRSS,ReqMem,AllocCPUS'
  while IFS=$'\t' read -r component jid; do
    [[ "$component" == "component" || -z "${jid:-}" ]] && continue
    sacct -X -j "$jid" -n -P --format=JobIDRaw,JobName,State,ExitCode,ElapsedRaw,TotalCPU,MaxRSS,ReqMem,AllocCPUS | \
      awk -F'|' -v c="$component" 'NF>=9 && $1!="" {gsub(/,/,"_",c); print c","$1","$2","$3","$4","$5","$6","$7","$8","$9}'
  done < "$JFILE"
} > "$OUT/slurm_resource_usage.csv"
echo "Wrote $OUT/slurm_resource_usage.csv"
