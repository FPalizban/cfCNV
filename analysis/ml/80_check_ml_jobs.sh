#!/bin/bash
# =====================================================================
# cfCNV - monitor any ML-caller Slurm jobs by job-name prefix
#
# Usage:
#   bash 80_check_ml_jobs.sh            # dudeML jobs (prefix dude_)
#   bash 80_check_ml_jobs.sh meann_     # CNV-MEANN jobs
#   bash 80_check_ml_jobs.sh lyc_       # LYCEUM jobs
#   watch -n 60 "bash 80_check_ml_jobs.sh dude_"
# =====================================================================

PREFIX="${1:-dude_}"
SINCE="${2:-$(date -d '10 days ago' +%F)}"

echo "=================================================="
echo " ML job monitor  prefix=${PREFIX}  since=${SINCE}"
echo " $(date)"
echo "=================================================="

echo
echo "---- In queue now ----"
squeue -u "$USER" -h -o "%i|%j|%T|%M|%R" \
| awk -F'|' -v p="$PREFIX" '$2 ~ "^"p {printf "%-24s %-14s %-10s %-12s %s\n",$1,$2,$3,$4,$5}'

echo
echo "---- Counts per job name and state (finished + running) ----"
sacct -u "$USER" -S "$SINCE" -X -n -P --format=JobName,State \
| awk -F'|' -v p="$PREFIX" '$1 ~ "^"p {split($2,s," "); c[$1"\t"s[1]]++}
  END {for (k in c) print c[k] "\t" k}' \
| sort -t$'\t' -k2,2 -k3,3

echo
echo "---- Tasks that did NOT complete (FAILED / TIMEOUT / OOM / CANCELLED) ----"
sacct -u "$USER" -S "$SINCE" -X -n -P \
      --format=JobID,JobName,State,ExitCode,Elapsed \
| awk -F'|' -v p="$PREFIX" '$2 ~ "^"p && $3 !~ /^(COMPLETED|RUNNING|PENDING)/ {
      printf "%-24s %-14s %-22s %-8s %s\n",$1,$2,$3,$4,$5; n++}
  END {if (!n) print "none"}'

echo
echo "Tip: an array waiting with (JobArrayTaskLimit) is only throttled."
echo "To allow more tasks at once:  scontrol update JobId=<ARRAY_JOB_ID> ArrayTaskThrottle=12"
