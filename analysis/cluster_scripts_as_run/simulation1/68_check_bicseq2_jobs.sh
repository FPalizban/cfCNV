#!/bin/bash

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project

source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

JOBFILE="$SIM1_MAN/bicseq2_final_jobs.tsv"

echo "=========================================="
echo "BIC-seq2 JOB MONITOR"
echo "=========================================="
date
echo

if [[ ! -f "$JOBFILE" ]]; then
    echo "ERROR: job file not found:"
    echo "$JOBFILE"
    exit 1
fi

cat "$JOBFILE"

echo
echo "=========================================="
echo "CURRENTLY IN SLURM QUEUE"
echo "=========================================="

squeue -u "$USER" \
-o "%.20i %.22j %.10T %.12M %.6D %R" \
| grep -E 'JOBID|bic_mapfix|bic_normseg|bic_audit'

echo
echo "=========================================="
echo "ACCOUNTING STATUS"
echo "=========================================="

IDS=$(awk 'NR>1 {print $2}' "$JOBFILE" | paste -sd, -)

sacct -j "$IDS" \
--format=JobID%24,JobName%22,State%20,ExitCode%10,Elapsed%12,MaxRSS%12

echo
echo "=========================================="
echo "ARRAY TASK SUMMARY"
echo "=========================================="

JPROD=$(awk -F'\t' '$1=="production"{print $2}' "$JOBFILE")

if [[ -n "$JPROD" ]]; then

    sacct -n -X -j "$JPROD" \
    --format=State \
    | sed 's/[[:space:]]//g' \
    | grep -v '^$' \
    | sort \
    | uniq -c \
    | sort -nr

fi

echo
echo "=========================================="
echo "OUTPUT COUNTS"
echo "=========================================="

SIM1_OUT="${SIM1_OUT:-$PROJECT/07_tool_outputs/simulation1}"

echo -n "Completed .done configurations : "
find "$SIM1_OUT/BICseq2" \
-type f -name '.done' \
2>/dev/null | grep '/bin_' | wc -l

echo -n "CNV output files              : "
find "$SIM1_OUT/BICseq2" \
-type f -name '*.cnv' \
2>/dev/null | wc -l

echo -n "Unique-read samples completed : "
find "$SIM1_OUT/BICseq2/unique" \
-type f -name '.done' \
2>/dev/null | wc -l

echo
echo "Expected final values:"
echo "  CNV configurations = 684"
echo "  Unique samples     = 114"

echo
echo "=========================================="
