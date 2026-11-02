#!/bin/bash

set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project

source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

SIM1_OUT="${SIM1_OUT:-$PROJECT/07_tool_outputs/simulation1}"

NUNIQUE=$(
    find "$SIM1_OUT/BICseq2/unique" \
    -name '.done' \
    -type f \
    | wc -l
)

NFRAG=$(
    find "$SIM1_OUT/BICseq2/unique" \
    -name 'fragment_size.txt' \
    -type f \
    | wc -l
)

NTASK=$(( $(wc -l < "$SIM1_MAN/bicseq2_tasks.tsv") - 1 ))

echo "Unique samples = $NUNIQUE"
echo "Fragment files = $NFRAG"
echo "Tasks          = $NTASK"

[[ "$NUNIQUE" -eq 114 ]] || {
    echo "ERROR: expected 114 completed unique-read samples"
    exit 2
}

[[ "$NFRAG" -eq 114 ]] || {
    echo "ERROR: expected 114 fragment-size files"
    exit 3
}

[[ "$NTASK" -eq 684 ]] || {
    echo "ERROR: expected 684 BIC-seq2 tasks"
    exit 4
}


echo
echo "Submitting resource repair..."

JFIX=$(sbatch --parsable \
"$PROJECT/00_metadata/61b_finish_bicseq2_resources.slurm")


echo "Submitting one 1-Mb smoke test..."

JTEST=$(sbatch --parsable \
--dependency=afterok:$JFIX \
--array=5 \
"$PROJECT/00_metadata/63_bicseq2_normseg.slurm")


echo "Submitting remaining 683 configurations..."

JPROD=$(sbatch --parsable \
--dependency=afterok:$JTEST \
--array=0-4,6-683%12 \
"$PROJECT/00_metadata/63_bicseq2_normseg.slurm")


echo "Submitting final audit..."

JAUDIT=$(sbatch --parsable \
--dependency=afterany:$JPROD \
"$PROJECT/00_metadata/64_bicseq2_audit.slurm")


STATUS="$SIM1_MAN/bicseq2_resume_jobs.tsv"

cat > "$STATUS" <<TABLE
component	job_id
resource_fix	$JFIX
smoke_test	$JTEST
production	$JPROD
audit	$JAUDIT
TABLE

echo
echo "======================================"
echo "BICSEQ2 RESUME PIPELINE SUBMITTED"
echo "======================================"

cat "$STATUS"
