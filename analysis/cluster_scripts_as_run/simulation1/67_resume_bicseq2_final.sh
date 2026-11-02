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

[[ "$NUNIQUE" -eq 114 ]]
[[ "$NFRAG" -eq 114 ]]
[[ "$NTASK" -eq 684 ]]

echo
echo "Submitting WIG/resource repair..."

JFIX=$(sbatch --parsable \
"$PROJECT/00_metadata/61c_finish_bicseq2_resources.slurm")

echo "resource_fix = $JFIX"


echo
echo "Submitting 1-Mb smoke test..."

JTEST=$(sbatch --parsable \
--dependency=afterok:$JFIX \
--array=5 \
"$PROJECT/00_metadata/63_bicseq2_normseg.slurm")

echo "smoke_test = $JTEST"


echo
echo "Submitting remaining 683 configurations..."

JPROD=$(sbatch --parsable \
--dependency=afterok:$JTEST \
--array=0-4,6-683%12 \
"$PROJECT/00_metadata/63_bicseq2_normseg.slurm")

echo "production = $JPROD"


echo
echo "Submitting final audit..."

JAUDIT=$(sbatch --parsable \
--dependency=afterany:$JPROD \
"$PROJECT/00_metadata/64_bicseq2_audit.slurm")

echo "audit = $JAUDIT"


STATUS="$SIM1_MAN/bicseq2_final_jobs.tsv"

cat > "$STATUS" <<TABLE
component	job_id
resource_fix	$JFIX
smoke_test	$JTEST
production	$JPROD
audit	$JAUDIT
TABLE

echo
echo "======================================"
echo "BICSEQ2 FINAL PIPELINE SUBMITTED"
echo "======================================"

cat "$STATUS"
