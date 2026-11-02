#!/bin/bash

set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project

source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

echo "======================================"
echo "Submitting BIC-seq2 Simulation-1"
echo "======================================"

[[ $(wc -l < "$SIM1_MAN/master.tsv") -eq 115 ]] || {
    echo "ERROR: master.tsv should have 115 lines"
    exit 2
}

[[ $(wc -l < "$SIM1_MAN/bicseq2_tasks.tsv") -eq 685 ]] || {
    echo "ERROR: bicseq2_tasks.tsv should have 685 lines"
    exit 3
}


JRES=$(sbatch --parsable \
"$PROJECT/00_metadata/61_bicseq2_resources.slurm")

JUNIQ=$(sbatch --parsable \
"$PROJECT/00_metadata/62_bicseq2_unique.slurm")


# first sample, 1-Mb configuration
JTEST=$(sbatch --parsable \
--dependency=afterok:$JRES:$JUNIQ \
--array=5 \
"$PROJECT/00_metadata/63_bicseq2_normseg.slurm")


JPROD=$(sbatch --parsable \
--dependency=afterok:$JTEST \
--array=0-4,6-683%12 \
"$PROJECT/00_metadata/63_bicseq2_normseg.slurm")


JAUDIT=$(sbatch --parsable \
--dependency=afterany:$JPROD \
"$PROJECT/00_metadata/64_bicseq2_audit.slurm")


cat > "$SIM1_MAN/bicseq2_jobs.tsv" <<TABLE
component	job_id
resources	$JRES
unique_reads	$JUNIQ
smoke_test	$JTEST
production	$JPROD
audit	$JAUDIT
TABLE


echo
echo "======================================"
echo "SUBMITTED"
echo "======================================"

cat "$SIM1_MAN/bicseq2_jobs.tsv"
