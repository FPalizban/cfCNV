#!/bin/bash
#SBATCH --job-name=P3_tasks
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/phase3/tasks_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/phase3/tasks_%j.err
# 13  After all samples are aligned and derived: BAM manifest, QC table, read length, tool task lists.
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/phase3/full/p3_env.sh
"$SYSPY" "$S2/p2_make_tasks.py"
echo STAGE1_COMPLETE > "$S2_META/STAGE1_COMPLETE"
