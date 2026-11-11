#!/bin/bash
#SBATCH --job-name=ML_B_count
#SBATCH --cpus-per-task=3
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/ml_dataset/B_count_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/ml_dataset/B_count_%A_%a.err
# B  10-kb read/fragment counts, one array task per BAM (or per few BAMs if there are > 500).
#    The BAMs are opened read-only. Finished samples are skipped when re-submitted.
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/ml_dataset/ml_env.sh
pyd "$ML/b_count_bins.py" "$SLURM_ARRAY_TASK_ID" "$NTASKS_B"
