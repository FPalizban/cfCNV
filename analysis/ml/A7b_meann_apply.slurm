#!/bin/bash
#SBATCH --job-name=S2_meann_apply
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --array=0-5
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/meann_apply_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/meann_apply_%A_%a.err
# A7b  CNV-MEANN v2: train on all Simulation 1 positives, apply to all HG008 mixtures (6 tasks, one per bin).
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/simulation2/full/s2_env.sh
export OMP_NUM_THREADS=4
BB=($BINS)
pyd "$S2/a7_meann_train_apply.py" --bin "${BB[$SLURM_ARRAY_TASK_ID]}"
