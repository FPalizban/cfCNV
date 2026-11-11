#!/bin/bash
#SBATCH --job-name=ML_D_baseline
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/ml_dataset/D_baseline_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/ml_dataset/D_baseline_%j.err
# D  baseline models (threshold + gradient boosting): train, validate, test; writes the small bundle to send back
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/ml_dataset/ml_env.sh
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
pyd "$ML/d_baseline.py"
