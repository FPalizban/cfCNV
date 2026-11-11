#!/bin/bash
#SBATCH --job-name=ML_E_aug
#SBATCH --cpus-per-task=6
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/ml_dataset/E_aug_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/ml_dataset/E_aug_%j.err
# E  model v2: tumour-fraction + depth augmentation from the training negatives, retrain, score all splits
#    reads $MLO/hf_dataset (step C); writes only $MLO/augmented, $MLO/results_v2, $MLO/ml_results_v2_for_claude.tar.gz
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/ml_dataset/ml_env.sh
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
pyd "$ML/e_augmented.py"
