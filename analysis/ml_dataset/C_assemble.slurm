#!/bin/bash
#SBATCH --job-name=ML_C_assemble
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/ml_dataset/C_assemble_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/ml_dataset/C_assemble_%j.err
# C  features + labels + splits -> $MLO/hf_dataset (Hugging Face layout + dataset card)
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/ml_dataset/ml_env.sh
pyd "$ML/c_assemble.py"
