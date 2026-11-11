#!/bin/bash
# 3_run_v2.sh - submits model v2 (augmented training) on the LOGIN node:   bash 3_run_v2.sh
# Needs the finished dataset from step C. Does not change the dataset or the v1 results.
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/ml_dataset/ml_env.sh
[[ -f "$MLO/hf_dataset/README.md" ]] || { echo "dataset not found - run 1_run_all.sh first"; exit 1; }
cd "$ML"
e=$(sbatch --parsable E_augmented.slurm)
echo "model_v2=$e" | tee -a "$MLLOG/ml_jobids.txt"
echo "Submitted job $e. Check with:  squeue -u $USER   (finished when it disappears)"
