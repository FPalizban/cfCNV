#!/bin/bash
#SBATCH --job-name=ML_A_gc
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=03:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/ml_dataset/A_gc_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/ml_dataset/A_gc_%j.err
# A  GC / N / blacklist annotation of 10-kb bins for GRCh37 and GRCh38 (reads the FASTA files only).
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/ml_dataset/ml_env.sh
pyd "$ML/a_gc_bins.py"
