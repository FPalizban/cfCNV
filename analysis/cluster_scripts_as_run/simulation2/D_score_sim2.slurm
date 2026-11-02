#!/bin/bash
#SBATCH --job-name=S2_score
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/score_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/score_%j.err
# D  Simulation 2 scoring: tumour fraction, profile correlation + LoD, event level. Scores
#    whatever has finished; safe to re-run any time (e.g. after BIC-seq2 completes).
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/simulation2/full/s2_env.sh
pyd "$S2/d_score_sim2.py"
