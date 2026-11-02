#!/bin/bash
#SBATCH --job-name=cfcnv_score
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --time=20:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/score_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/score_%j.err
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/00_env.sh
python3 "$PROJECT/00_metadata/20_collect_calls.py"
python3 "$PROJECT/00_metadata/21_score_sim1.py"
python3 "$PROJECT/00_metadata/22_bin_metrics.py"
echo "RESULTS";ls -lh "$SIM1_MET"
