#!/bin/bash
#SBATCH --job-name=ML_F_report
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/ml_dataset/F_report_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/ml_dataset/F_report_%j.err
# F  manuscript tables + figure (v1 and v2 models next to the 8 callers, manuscript scoring definitions)
#    writes only $MLO/manuscript and $MLO/cfCNV_ML_manuscript_results.tar.gz
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/ml_dataset/ml_env.sh
pyd "$ML/f_report.py"
if pyd -c "import matplotlib" 2>/dev/null; then
  pyd "$ML/g_figure.py" "$MLO/manuscript"
else
  echo "matplotlib not in this python - figure skipped; run g_figure.py on your laptop on the manuscript folder"
fi
cp "$ML/g_figure.py" "$MLO/manuscript/"
cd "$MLO" && tar czf cfCNV_ML_manuscript_results.tar.gz manuscript results/run_info.json results_v2/run_info.json
echo "DONE -> $MLO/cfCNV_ML_manuscript_results.tar.gz"
