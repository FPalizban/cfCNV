#!/bin/bash
# =====================================================================
# 4_final.sh - THE ONE COMMAND for the final ML results (login node):   bash 4_final.sh
#   E  model v2 (augmented training) -> F  manuscript tables + figure  (F waits for E automatically)
# Uses the dataset and v1 results already built (steps A-D). Writes only inside $P/13_ml_dataset.
# Result: $P/13_ml_dataset/cfCNV_ML_manuscript_results.tar.gz
# =====================================================================
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/ml_dataset/ml_env.sh
cd "$ML"
for f in "$MLO/hf_dataset/README.md" "$MLO/results/run_info.json" "$MLO/results/predictions/test_hg008.csv.gz"; do
  [[ -f "$f" ]] || { echo "missing $f - steps A-D (bash 1_run_all.sh) must finish first"; exit 1; }
done
if [[ -f "$MLO/results_v2/run_info.json" ]]; then
  echo "model v2 already finished - only the report is submitted"
  f=$(sbatch --parsable F_report.slurm); echo "report=$f" | tee -a "$MLLOG/ml_jobids.txt"
else
  e=$(sbatch --parsable E_augmented.slurm)
  f=$(sbatch --parsable --dependency=afterok:$e F_report.slurm)
  echo "model_v2=$e report=$f" | tee -a "$MLLOG/ml_jobids.txt"
fi
echo "Submitted. Finished when 'squeue -u $USER' shows no ML_ jobs and this file exists:"
echo "   $MLO/cfCNV_ML_manuscript_results.tar.gz"
