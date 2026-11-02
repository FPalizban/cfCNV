#!/bin/bash
#SBATCH --job-name=S2_dude_pred
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/dude_pred_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/dude_pred_%A_%a.err
# A6b  dudeML prediction with the Simulation 1 models of the SAME depth and bin size
#      (all 5 fold models together, as for the Simulation 1 negative controls). dudeML
#      features are coverage relative to the sample, so models carry over between builds.
#      This tests a model trained on simulated germline-like CNVs on real tumour CNVs. 624 tasks.
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/simulation2/full/s2_env.sh
module load BEDTools/2.30.0-GCC-11.3.0
IFS=$'\t' read -r S BAM D B < <(task_line "$S2_META/tasks_dudeml_pred.tsv" "$SLURM_ARRAY_TASK_ID")
FEAT="$S2_OUT/dudeML/features/$S/bin_$B/$S.sample.bed"; test -s "$FEAT"
MODELS="$SIM1_DUDE_MODELS/depth_$D/bin_$B"
[[ $(ls "$MODELS"/fold*.sav 2>/dev/null | wc -l) -eq 5 ]] || { echo "no 5 Sim1 models in $MODELS"; exit 3; }
OD="$S2_OUT/dudeML/predictions/$S/bin_$B"; mkdir -p "$OD"
[[ -f "$OD/.done" ]] && { echo "SKIP $S $B"; exit 0; }
cd "$OD"
"$TOOLS/bin/run_dudeml" predict -i "$FEAT" -t "$MODELS" -o "$OD/$S.predict.bed"
"$TOOLS/bin/run_dudeml" summarize -i "$OD/$S.predict.bed" -o "$OD/$S.segments.bed" -c 0.5 -w "$B" -id "$S"
gzip -f "$OD/$S.predict.bed"
touch "$OD/.done"; echo "DUDE_PRED_DONE $S $B"
