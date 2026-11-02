#!/bin/bash
#SBATCH --job-name=tf_ace
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/tf/ace_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/tf/ace_%A_%a.err
#
# One task per row of the tumour-fraction manifest (see RUN_ORDER.md, T1).
# Manifest columns: sample_id bam dataset build true_tf target_depth

set -euo pipefail
PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project
R_MODULE="${R_MODULE:-R}"
BIN_KB="${BIN_KB:-1000}"
MAN="${MAN:-$PROJECT/00_metadata/manifests/tf/tf_samples.tsv}"

LINE=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+2)) 'NR==n{print;exit}' "$MAN")
[[ -n "$LINE" ]] || { echo "no manifest line"; exit 2; }
IFS=$'\t' read -r SAMPLE BAM DATASET BUILD TRUE_TF DEPTH <<< "$LINE"

OD="$PROJECT/07_tool_outputs/tf/ACE/bin${BIN_KB}kb/$DATASET/$SAMPLE"
[[ -f "$OD/.done" ]] && { echo "already done"; exit 0; }

module load "$R_MODULE"
/usr/bin/time -v Rscript "$PROJECT/00_metadata/98_ace.R" \
    "$BAM" "$SAMPLE" "$BUILD" "$BIN_KB" "$OD" 2
