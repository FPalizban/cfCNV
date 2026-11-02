#!/bin/bash
#SBATCH --job-name=tf_fragle
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/tf/fragle_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/tf/fragle_%A_%a.err
#
# Fragle, raw-WGS mode (R), one BAM per task.
# Valid ONLY for real plasma cfDNA (incl. Sim-1 healthy BAMs as true-TF=0
# negative controls). Do NOT run on HG008 mixtures: sheared cell-line DNA
# has no plasma fragment-length profile. The script refuses dataset=HG008.

set -euo pipefail
PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project
FR="$PROJECT/06_tools/ML_based/FRAGLE"
MAN="${MAN:-$PROJECT/00_metadata/manifests/tf/tf_samples.tsv}"

LINE=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+2)) 'NR==n{print;exit}' "$MAN")
[[ -n "$LINE" ]] || { echo "no manifest line"; exit 2; }
IFS=$'\t' read -r SAMPLE BAM DATASET BUILD TRUE_TF DEPTH <<< "$LINE"

if [[ "$DATASET" == HG008* ]]; then
    echo "SKIP: Fragle is not run on HG008 mixtures (no plasma fragment profile)"
    exit 0
fi
[[ -s "${BAM}.bai" || -s "${BAM%.bam}.bai" ]] || { echo "ERROR: no BAM index for $BAM"; exit 3; }

OD="$PROJECT/07_tool_outputs/tf/Fragle/$DATASET/$SAMPLE"
[[ -f "$OD/.done" ]] && { echo "already done"; exit 0; }
mkdir -p "$OD"

source "$FR/env/bin/activate"
cd "$FR"
/usr/bin/time -v python main.py --input "$BAM" --output "$OD" \
    --mode R --genome_build "$BUILD" --cpu "$SLURM_CPUS_PER_TASK"
source "$FR/env/bin/deactivate"

[[ -s "$OD/Fragle.csv" ]] || { echo "ERROR: Fragle.csv missing"; exit 4; }
touch "$OD/.done"
echo "FRAGLE DONE $SAMPLE"; cat "$OD/Fragle.csv"
