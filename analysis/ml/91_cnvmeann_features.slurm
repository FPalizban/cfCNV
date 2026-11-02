#!/bin/bash
#SBATCH --job-name=meann_feat
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/meann_feat_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/meann_feat_%A_%a.err
#
# One array task per Simulation-1 BAM (114 tasks: --array=0-113).
# Requires the per-sample manifest made in step C2 of RUN_ORDER.md.

set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project
source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

module load Python/3.10.8-GCCcore-12.2.0
source "$PROJECT/06_tools/envs/cnvmeann/bin/activate"

BIN="${BIN:-1000}"
EXTRA="${EXTRA:-}"            # e.g. EXTRA="--drop-zero" for the faithful variant
REF="$PROJECT/02_reference/GRCh37/fasta/human_g1k_v37.fasta"
MAN="$SIM1_MAN/sim1_samples.tsv"
OUT="$PROJECT/07_tool_outputs/simulation1/CNVMEANN/features"
TRUTHDIR="$SIM/healthy_cfDNA_spikein/truth"

LINE=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+2)) 'NR==n{print;exit}' "$MAN")
[[ -n "$LINE" ]] || { echo "no manifest line for task $SLURM_ARRAY_TASK_ID"; exit 2; }
IFS=$'\t' read -r SAMPLE BAM BACKGROUND DEPTH REPLICATE CLASS <<< "$LINE"

OD="$OUT/bin${BIN}/$SAMPLE"
mkdir -p "$OD"
if [[ -f "$OD/.done.bin${BIN}" ]]; then echo "already done"; exit 0; fi

TRUTHARG=""
if [[ "$CLASS" == "positive" ]]; then
    T="$TRUTHDIR/rep${REPLICATE}.truth.tsv"
    [[ -s "$T" ]] || { echo "missing truth $T"; exit 3; }
    TRUTHARG="--truth $T"
fi

echo "sample=$SAMPLE bg=$BACKGROUND depth=$DEPTH rep=$REPLICATE class=$CLASS bin=$BIN"
/usr/bin/time -v python "$PROJECT/00_metadata/90_cnvmeann_features.py" \
    --bam "$BAM" --ref "$REF" --sample "$SAMPLE" \
    --outdir "$OD" --bin "$BIN" $TRUTHARG $EXTRA
