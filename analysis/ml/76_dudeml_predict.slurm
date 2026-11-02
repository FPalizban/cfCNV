#!/bin/bash
#SBATCH --job-name=dude_pred
#SBATCH --cpus-per-task=12
#SBATCH --mem=160G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/dudeML/dude_pred_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/dudeML/dude_pred_%A_%a.err

set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project

source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

module load BEDTools/2.30.0-GCC-11.3.0

TOOLS="$PROJECT/06_tools"
RUN="$TOOLS/bin/run_dudeml"

DUDE="$PROJECT/07_tool_outputs/simulation1/dudeML"

TASKS="$DUDE/metadata/prediction_tasks.tsv"

LINE=$(awk -F'\t' \
-v n=$((SLURM_ARRAY_TASK_ID+2)) \
'NR==n {print; exit}' \
"$TASKS")

[[ -n "$LINE" ]] || exit 2

IFS=$'\t' read -r \
TASK SAMPLE BAM BACKGROUND DEPTH REP CLASS BIN \
<<< "$LINE"

FEATURE="$DUDE/features/$SAMPLE/bin_${BIN}/${SAMPLE}.sample.bed"

OD="$DUDE/predictions/$SAMPLE/bin_${BIN}"

mkdir -p "$OD"

if [[ -f "$OD/.done" ]]
then
    echo "Already complete"
    exit 0
fi

[[ -s "$FEATURE" ]] || {
    echo "ERROR feature missing:"
    echo "$FEATURE"
    exit 10
}


MODEL_DIR="$DUDE/models/depth_${DEPTH}/bin_${BIN}"


if [[ "$CLASS" == "positive" ]]
then

    TRAIN="$MODEL_DIR/fold${REP}.sav"

    [[ -s "$TRAIN" ]] || {
        echo "ERROR held-out model missing:"
        echo "$TRAIN"
        exit 11
    }

else

    # Negative controls:
    # use dudeML's native multi-model consensus.
    TRAIN="$MODEL_DIR"

    NMODELS=$(find "$MODEL_DIR" \
        -maxdepth 1 \
        -type f \
        -name 'fold*.sav' \
        | wc -l)

    [[ "$NMODELS" -eq 5 ]] || {
        echo "ERROR expected five fold models"
        exit 12
    }

fi


PRED="$OD/${SAMPLE}.predict.bed"
SEG="$OD/${SAMPLE}.segments.bed"


echo "======================================"
echo "dudeML prediction"
echo "sample = $SAMPLE"
echo "depth  = $DEPTH"
echo "bin    = $BIN"
echo "class  = $CLASS"
echo "model  = $TRAIN"
echo "======================================"


cd "$OD"


"$RUN" predict \
    -i "$FEATURE" \
    -t "$TRAIN" \
    -o "$PRED"


[[ -s "$PRED" ]] || {
    echo "ERROR prediction output missing"
    exit 20
}


# --------------------------------------------------
# Merge adjacent Dup/Del windows into CNV segments.
#
# Confidence threshold = 0.5.
# This is NOT optimized against benchmark truth.
# --------------------------------------------------

"$RUN" summarize \
    -i "$PRED" \
    -o "$SEG" \
    -c 0.5 \
    -w "$BIN" \
    -id "$SAMPLE"


# Zero-call samples are valid:
[[ -f "$SEG" ]] || {
    echo "ERROR summarized segment file absent"
    exit 21
}


# Preserve result but compress very large per-window prediction.
gzip -f "$PRED"


touch "$OD/.done"

echo
echo "DUDEML PREDICTION COMPLETE"
echo "$SEG"
