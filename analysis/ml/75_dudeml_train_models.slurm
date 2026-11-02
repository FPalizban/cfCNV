#!/bin/bash
#SBATCH --job-name=dude_train
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/dudeML/dude_train_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/dudeML/dude_train_%A_%a.err

set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project

source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

TOOLS="$PROJECT/06_tools"
RUN="$TOOLS/bin/run_dudeml"

DUDE="$PROJECT/07_tool_outputs/simulation1/dudeML"

TASKS="$DUDE/metadata/model_tasks.tsv"

LINE=$(awk -F'\t' \
-v n=$((SLURM_ARRAY_TASK_ID+2)) \
'NR==n {print; exit}' \
"$TASKS")

[[ -n "$LINE" ]] || exit 2

IFS=$'\t' read -r \
TASK DEPTH BIN FOLD \
<<< "$LINE"

MODEL_DIR="$DUDE/models/depth_${DEPTH}/bin_${BIN}"
TRAIN_DIR="$DUDE/training/task_${TASK}"

mkdir -p \
"$MODEL_DIR" \
"$TRAIN_DIR"

MODEL="$MODEL_DIR/fold${FOLD}.sav"

if [[ -s "$MODEL" && -s "${MODEL}2" ]]
then
    echo "Model already complete"
    exit 0
fi


echo "======================================"
echo "dudeML MODEL TRAINING"
echo "task        = $TASK"
echo "depth       = $DEPTH"
echo "bin         = $BIN"
echo "heldout rep = $FOLD"
echo "======================================"


FILES="$TRAIN_DIR/train_files.txt"

: > "$FILES"


# --------------------------------------------------
# Select POSITIVE feature vectors:
#
# same depth
# same bin
# all backgrounds available at that depth
# EXCLUDING held-out replicate
# --------------------------------------------------

awk -F'\t' \
-v depth="$DEPTH" \
-v fold="$FOLD" '
NR>1 &&
$6=="positive" &&
$4==depth &&
$5!=fold {
    print $1
}' \
"$SIM1_MAN/master.tsv" \
| while read -r SAMPLE
do

    F="$DUDE/features/$SAMPLE/bin_${BIN}/${SAMPLE}.train.bed"

    [[ -s "$F" ]] || {
        echo "ERROR missing training file: $F" >&2
        exit 20
    }

    echo "$F" >> "$FILES"

done


NFILES=$(wc -l < "$FILES")

echo "Training samples = $NFILES"

[[ "$NFILES" -ge 4 ]] || {
    echo "ERROR: insufficient training samples"
    exit 21
}


BAL="$TRAIN_DIR/balanced_train.bed"

python \
"$PROJECT/00_metadata/74_balance_dudeml.py" \
--files-from "$FILES" \
--output "$BAL" \
--per-class 10000 \
--seed $((2026 + TASK))


echo
echo "Balanced classes:"

cut -f4 "$BAL" \
| sort \
| uniq -c


# --------------------------------------------------
# Train dudeML Random Forest
#
# RFC100 is dudeML's default model.
# --------------------------------------------------

cd "$TRAIN_DIR"

"$RUN" classify \
    -i "$BAL" \
    -o "$MODEL" \
    -m RFC100


[[ -s "$MODEL" ]] || {
    echo "ERROR primary model not generated"
    exit 30
}

[[ -s "${MODEL}2" ]] || {
    echo "ERROR ploidy model not generated"
    exit 31
}

echo
echo "DUDEML MODEL COMPLETE"
echo "$MODEL"
