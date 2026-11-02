#!/usr/bin/env bash
# Train one dudeML random forest (RFC100, dudeML's default) on labelled feature vectors, balanced to
# 10,000 windows per class (normal / Dup / Del), as in the cfCNV paper.
# usage: dudeml_train.sh <model.sav> <seed> <train.bed> [<train.bed> ...]
set -euo pipefail
OUT=$1; SEED=$2; shift 2; D=${DUDEML_DIR:-/opt/dudeML}
printf '%s\n' "$@" > train_files.txt
dudeml_balance.py --files-from train_files.txt --output balanced_train.bed --per-class 10000 --seed "$SEED"
python "$D/dudeML.py" classify -i balanced_train.bed -o "$OUT" -m RFC100
test -s "$OUT" && test -s "${OUT}2"
rm -f balanced_train.bed
