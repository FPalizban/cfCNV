#!/bin/bash

set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project

source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

DUDE="$PROJECT/07_tool_outputs/simulation1/dudeML"


echo "======================================"
echo "Submitting dudeML Simulation-1"
echo "======================================"


# --------------------------------------------------
# FEATURE SMOKE TEST
# sample 0
# --------------------------------------------------

JFTEST=$(sbatch --parsable \
--array=0 \
"$PROJECT/00_metadata/73_dudeml_features.slurm")

echo "feature smoke = $JFTEST"


# --------------------------------------------------
# OTHER 113 SAMPLES
#
# Limit to four concurrent genome-coverage jobs.
# --------------------------------------------------

JFEAT=$(sbatch --parsable \
--dependency=afterok:$JFTEST \
--array=1-113%4 \
"$PROJECT/00_metadata/73_dudeml_features.slurm")

echo "features = $JFEAT"


# --------------------------------------------------
# MODEL SMOKE TEST
# --------------------------------------------------

JMTEST=$(sbatch --parsable \
--dependency=afterok:$JFTEST:$JFEAT \
--array=0 \
"$PROJECT/00_metadata/75_dudeml_train_models.slurm")

echo "model smoke = $JMTEST"


# --------------------------------------------------
# REMAINING 239 MODELS
# --------------------------------------------------

JMODEL=$(sbatch --parsable \
--dependency=afterok:$JMTEST \
--array=1-239%12 \
"$PROJECT/00_metadata/75_dudeml_train_models.slurm")

echo "models = $JMODEL"


# --------------------------------------------------
# PREDICTION SMOKE TEST
# --------------------------------------------------

JPTEST=$(sbatch --parsable \
--dependency=afterok:$JMODEL \
--array=0 \
"$PROJECT/00_metadata/76_dudeml_predict.slurm")

echo "prediction smoke = $JPTEST"


# --------------------------------------------------
# REMAINING 683 PREDICTIONS
# --------------------------------------------------

JPRED=$(sbatch --parsable \
--dependency=afterok:$JPTEST \
--array=1-683%16 \
"$PROJECT/00_metadata/76_dudeml_predict.slurm")

echo "predictions = $JPRED"


# --------------------------------------------------
# AUDIT
# --------------------------------------------------

JAUDIT=$(sbatch --parsable \
--dependency=afterany:$JPRED \
"$PROJECT/00_metadata/77_dudeml_audit.slurm")

echo "audit = $JAUDIT"


cat > "$DUDE/metadata/submitted_jobs.tsv" <<TABLE
component	job_id
feature_smoke	$JFTEST
features	$JFEAT
model_smoke	$JMTEST
models	$JMODEL
prediction_smoke	$JPTEST
predictions	$JPRED
audit	$JAUDIT
TABLE


echo
echo "======================================"
echo "DUDEML PIPELINE SUBMITTED"
echo "======================================"

cat "$DUDE/metadata/submitted_jobs.tsv"
