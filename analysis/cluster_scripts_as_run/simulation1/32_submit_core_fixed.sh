#!/bin/bash

set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project

source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

echo "============================================"
echo "cfCNV SIMULATION 1 — CLEAN SUBMISSION"
echo "============================================"


# --------------------------------------------------
# Validate manifests
# --------------------------------------------------

MASTER_N=$(($(wc -l < "$SIM1_MAN/master.tsv") - 1))
Q_N=$(($(wc -l < "$SIM1_MAN/qdnaseq_tasks.tsv") - 1))
C_N=$(($(wc -l < "$SIM1_MAN/cnvkit_tasks.tsv") - 1))
I_N=$(($(wc -l < "$SIM1_MAN/ichor_tasks.tsv") - 1))
S_N=$(($(wc -l < "$SIM1_MAN/spectre_tasks.tsv") - 1))

echo "Master  = $MASTER_N"
echo "QDNAseq = $Q_N"
echo "CNVkit  = $C_N"
echo "ichor   = $I_N"
echo "Spectre = $S_N"

[[ "$MASTER_N" -eq 114 ]]
[[ "$Q_N" -eq 684 ]]
[[ "$C_N" -eq 684 ]]
[[ "$I_N" -eq 114 ]]
[[ "$S_N" -eq 114 ]]


# --------------------------------------------------
# Ensure CNVkit tumor-only WGS uses flat reference
# --------------------------------------------------

CNVK="$PROJECT/00_metadata/11_cnvkit_array.sh"

if grep -q -- '-m wgs -f' "$CNVK"; then
    echo "Adding CNVkit flat-reference option (-n)..."
    sed -i 's/-m wgs -f/-m wgs -n -f/' "$CNVK"
fi

bash -n "$CNVK"


# --------------------------------------------------
# Submit reference
# --------------------------------------------------

JREF=$(sbatch --parsable \
"$PROJECT/00_metadata/02_prepare_reference.sh")

echo "Reference = $JREF"


# --------------------------------------------------
# Submit caller arrays
# --------------------------------------------------

JQ=$(sbatch --parsable \
--dependency=afterok:$JREF \
--array=0-$((Q_N-1))%18 \
"$PROJECT/00_metadata/10_qdnaseq_array.sh")

JC=$(sbatch --parsable \
--dependency=afterok:$JREF \
--array=0-$((C_N-1))%18 \
"$PROJECT/00_metadata/11_cnvkit_array.sh")

JI=$(sbatch --parsable \
--dependency=afterok:$JREF \
--array=0-$((I_N-1))%12 \
"$PROJECT/00_metadata/12_ichorcna_array.sh")

JS=$(sbatch --parsable \
--dependency=afterok:$JREF \
--array=0-$((S_N-1))%12 \
"$PROJECT/00_metadata/13_spectre_array.sh")


# --------------------------------------------------
# Record IDs
# --------------------------------------------------

STATUS="$SIM1_MAN/submitted_jobs.tsv"

cat > "$STATUS" <<TABLE
component	job_id
reference	$JREF
QDNAseq	$JQ
CNVkit	$JC
ichorCNA	$JI
Spectre	$JS
TABLE


echo
echo "============================================"
echo "SUBMITTED SUCCESSFULLY"
echo "============================================"
cat "$STATUS"

echo
echo "NOTE: scoring intentionally NOT submitted yet."
echo "We will validate outputs and then run corrected scoring."
echo "============================================"
