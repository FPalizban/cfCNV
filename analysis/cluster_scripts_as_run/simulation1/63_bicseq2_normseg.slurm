#!/bin/bash
#SBATCH --job-name=bic_normseg
#SBATCH --cpus-per-task=20
#SBATCH --mem=160G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/bic_normseg_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/bic_normseg_%A_%a.err

set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project

source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

TOOLS="$PROJECT/06_tools"
WRAP="$TOOLS/bin/run_bicseq2"

SIM1_OUT="${SIM1_OUT:-$PROJECT/07_tool_outputs/simulation1}"

BICOUT="$SIM1_OUT/BICseq2"
BICRES="$TOOLS/BICseq2_resources/GRCh37_101mer"

CHRFILE="$BICRES/chromosomes.txt"
REFDIR="$BICRES/reference"

MER="human_g1k_v37.101mer"
MAPDIR="$BICRES/mappability/$MER"

TASKS="$SIM1_MAN/bicseq2_tasks.tsv"

LINE=$(awk -F'\t' \
-v n=$((SLURM_ARRAY_TASK_ID+2)) \
'NR==n {print; exit}' \
"$TASKS")

[[ -n "$LINE" ]] || {
    echo "ERROR: task $SLURM_ARRAY_TASK_ID not found"
    exit 2
}

IFS=$'\t' read -r \
TASK SAMPLE BAM BACKGROUND DEPTH REPLICATE CLASS BIN \
<<< "$LINE"

OD="$BICOUT/$SAMPLE/bin_${BIN}"
NORMD="$OD/norm"
SEGD="$OD/seg"
TMPD="$OD/tmp"

mkdir -p "$NORMD" "$SEGD" "$TMPD"

echo "=========================================="
echo "BIC-seq2 Simulation-1"
echo "=========================================="
echo "task       = $TASK"
echo "sample     = $SAMPLE"
echo "background = $BACKGROUND"
echo "depth      = $DEPTH"
echo "replicate  = $REPLICATE"
echo "class      = $CLASS"
echo "bin bp     = $BIN"
echo "=========================================="

if [[ -f "$OD/.done" ]]; then
    echo "Already complete"
    exit 0
fi

[[ -f "$BICRES/RESOURCES_READY" ]] || {
    echo "ERROR: resources are not ready"
    exit 3
}

[[ -f "$BICOUT/unique/$SAMPLE/.done" ]] || {
    echo "ERROR: unique-read data missing for $SAMPLE"
    exit 4
}

FRAGFILE="$BICOUT/unique/$SAMPLE/fragment_size.txt"

[[ -s "$FRAGFILE" ]] || {
    echo "ERROR: fragment-size file missing for $SAMPLE"
    exit 5
}

FRAG=$(cat "$FRAGFILE")

echo "fragment size = $FRAG"


# ==========================================
# Normalization configuration
# ==========================================

NCFG="$OD/${SAMPLE}.norm-config.tsv"

printf \
"chromName\tfaFile\tMapFile\treadPosFile\tbinFileNorm\n" \
> "$NCFG"

while read -r CHR
do

    FA="$REFDIR/${CHR}.fa"
    MAP="$MAPDIR/${MER}.${CHR}.txt"
    SEQ="$BICOUT/unique/$SAMPLE/${SAMPLE}.${CHR}.seq"
    NORM="$NORMD/${SAMPLE}.${CHR}.norm.bin"

    [[ -s "$FA" ]] || {
        echo "ERROR: missing reference $FA"
        exit 10
    }

    [[ -s "$MAP" ]] || {
        echo "ERROR: missing mappability $MAP"
        exit 11
    }

    [[ -f "$SEQ" ]] || {
        echo "ERROR: missing read-position file $SEQ"
        exit 12
    }

    printf \
    "%s\t%s\t%s\t%s\t%s\n" \
    "$CHR" \
    "$FA" \
    "$MAP" \
    "$SEQ" \
    "$NORM" \
    >> "$NCFG"

done < "$CHRFILE"


# ==========================================
# Normalization
# ==========================================

PDF="$OD/${SAMPLE}.GC.pdf"
PARS="$OD/${SAMPLE}.norm_parameters.txt"

echo
echo "===== NORMALIZATION ====="

"$WRAP" \
perl /NBICseq-norm_v0.2.4/NBICseq-norm.pl \
--tmp="$TMPD" \
-l 101 \
-s "$FRAG" \
-b "$BIN" \
--fig "$PDF" \
"$NCFG" \
"$PARS"


while read -r CHR
do

    F="$NORMD/${SAMPLE}.${CHR}.norm.bin"

    [[ -f "$F" ]] || {
        echo "ERROR: normalization output missing: $F"
        exit 20
    }

done < "$CHRFILE"


# ==========================================
# Case-only segmentation configuration
# ==========================================

SCFG="$OD/${SAMPLE}.seg-config.tsv"

printf \
"chromName\tbinFileNorm\n" \
> "$SCFG"

while read -r CHR
do

    printf \
    "%s\t%s\n" \
    "$CHR" \
    "$NORMD/${SAMPLE}.${CHR}.norm.bin" \
    >> "$SCFG"

done < "$CHRFILE"


mkdir -p "$SEGD/tmp"

CNV="$SEGD/${SAMPLE}.cnv"
FIG="$SEGD/${SAMPLE}.png"


# ==========================================
# Case-only segmentation
# lambda=2 is BIC-seq2 default.
# No --control because this is case-only.
# ==========================================

echo
echo "===== SEGMENTATION ====="

"$WRAP" \
perl /NBICseq-seg_v0.7.2/NBICseq-seg.pl \
--lambda=2 \
--tmp="$SEGD/tmp" \
--fig="$FIG" \
"$SCFG" \
"$CNV"


# A zero-call sample must NOT be classified as a
# pipeline failure merely because it has no CNV calls.
[[ -f "$CNV" ]] || {
    echo "ERROR: segmentation output file missing"
    exit 30
}

touch "$OD/.done"

echo
echo "=========================================="
echo "BICSEQ2 TASK COMPLETE"
echo "CNV: $CNV"
echo "=========================================="
