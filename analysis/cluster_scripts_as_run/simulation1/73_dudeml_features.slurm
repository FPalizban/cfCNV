#!/bin/bash
#SBATCH --job-name=dude_feat
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=36:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/dudeML/dude_feat_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/dudeML/dude_feat_%A_%a.err

set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project

source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

module load BEDTools/2.30.0-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0

TOOLS="$PROJECT/06_tools"
RUN="$TOOLS/bin/run_dudeml"

DUDE="$PROJECT/07_tool_outputs/simulation1/dudeML"

MASTER="$SIM1_MAN/master.tsv"

LINE=$(awk -F'\t' \
-v n=$((SLURM_ARRAY_TASK_ID+2)) \
'NR==n {print; exit}' \
"$MASTER")

[[ -n "$LINE" ]] || {
    echo "ERROR: no manifest row"
    exit 2
}

IFS=$'\t' read -r \
SAMPLE BAM BACKGROUND DEPTH REP CLASS \
<<< "$LINE"

OD="$DUDE/features/$SAMPLE"

mkdir -p "$OD"

if [[ -f "$OD/.features_done" ]]; then
    echo "$SAMPLE already complete"
    exit 0
fi

echo "======================================"
echo "dudeML feature generation"
echo "sample     = $SAMPLE"
echo "background = $BACKGROUND"
echo "depth      = $DEPTH"
echo "replicate  = $REP"
echo "class      = $CLASS"
echo "bam        = $BAM"
echo "======================================"

[[ -s "$BAM" ]] || {
    echo "ERROR BAM missing"
    exit 3
}


# ---------------------------------------------
# ONE genome-wide coverage pass,
# SIX output bin sizes
# ---------------------------------------------

NEED_WINSTAT=0

for BIN in \
10000 30000 50000 100000 500000 1000000
do
    WIN="$OD/${SAMPLE}.${BIN}.winStat.bed"

    if [[ ! -s "$WIN" ]]
    then
        NEED_WINSTAT=1
    fi
done

if [[ "$NEED_WINSTAT" -eq 1 ]]
then

    echo
    echo "Generating genome-wide window statistics..."

    python \
    "$PROJECT/00_metadata/72_dudeml_winstat_bga.py" \
    --bam "$BAM" \
    --sample "$SAMPLE" \
    --outdir "$OD" \
    --chroms "$DUDE/metadata/chromosomes.txt" \
    --bins \
    10000 30000 50000 100000 500000 1000000

else

    echo
    echo "All six winStat files already exist."
    echo "Reusing existing window statistics."

fi


# ---------------------------------------------
# Generate dudeML feature vectors
# ---------------------------------------------

for BIN in \
10000 30000 50000 100000 500000 1000000
do

    BD="$OD/bin_${BIN}"

    mkdir -p "$BD"

    WIN="$OD/${SAMPLE}.${BIN}.winStat.bed"

    [[ -s "$WIN" ]] || {
        echo "ERROR missing $WIN"
        exit 10
    }

    cd "$BD"


    # TEST/SAMPLE FEATURE VECTOR
    "$RUN" fvecSample \
        -i "$WIN" \
        -o "$BD/${SAMPLE}.sample.bed" \
        -id "$SAMPLE" \
        -d "$BD" \
        -w "$BIN" \
        -s "$BIN" \
        -windows 5


    [[ -s "$BD/${SAMPLE}.sample.bed" ]] || {
        echo "ERROR sample feature vector missing"
        exit 11
    }


    # TRAINING FEATURE VECTOR
    # Only positive samples have known CNVs.
    if [[ "$CLASS" == "positive" ]]
    then

        DUP="$DUDE/metadata/labels/rep${REP}.dup.bed"
        DEL="$DUDE/metadata/labels/rep${REP}.del.bed"

        "$RUN" fvecTrain \
            -i "$WIN" \
            -o "$BD/${SAMPLE}.train.bed" \
            -d "$BD" \
            -w "$BIN" \
            -s "$BIN" \
            -windows 5 \
            -dups "$DUP" \
            -dels "$DEL"

        [[ -s "$BD/${SAMPLE}.train.bed" ]] || {
            echo "ERROR training feature vector missing"
            exit 12
        }

        echo "Training classes:"
        cut -f4 "$BD/${SAMPLE}.train.bed" \
        | sort \
        | uniq -c

    fi

done

touch "$OD/.features_done"

echo
echo "DUDEML FEATURES COMPLETE: $SAMPLE"
