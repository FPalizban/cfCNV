#!/bin/bash
#SBATCH --job-name=bic_unique
#SBATCH --array=0-113%8
#SBATCH --cpus-per-task=20
#SBATCH --mem=160G
#SBATCH --time=18:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/bic_unique_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/bic_unique_%A_%a.err

set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project

source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

TOOLS="$PROJECT/06_tools"
WRAP="$TOOLS/bin/run_bicseq2"

SIM1_OUT="${SIM1_OUT:-$PROJECT/07_tool_outputs/simulation1}"
BICOUT="$SIM1_OUT/BICseq2"

CHRFILE="$TOOLS/BICseq2_resources/GRCh37_101mer/chromosomes.txt"

LINE=$(awk -F'\t' \
-v n=$((SLURM_ARRAY_TASK_ID+2)) \
'NR==n {print; exit}' \
"$SIM1_MAN/master.tsv")

[[ -n "$LINE" ]] || exit 2

IFS=$'\t' read -r \
SAMPLE BAM BACKGROUND DEPTH REPLICATE CLASS \
<<< "$LINE"

OUT="$BICOUT/unique/$SAMPLE"

mkdir -p "$OUT"

echo "===================================="
echo "sample = $SAMPLE"
echo "depth  = $DEPTH"
echo "bam    = $BAM"
echo "===================================="

[[ -s "$BAM" ]] || {
    echo "ERROR: BAM missing"
    exit 3
}

if [[ -f "$OUT/.done" ]]; then
    echo "Already completed"
    exit 0
fi


"$WRAP" bash -c '
set -euo pipefail

BAM="$1"
SAMPLE="$2"
OUT="$3"
CHRFILE="$4"

HELPER=/samtools-0.1.7a_getUnique-0.1.3/misc/samtools.pl

[[ -f "$HELPER" ]] || {
    echo "ERROR: getUnique helper not found"
    exit 10
}

PIDS=()

wait_group () {
    for P in "${PIDS[@]}"
    do
        wait "$P"
    done
    PIDS=()
}

while read -r CHR
do

    (
        set -o pipefail

        samtools view "$BAM" "$CHR" \
        | perl "$HELPER" unique - \
        | cut -f4 \
        > "$OUT/${SAMPLE}.${CHR}.seq"
    ) &

    PIDS+=("$!")

    if [[ "${#PIDS[@]}" -eq 4 ]]; then
        wait_group
    fi

done < "$CHRFILE"

if [[ "${#PIDS[@]}" -gt 0 ]]; then
    wait_group
fi

' _ "$BAM" "$SAMPLE" "$OUT" "$CHRFILE"


# -------------------------------------------------
# VALIDATE 22 OUTPUTS
# -------------------------------------------------

N=0

while read -r CHR
do

    F="$OUT/${SAMPLE}.${CHR}.seq"

    [[ -f "$F" ]] || {
        echo "ERROR: missing $F"
        exit 20
    }

    N=$((N+1))

done < "$CHRFILE"

[[ "$N" -eq 22 ]]


# -------------------------------------------------
# ESTIMATE MEDIAN FRAGMENT SIZE ONCE PER SAMPLE
# -------------------------------------------------

TMP="$OUT/fragment_lengths.tmp"

"$WRAP" bash -c '
set +o pipefail

samtools view \
-f 2 \
-F 3852 \
"$1" \
| awk '"'"'
$9>0 && $9<2000 {
    print $9
    n++
    if(n>=100000)
        exit
}
'"'"'
' _ "$BAM" > "$TMP"

if [[ $(wc -l < "$TMP") -ge 100 ]]; then

    sort -n "$TMP" \
    | awk '
    {
        a[NR]=$1
    }
    END {
        if(NR%2)
            print a[(NR+1)/2]
        else
            print int((a[NR/2]+a[NR/2+1])/2+0.5)
    }
    ' > "$OUT/fragment_size.txt"

else

    echo 167 > "$OUT/fragment_size.txt"

fi

rm -f "$TMP"

echo "fragment size = $(cat "$OUT/fragment_size.txt")"

touch "$OUT/.done"

echo "BICSEQ2 UNIQUE COMPLETE"
