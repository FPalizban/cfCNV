#!/bin/bash
#SBATCH --job-name=cfCNV_depth
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=30:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/depth_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/depth_%A_%a.err

set -euo pipefail

source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/config.sh
module load SAMtools/1.16.1-GCC-11.3.0

MANIFEST=$1

LINE=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+2)) 'NR==n' "$MANIFEST")

IFS=$'\t' read -r \
SAMPLE REP TARGET SOURCE_DEPTH SOURCE OUTPUT SEED TRUTH <<< "$LINE"

mkdir -p "$(dirname "$OUTPUT")"

echo "Sample:       $SAMPLE"
echo "Replicate:    $REP"
echo "Target depth: ${TARGET}x"
echo "Source depth: ${SOURCE_DEPTH}x"
echo "Input:        $SOURCE"
echo "Output:       $OUTPUT"

if [[ -s "$OUTPUT" ]] && samtools quickcheck "$OUTPUT" 2>/dev/null
then
    echo "Already complete: $OUTPUT"
    exit 0
fi

# If requested depth equals source depth, link instead of duplicating data.
SAME=$(python3 - <<PY
a=float("$TARGET")
b=float("$SOURCE_DEPTH")
print(1 if abs(a-b)<1e-10 else 0)
PY
)

if [[ "$SAME" == "1" ]]
then
    ln -sfn "$SOURCE" "$OUTPUT"

    if [[ -f "${SOURCE}.bai" ]]
    then
        ln -sfn "${SOURCE}.bai" "${OUTPUT}.bai"
    else
        samtools index "$SOURCE"
        ln -sfn "${SOURCE}.bai" "${OUTPUT}.bai"
    fi

    echo "Full-depth file linked rather than copied."
    exit 0
fi

FRACTION=$(python3 - <<PY
source=float("$SOURCE_DEPTH")
target=float("$TARGET")
f=target/source

if not (0 < f < 1):
    raise SystemExit(
        f"Invalid downsampling fraction: {f}"
    )

print(f)
PY
)

SAMARG=$(python3 - <<PY
seed=int("$SEED")
f=float("$FRACTION")
digits=f"{f:.12f}".split(".")[1].rstrip("0")

if not digits:
    digits="0"

print(f"{seed}.{digits}")
PY
)

echo "Fraction: $FRACTION"
echo "samtools -s: $SAMARG"

/usr/bin/time -v \
samtools view \
    -@ ${SLURM_CPUS_PER_TASK:-8} \
    -b \
    -s "$SAMARG" \
    "$SOURCE" \
    -o "$OUTPUT" \
    2> "${OUTPUT}.time.txt"

samtools index \
    -@ ${SLURM_CPUS_PER_TASK:-8} \
    "$OUTPUT"

samtools quickcheck -v "$OUTPUT"

echo "Completed:"
ls -lh "$OUTPUT" "${OUTPUT}.bai"
