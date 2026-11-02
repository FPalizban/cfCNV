#!/bin/bash
#SBATCH --job-name=cfCNV_down
#SBATCH --cpus-per-task=25
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/downsample_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/downsample_%j.err

set -euo pipefail

source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/config.sh

module load SAMtools/1.16.1-GCC-11.3.0

BAM=$1
SOURCE_DEPTH=$2
TARGET_DEPTH=$3
OUTPUT=$4
SEED=${5:-42}

mkdir -p "$(dirname "$OUTPUT")"

FRACTION=$(python3 - <<PY
source=float("$SOURCE_DEPTH")
target=float("$TARGET_DEPTH")
print(target/source)
PY
)

echo "Source BAM:     $BAM"
echo "Source depth:   $SOURCE_DEPTH"
echo "Target depth:   $TARGET_DEPTH"
echo "Fraction:       $FRACTION"
echo "Output:         $OUTPUT"

VALID=$(python3 - <<PY
f=float("$FRACTION")
print(1 if 0 < f < 1 else 0)
PY
)

if [ "$VALID" != "1" ]; then
    echo "ERROR: sampling fraction must be between 0 and 1."
    exit 1
fi

SAMARG=$(python3 - <<PY
seed=int("$SEED")
f=float("$FRACTION")
digits=f"{f:.12f}".split(".")[1].rstrip("0")
if not digits:
    digits="0"
print(f"{seed}.{digits}")
PY
)

echo "samtools -s argument: $SAMARG"

/usr/bin/time -v \
samtools view \
    -@ ${SLURM_CPUS_PER_TASK:-8} \
    -b \
    -s "$SAMARG" \
    "$BAM" \
    -o "$OUTPUT" \
    2> "${OUTPUT}.time.txt"

samtools index \
    -@ ${SLURM_CPUS_PER_TASK:-8} \
    "$OUTPUT"

samtools quickcheck -v "$OUTPUT"

echo "Finished."
ls -lh "$OUTPUT" "${OUTPUT}.bai"
