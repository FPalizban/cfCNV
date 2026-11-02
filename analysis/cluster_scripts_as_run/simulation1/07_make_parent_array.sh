#!/bin/bash
#SBATCH --job-name=cfCNV_parent
#SBATCH --cpus-per-task=25
#SBATCH --mem=180G
#SBATCH --time=30:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/parent_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/parent_%A_%a.err

set -euo pipefail

source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/config.sh
module load SAMtools/1.16.1-GCC-11.3.0

MANIFEST=$SIM/healthy_cfDNA_spikein/manifests/parent_tasks.tsv

LINE=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+2)) 'NR==n' "$MANIFEST")

IFS=$'\t' read -r SAMPLE NATIVE PARENT SOURCE OUTPUT SEED <<< "$LINE"

mkdir -p "$(dirname "$OUTPUT")"

echo "Sample:       $SAMPLE"
echo "Native depth: $NATIVE"
echo "Target depth: $PARENT"
echo "Source:       $SOURCE"
echo "Output:       $OUTPUT"

if [[ -s "$OUTPUT" ]] && samtools quickcheck "$OUTPUT" 2>/dev/null
then
    echo "Already complete: $OUTPUT"
    exit 0
fi

FRACTION=$(python3 - <<PY
print(float("$PARENT") / float("$NATIVE"))
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

echo "Sampling fraction: $FRACTION"
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
