#!/bin/bash
#SBATCH --job-name=cfCNV_spike
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=30:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/spike_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/spike_%A_%a.err

set -euo pipefail

source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/config.sh

module load SAMtools/1.16.1-GCC-11.3.0
module load BEDTools/2.30.0-GCC-11.3.0

MANIFEST=$SIM/healthy_cfDNA_spikein/manifests/spike_tasks.tsv

LINE=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+2)) 'NR==n' "$MANIFEST")

IFS=$'\t' read -r \
SAMPLE REP PARENT_DEPTH BAM TRUTH OUTPUT SEED <<< "$LINE"

mkdir -p "$(dirname "$OUTPUT")"

echo "Sample:    $SAMPLE"
echo "Replicate: $REP"
echo "Input:     $BAM"
echo "Truth:     $TRUTH"
echo "Output:    $OUTPUT"

if [[ -s "$OUTPUT" ]] && samtools quickcheck "$OUTPUT" 2>/dev/null
then
    echo "Already complete: $OUTPUT"
    exit 0
fi

TMPROOT=$SIM/healthy_cfDNA_spikein/tmp
mkdir -p "$TMPROOT"

TMP=$(mktemp -d "$TMPROOT/${SAMPLE}_rep${REP}.XXXXXX")

cleanup() {
    rm -rf "$TMP"
}
trap cleanup EXIT

# ------------------------------------------------------------
# Build sorted BED files
# ------------------------------------------------------------

awk 'BEGIN{OFS="\t"} !/^#/ && $4=="loss" {print $1,$2,$3}' \
"$TRUTH" > "$TMP/loss.raw.bed"

awk 'BEGIN{OFS="\t"} !/^#/ && $4=="gain" {print $1,$2,$3}' \
"$TRUTH" > "$TMP/gain.raw.bed"

awk 'BEGIN{OFS="\t"} !/^#/ {print $1,$2,$3}' \
"$TRUTH" > "$TMP/all.raw.bed"

bedtools sort \
-g "$REF37_DIR/Shendure_GRCh37.chrom.sizes" \
-i "$TMP/loss.raw.bed" \
> "$TMP/loss.bed"

bedtools sort \
-g "$REF37_DIR/Shendure_GRCh37.chrom.sizes" \
-i "$TMP/gain.raw.bed" \
> "$TMP/gain.bed"

bedtools sort \
-g "$REF37_DIR/Shendure_GRCh37.chrom.sizes" \
-i "$TMP/all.raw.bed" \
> "$TMP/all.bed"

bedtools complement \
-i "$TMP/all.bed" \
-g "$REF37_DIR/Shendure_GRCh37.chrom.sizes" \
> "$TMP/neutral.bed"

# ------------------------------------------------------------
# Neutral genome
# ------------------------------------------------------------

samtools view \
    -@ ${SLURM_CPUS_PER_TASK:-8} \
    -b \
    -L "$TMP/neutral.bed" \
    "$BAM" \
    -o "$TMP/neutral.bam"

# ------------------------------------------------------------
# LOSS: CN2 -> CN1
# retain ~50% of templates
# ------------------------------------------------------------

samtools view \
    -@ ${SLURM_CPUS_PER_TASK:-8} \
    -b \
    -L "$TMP/loss.bed" \
    "$BAM" \
    -o "$TMP/loss.original.bam"

samtools view \
    -@ ${SLURM_CPUS_PER_TASK:-8} \
    -b \
    -s "${SEED}.500000" \
    "$TMP/loss.original.bam" \
    -o "$TMP/loss.half.bam"

# ------------------------------------------------------------
# GAIN: CN2 -> CN3
# original reads + another ~50%
# ------------------------------------------------------------

samtools view \
    -@ ${SLURM_CPUS_PER_TASK:-8} \
    -b \
    -L "$TMP/gain.bed" \
    "$BAM" \
    -o "$TMP/gain.original.bam"

samtools view \
    -h \
    -s "$((SEED+100)).500000" \
    "$TMP/gain.original.bam" \
| awk '
BEGIN {OFS="\t"}
/^@/ {print; next}
{$1=$1"_SIMGAIN"; print}
' \
| samtools view -b - \
> "$TMP/gain.extra.bam"

# ------------------------------------------------------------
# Reconstruct simulated BAM
# ------------------------------------------------------------

samtools merge \
    -@ ${SLURM_CPUS_PER_TASK:-8} \
    -f \
    "$TMP/merged.bam" \
    "$TMP/neutral.bam" \
    "$TMP/loss.half.bam" \
    "$TMP/gain.original.bam" \
    "$TMP/gain.extra.bam"

samtools sort \
    -@ ${SLURM_CPUS_PER_TASK:-8} \
    -o "$OUTPUT" \
    "$TMP/merged.bam"

samtools index \
    -@ ${SLURM_CPUS_PER_TASK:-8} \
    "$OUTPUT"

samtools quickcheck -v "$OUTPUT"

echo "Finished:"
ls -lh "$OUTPUT" "${OUTPUT}.bai"
