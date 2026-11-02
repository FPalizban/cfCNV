#!/bin/bash
#SBATCH --job-name=cfCNV_spike
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/spike_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/spike_%j.err

set -euo pipefail

source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/config.sh

module load SAMtools/1.16.1-GCC-11.3.0
module load BEDTools/2.30.0-GCC-11.3.0

BAM=$1
TRUTH=$2
OUTPUT=$3
SEED=${4:-42}

TMPROOT=${TMPDIR:-$SIM/tmp}
mkdir -p "$TMPROOT"

TMP=$(mktemp -d "$TMPROOT/cfCNV.XXXXXX")

echo "Input BAM: $BAM"
echo "Truth:     $TRUTH"
echo "Output:    $OUTPUT"
echo "Temp:      $TMP"

mkdir -p "$(dirname "$OUTPUT")"

# Extract truth intervals.
awk 'BEGIN{OFS="\t"} !/^#/ && $4=="loss" {print $1,$2,$3}' \
"$TRUTH" > "$TMP/loss.unsorted.bed"

awk 'BEGIN{OFS="\t"} !/^#/ && $4=="gain" {print $1,$2,$3}' \
"$TRUTH" > "$TMP/gain.unsorted.bed"

awk 'BEGIN{OFS="\t"} !/^#/ {print $1,$2,$3}' \
"$TRUTH" > "$TMP/all.unsorted.bed"

bedtools sort \
-g "$REF37_DIR/Shendure_GRCh37.chrom.sizes" \
-i "$TMP/loss.unsorted.bed" \
> "$TMP/loss.bed"

bedtools sort \
-g "$REF37_DIR/Shendure_GRCh37.chrom.sizes" \
-i "$TMP/gain.unsorted.bed" \
> "$TMP/gain.bed"

bedtools sort \
-g "$REF37_DIR/Shendure_GRCh37.chrom.sizes" \
-i "$TMP/all.unsorted.bed" \
> "$TMP/all_cnv.bed"

# Regions not affected by a simulated CNV.
bedtools complement \
-i "$TMP/all_cnv.bed" \
-g "$REF37_DIR/Shendure_GRCh37.chrom.sizes" \
> "$TMP/neutral.bed"

# Retain neutral genome.
samtools view \
-@ ${SLURM_CPUS_PER_TASK:-8} \
-b \
-L "$TMP/neutral.bed" \
"$BAM" \
-o "$TMP/neutral.bam"

# CN2 -> CN1:
# retain approximately half the reads.
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

# CN2 -> CN3:
# original region + ~50% additional reads.
samtools view \
-@ ${SLURM_CPUS_PER_TASK:-8} \
-b \
-L "$TMP/gain.bed" \
"$BAM" \
-o "$TMP/gain.original.bam"

samtools view \
-h \
-s "$((SEED+1)).500000" \
"$TMP/gain.original.bam" \
| awk '
BEGIN {OFS="\t"}

(/^@/) {
    print
    next
}

{
    $1=$1"_SIMGAIN"
    print
}
' \
| samtools view -b - \
> "$TMP/gain.extra.bam"

# Reassemble genome.
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

rm -rf "$TMP"

echo "Simulation finished."
ls -lh "$OUTPUT" "${OUTPUT}.bai"
