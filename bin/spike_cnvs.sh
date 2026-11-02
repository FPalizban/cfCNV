#!/usr/bin/env bash
# Simulation 1 spike-in, identical to the cfCNV paper (05_spike_CNVs.sh):
#   loss (CN2->CN1): keep ~50% of reads in the event; gain (CN2->CN3): add a ~50% copy of the
#   reads in the event (read names suffixed _SIMGAIN). The rest of the genome is untouched.
# usage: spike_cnvs.sh <in.bam> <truth.tsv> <chrom_sizes> <out.bam> <seed> <threads>
set -euo pipefail
BAM=$1; TRUTH=$2; SIZES=$3; OUT=$4; SEED=${5:-42}; T=${6:-2}
TMP=$(mktemp -d ./spike.XXXXXX)
cut -f1,2 "$SIZES" > "$TMP/g"
awk 'BEGIN{OFS="\t"} !/^#/ && $4=="loss"{print $1,$2,$3}' "$TRUTH" | bedtools sort -g "$TMP/g" -i - > "$TMP/loss.bed"
awk 'BEGIN{OFS="\t"} !/^#/ && $4=="gain"{print $1,$2,$3}' "$TRUTH" | bedtools sort -g "$TMP/g" -i - > "$TMP/gain.bed"
awk 'BEGIN{OFS="\t"} !/^#/{print $1,$2,$3}' "$TRUTH" | bedtools sort -g "$TMP/g" -i - > "$TMP/all.bed"
bedtools complement -i "$TMP/all.bed" -g "$TMP/g" > "$TMP/neutral.bed"
samtools view -@ "$T" -b -L "$TMP/neutral.bed" "$BAM" -o "$TMP/neutral.bam"
samtools view -@ "$T" -b -L "$TMP/loss.bed" "$BAM" -o "$TMP/loss.original.bam"
samtools view -@ "$T" -b -s "${SEED}.500000" "$TMP/loss.original.bam" -o "$TMP/loss.half.bam"
samtools view -@ "$T" -b -L "$TMP/gain.bed" "$BAM" -o "$TMP/gain.original.bam"
samtools view -h -s "$((SEED+1)).500000" "$TMP/gain.original.bam" \
  | awk 'BEGIN{OFS="\t"} /^@/{print; next} {$1=$1"_SIMGAIN"; print}' \
  | samtools view -b - > "$TMP/gain.extra.bam"
samtools merge -@ "$T" -f "$TMP/merged.bam" "$TMP/neutral.bam" "$TMP/loss.half.bam" "$TMP/gain.original.bam" "$TMP/gain.extra.bam"
samtools sort -@ "$T" -o "$OUT" "$TMP/merged.bam"
samtools index -@ "$T" "$OUT"
samtools quickcheck -v "$OUT"
rm -rf "$TMP"
