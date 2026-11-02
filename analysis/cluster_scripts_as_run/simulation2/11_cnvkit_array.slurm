#!/bin/bash
#SBATCH --job-name=S2_CNVkit
#SBATCH --cpus-per-task=15
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/CNVkit/sample_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/CNVkit/sample_%A_%a.err

set -euo pipefail

P=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project
cd "$P"

if [[ -f "$P/00_metadata/config.sh" ]]; then
    source "$P/00_metadata/config.sh"
fi

command -v cnvkit.py

MAN="$P/00_metadata/simulation2/HG008_mixtures.tsv"
MIX="$P/05_simulation/simulation2_HG008/mixtures"

CNN="$P/07_tool_outputs/simulation2/CNVkit/reference/flat_GRCh38_1000000.cnn"

LINE=$(awk -F '\t' -v n="$((SLURM_ARRAY_TASK_ID+2))" \
    'NR==n {print}' "$MAN")

IFS=$'\t' read -r ID DEPTH FRACTION <<< "$LINE"

test -n "$ID"

BAM="$MIX/$ID.bam"
OD="$P/07_tool_outputs/simulation2/CNVkit/$ID"

mkdir -p "$OD"

test -s "$BAM"
test -s "$CNN"

if [[ -s "$OD/.done" && -s "$OD/$ID.cns" ]]; then
    echo "SKIP_COMPLETED $ID"
    exit 0
fi

/usr/bin/time -v cnvkit.py batch "$BAM" \
    -m wgs \
    -r "$CNN" \
    -d "$OD" \
    -p "${SLURM_CPUS_PER_TASK}"

test -s "$OD/$ID.cns"

# Preserve all segment log2 values for subsequent threshold analysis.
awk -v sample="$ID" '
BEGIN {
    OFS="\t"
    print "sample_id","chrom","start","end","log2"
}
NR>1 {
    print sample,$1,$2,$3,$5
}' "$OD/$ID.cns" > "$OD/$ID.segments.tsv"

test -s "$OD/$ID.segments.tsv"

printf '%s\t%s\t%s\n' \
    "$ID" "$DEPTH" "$FRACTION" > "$OD/.done"

echo "CNVKIT_SAMPLE_COMPLETE $ID"
