#!/bin/bash
#SBATCH --job-name=S2_CNVref
#SBATCH --cpus-per-task=14
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/CNVkit/reference_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/CNVkit/reference_%j.err

set -euo pipefail

P=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project
cd "$P"

if [[ -f "$P/00_metadata/config.sh" ]]; then
    source "$P/00_metadata/config.sh"
fi

command -v cnvkit.py
command -v Rscript

REF="$P/02_reference/hg38_GIABv3/GRCh38_GIABv3.fa"
OD="$P/07_tool_outputs/simulation2/CNVkit/reference"

ACCESS="$OD/access_autosomes.bed"
BINS="$OD/targets_1000000.bed"
CNN="$OD/flat_GRCh38_1000000.cnn"

if [[ -s "$CNN" ]]; then
    echo "REFERENCE_ALREADY_READY"
    exit 0
fi

cnvkit.py access "$REF" -o "$OD/access_all.bed"

awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/' \
    "$OD/access_all.bed" > "$ACCESS"

cnvkit.py target "$ACCESS" \
    --split --avg-size 1000000 \
    -o "$BINS"

cnvkit.py reference \
    -t "$BINS" \
    -f "$REF" \
    -o "$CNN"

test -s "$CNN"

echo "CNVKIT_REFERENCE_READY"
