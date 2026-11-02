#!/bin/bash
#SBATCH --job-name=lyc_call
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/hg008/lyc_call_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/hg008/lyc_call_%A_%a.err
#
# LYCEUM, out-of-the-box, exon level, on HG008 tumour-normal mixtures (GRCh38).
# One array task per mixture BAM. Manifest columns (tab, with header):
#   sample_id  bam  tumor_fraction  target_depth
# Each sample is processed in its own folder so the per-sample mean/std
# lookup that LYCEUM uses for normalisation is never computed across samples.

set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project
LYC="$PROJECT/06_tools/ML_based/LYCEUM"
BED="$LYC/hg38_hglft_genome_64dc_dcbaa0_unique.bed"
MAN="${MAN:-$PROJECT/00_metadata/manifests/hg008/hg008_mixtures.tsv}"
OUTROOT="$PROJECT/07_tool_outputs/hg008/LYCEUM"

LINE=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+2)) 'NR==n{print;exit}' "$MAN")
[[ -n "$LINE" ]] || { echo "no manifest line"; exit 2; }
IFS=$'\t' read -r SAMPLE BAM TF DEPTH <<< "$LINE"

OD="$OUTROOT/$SAMPLE"
if [[ -f "$OD/.done" ]]; then echo "already done"; exit 0; fi
rm -rf "$OD"; mkdir -p "$OD"/{read_depths,processed,stats,calls}

# Chromosome-name check: BAM and BED must agree (chr1 vs 1)
BAMCHR=$(samtools view -H "$BAM" | awk '$1=="@SQ"{sub("SN:","",$2); print $2; exit}')
BEDCHR=$(head -1 "$BED" | cut -f1)
if [[ "${BAMCHR:0:3}" != "${BEDCHR:0:3}" ]]; then
    echo "ERROR: chromosome naming differs (BAM=$BAMCHR, BED=$BEDCHR)"; exit 3
fi

eval "$(conda shell.bash hook)"
conda activate lyceum_env
cd "$LYC"   # released scripts expect to run from the repo root

echo "sample=$SAMPLE TF=$TF depth=$DEPTH"
/usr/bin/time -v bash -c "
  sambamba depth base -t 4 -L '$BED' '$BAM' > '$OD/read_depths/$SAMPLE.bam.txt'
  python ./scripts/preprocess_sample.py --readdepth '$OD/read_depths' \
         --output '$OD/processed' --target '$BED'
  python ./scripts/mean_std_calculator.py --processed_samples '$OD/processed' \
         --setname '$SAMPLE' --output '$OD/stats/'
  python ./scripts/LYCEUM_call.py --model lyceum --input '$OD/processed' \
         --output '$OD/calls' --cnv exonlevel --batch_size 16 \
         --stats_lookup '$OD/stats/${SAMPLE}_stats_lookup.npy'
"
ls "$OD/calls" | grep -q . || { echo "ERROR: no LYCEUM output"; exit 4; }
touch "$OD/.done"
echo "LYCEUM DONE $SAMPLE"
