#!/bin/bash
#SBATCH --job-name=S2_meann_feat
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/meann_feat_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/meann_feat_%A_%a.err
# A7a  CNV-MEANN raw features (same script as Simulation 1), then the v2 depth/GC
#      normalisation. One task per mixture x bin (624 tasks).
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/simulation2/full/s2_env.sh
IFS=$'\t' read -r S BAM D B < <(task_line "$S2_META/tasks_meann_feat.tsv" "$SLURM_ARRAY_TASK_ID")
RAW="$S2_OUT/CNV-MEANN/features_raw/bin_$B/$S"; NRM="$S2_OUT/CNV-MEANN/features_v2/bin_$B/$S"
[[ -f "$NRM/.done" ]] && { echo "SKIP $S $B"; exit 0; }
mkdir -p "$RAW"
[[ -f "$RAW/.done.bin$B" ]] || pyd "$P/00_metadata/cfcnv_ml/90_cnvmeann_features.py" \
     --bam "$BAM" --ref "$REF38" --sample "$S" --outdir "$RAW" --bin "$B"
pyd "$S2/a7_meann_normalize.py" --raw "$RAW" --out "$NRM" --sample "$S" --bin "$B"
echo "MEANN_FEATURES_DONE $S $B"
