#!/bin/bash
# =====================================================================
# 5_prepare_hf_release.sh - builds the folder to upload to Hugging Face (login node):  bash 5_prepare_hf_release.sh
# Copies (never moves or deletes) into $P/13_ml_dataset/hf_release :
#   README.md                 the dataset card (from hf_card/README.md of this package)
#   data/bins_100kb, bins_1Mb the six splits at both resolutions
#   data/augmented_100kb      the augmented training / validation samples of model v2
#   samples.tsv               one row per sample
#   baselines/                result tables of the baseline models + run settings
# An older hf_release folder is renamed, not removed.
# =====================================================================
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/ml_dataset/ml_env.sh
REL="$MLO/hf_release"
for f in "$MLO/hf_dataset/data/bins_100kb/train.csv.gz" "$MLO/augmented/train_augmented_100kb.csv.gz" "$MLO/manuscript/key_numbers.md"; do
  [[ -f "$f" ]] || { echo "missing $f - run 4_final.sh first"; exit 1; }
done
[[ -e "$REL" ]] && mv "$REL" "$MLO/hf_release.prev_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$REL/data/augmented_100kb" "$REL/baselines"
cp -r "$MLO/hf_dataset/data/bins_100kb" "$MLO/hf_dataset/data/bins_1Mb" "$REL/data/"
cp "$MLO/augmented/train_augmented_100kb.csv.gz"      "$REL/data/augmented_100kb/train.csv.gz"
cp "$MLO/augmented/validation_augmented_100kb.csv.gz" "$REL/data/augmented_100kb/validation.csv.gz"
cp "$MLO/hf_dataset/samples.tsv" "$REL/"
cp "$ML/hf_card/README.md" "$REL/README.md"
cp "$MLO"/manuscript/table_*.tsv "$MLO/manuscript/hg008_profile_r_per_sample.tsv" "$MLO/manuscript/hg008_lod_by_depth.tsv" "$REL/baselines/"
cp "$MLO/results/run_info.json"    "$REL/baselines/run_info_v1.json"
cp "$MLO/results_v2/run_info.json" "$REL/baselines/run_info_v2.json"
chmod -R u+rwX,go-w "$REL"
echo "== release folder: $REL"; du -sh "$REL"; find "$REL" -type f | sort | sed "s#$REL/##"
echo "== samples per split:"; awk -F'\t' 'NR>1{n[$2]++} END{for(k in n) print "   "k": "n[k]}' "$REL/samples.tsv" | sort
cd "$MLO" && tar cf hf_release.tar hf_release && echo "== one file to copy to your laptop: $MLO/hf_release.tar" && ls -la hf_release.tar
