#!/bin/bash
#SBATCH --job-name=S2_dude_feat
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/dude_feat_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/dude_feat_%A_%a.err
# A6a  dudeML window statistics + feature vectors, the same scripts and settings as
#      Simulation 1 (6 bins, 5 windows), one task per mixture (104 tasks).
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/simulation2/full/s2_env.sh
module load BEDTools/2.30.0-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0
IFS=$'\t' read -r S BAM D X < <(task_line "$S2_META/tasks_dudeml_feat.tsv" "$SLURM_ARRAY_TASK_ID")
OD="$S2_OUT/dudeML/features/$S"; mkdir -p "$OD"
[[ -f "$OD/.features_done" ]] && { echo "SKIP $S"; exit 0; }
python3 "$P/00_metadata/72_dudeml_winstat_bga.py" --bam "$BAM" --sample "$S" --outdir "$OD" \
   --chroms "$S2_META/chromosomes_hg38.txt" --bins $BINS
for B in $BINS; do
  BD="$OD/bin_$B"; mkdir -p "$BD"; cd "$BD"
  "$TOOLS/bin/run_dudeml" fvecSample -i "$OD/$S.$B.winStat.bed" -o "$BD/$S.sample.bed" \
      -id "$S" -d "$BD" -w "$B" -s "$B" -windows 5
  test -s "$BD/$S.sample.bed"
done
touch "$OD/.features_done"; echo "DUDE_FEATURES_DONE $S"
