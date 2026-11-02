#!/bin/bash
#SBATCH --job-name=P3_derive
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/phase3/derive_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/phase3/derive_%A_%a.err
# 12  One sample: measured autosomal depth (mosdepth), then the analysis BAMs:
#     full (as aligned, <= ~12x), and down-sampled copies at 1x and 0.1x where the sample is
#     deep enough (real-data depth series for concordance and stability).
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/phase3/full/p3_env.sh
IFS=$'\t' read -r S FILES LAYOUT FRAC < <(task_line "$S2_META/tasks_align.tsv" "$SLURM_ARRAY_TASK_ID")
OD="$P3_BAM/$S"; BAM="$OD/$S.bam"
[[ -f "$OD/.align_done" ]] || { echo "not aligned: $S"; exit 3; }
[[ -f "$OD/.derive_done" ]] && { echo "SKIP $S"; exit 0; }
aln_modules
DEP=$(awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {L+=$2; B+=$3} END{printf "%.4f", B/L}' "$OD/$S.mosdepth.summary.txt")
echo "autosomal depth $DEP"
{
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${S}_full" "$S" "$BAM" "$DEP" full "$LAYOUT"
  for T in 1 0.1; do
    if awk -v d="$DEP" -v t="$T" 'BEGIN{exit !(d > 1.3*t)}'; then
      L=$( [[ $T == 1 ]] && echo d1 || echo d0p1 ); X="$OD/${S}_${L}.bam"
      FR=$(awk -v d="$DEP" -v t="$T" 'BEGIN{printf "%.6f", t/d}')
      [[ -s "$X.bai" ]] || { samtools view -@ 4 -b -s "11${FR#0}" -o "$X" "$BAM" && samtools index "$X"; }
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${S}_${L}" "$S" "$X" "$T" "$L" "$LAYOUT"
    fi
  done
} > "$OD/derived.tsv"
touch "$OD/.derive_done"; cat "$OD/derived.tsv"
