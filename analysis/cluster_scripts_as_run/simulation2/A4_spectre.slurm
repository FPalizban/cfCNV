#!/bin/bash
#SBATCH --job-name=S2_spectre
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/spectre_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/spectre_%A_%a.err
# A4  Spectre exactly as in the Simulation 1 re-run: mosdepth coverage at each benchmark bin
#     size, ENCODE blacklist, --cancer (HG008 is a tumour; germline and cancer mode gave the
#     same result in Simulation 1). 624 tasks.
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/simulation2/full/s2_env.sh
IFS=$'\t' read -r S BAM D B < <(task_line "$S2_META/tasks_spectre.tsv" "$SLURM_ARRAY_TASK_ID")
OD="$S2_OUT/Spectre/$S/$B"; mkdir -p "$OD/mosdepth" "$OD/calls"
[[ -f "$OD/.done" ]] && { echo "SKIP $S $B"; exit 0; }
"$MOSDEPTH" -t 3 -n -x -b "$B" -Q 20 "$OD/mosdepth/$S" "$BAM"
EXTRA=(); grep -q -- "--blacklist" "$S2/spectre_help.txt" && EXTRA+=(--blacklist "$BLACKLIST38")
MINLEN=$(( B > 50000 ? B : 50000 ))
/usr/bin/time -v "$SPECTRE" CNVCaller --coverage "$OD/mosdepth/$S.regions.bed.gz" --sample-id "$S" \
   --output-dir "$OD/calls" --reference "$REF38" --min-cnv-len "$MINLEN" --threads 3 --cancer \
   --only-chr "$(printf 'chr%s,' {1..22} | sed 's/,$//')" ${EXTRA[@]+"${EXTRA[@]}"}
ls "$OD"/calls/*_cnv.bed.gz >/dev/null
touch "$OD/.done"; echo "SPECTRE_DONE $S $B"
