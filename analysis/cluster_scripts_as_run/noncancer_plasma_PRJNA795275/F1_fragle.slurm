#!/bin/bash
#SBATCH --job-name=P3_fragle
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/phase3/fragle_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/phase3/fragle_%A_%a.err
# F1  Fragle ctDNA fraction from the fragment-length distribution (paired-end BAMs only),
#     raw WGS mode, GRCh38.
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/phase3/full/p3_env.sh
IFS=$'\t' read -r S BAM D X < <(task_line "$S2_META/tasks_fragle.tsv" "$SLURM_ARRAY_TASK_ID")
OD="$S2_OUT/Fragle/$S"; mkdir -p "$OD/in"
[[ -f "$OD/.done" ]] && { echo "SKIP $S"; exit 0; }
aln_modules
ln -sf "$BAM" "$OD/in/$S.bam"; ln -sf "$BAM.bai" "$OD/in/$S.bam.bai"
set +u; source "$FRAGLE_ENV/bin/activate"; set -u
cd "$FRAGLE_SRC"
/usr/bin/time -v python main.py --input "$OD/in" --output "$OD/out" --mode R --genome_build hg38 --cpu 8
ls "$OD/out/Fragle.csv" && touch "$OD/.done"; echo "FRAGLE_DONE $S"
