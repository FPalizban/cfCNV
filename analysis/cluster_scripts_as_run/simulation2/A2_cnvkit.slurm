#!/bin/bash
#SBATCH --job-name=S2_cnvkit
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/cnvkit_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/cnvkit_%A_%a.err
# A2  CNVkit WGS mode, flat GRCh38 reference per bin size (built earlier), 624 tasks.
#     A pooled-normal reference is not possible here: HG008 has a single donor and its
#     normal pool is the source of every mixture, so it would normalise the sample
#     against its own reads. The flat reference is reported as such.
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/simulation2/full/s2_env.sh
module load CNVkit/0.9.10-foss-2022b-R-4.2.2
IFS=$'\t' read -r S BAM D B < <(task_line "$S2_META/tasks_cnvkit.tsv" "$SLURM_ARRAY_TASK_ID")
CNN="$S2_OUT/CNVkit/reference/flat_GRCh38_${B}.cnn"; test -s "$CNN"
OD="$S2_OUT/CNVkit/$S/$B"; mkdir -p "$OD"
[[ -f "$OD/.done" ]] && { echo "SKIP $S $B"; exit 0; }
/usr/bin/time -v cnvkit.py batch "$BAM" -m wgs -r "$CNN" -d "$OD" -p 4
ls "$OD"/*.call.cns >/dev/null
touch "$OD/.done"; echo "CNVKIT_DONE $S $B"
