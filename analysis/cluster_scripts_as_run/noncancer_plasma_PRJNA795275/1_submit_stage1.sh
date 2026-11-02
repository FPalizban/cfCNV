#!/bin/bash
# =====================================================================
# 1_submit_stage1.sh - Phase 3 stage 1 (login node):   bash 1_submit_stage1.sh
#   [X0 BWA index, only if missing] -> 10 download (4 at a time) -> 11 align (6 at a time)
#   -> 12 derive (depth + 1x / 0.1x copies) -> 13 task lists.
# Job numbers -> logs/phase3/stage1_jobids.txt.  Progress: bash 3_check_p3.sh
# =====================================================================
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/phase3/full/p3_env.sh
cd "$S2"
last() { echo $(( $(wc -l < "$S2_META/tasks_$1.tsv") - 2 )); }
idx=""; [[ -s "$REF38.bwt" ]] || idx=$(sbatch --parsable X0_bwa_index.slurm)
d=$(sbatch --parsable --array=0-$(last download)%4 10_download.slurm)
dep="afterany:$d"; [[ -n "$idx" ]] && dep="$dep,afterok:$idx"
a=$(sbatch --parsable --dependency=$dep --array=0-$(last align)%6 11_align.slurm)
v=$(sbatch --parsable --dependency=afterany:$a --array=0-$(last align) 12_derive.slurm)
t=$(sbatch --parsable --dependency=afterany:$v 13_tasks.slurm)
echo "index=${idx:-none} download=$d align=$a derive=$v tasks=$t" | tee "$S2_LOG/stage1_jobids.txt"
echo "Submitted stage 1. Check with: bash 3_check_p3.sh"
