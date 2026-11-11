#!/bin/bash
# =====================================================================
# 1_run_all.sh  - run on the LOGIN node:   bash 1_run_all.sh
#   1. makes $MLO/manifest.tsv from the four existing manifests (seconds)
#   2. submits A (GC) + B (counts, array) -> C (dataset) -> D (baseline); each step waits for the one before
# Re-running is safe: finished BAMs are skipped, an older dataset folder is renamed, nothing is deleted.
# =====================================================================
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/ml_dataset/ml_env.sh
cd "$ML"
pyd m0_make_manifest.py
N=$(( $(wc -l < "$MLO/manifest.tsv") - 1 ))
NT=$(( N < 500 ? N : 500 ))                 # array tasks (each handles rows i, i+NT, ...)
a=$(sbatch --parsable A_gc_bins.slurm)
b=$(sbatch --parsable --export=ALL,NTASKS_B=$NT --array=0-$(( NT - 1 ))%30 B_count_bins.slurm)
c=$(sbatch --parsable --dependency=afterok:$a,afterany:$b C_assemble.slurm)
d=$(sbatch --parsable --dependency=afterok:$c D_baseline.slurm)
echo "gc=$a counts=$b assemble=$c baseline=$d bams=$N" | tee "$MLLOG/ml_jobids.txt"
echo "Submitted. Check progress any time with:   bash $ML/2_check.sh"
