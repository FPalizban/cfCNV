#!/bin/bash
# 1_submit_bams.sh - login node:  bash 1_submit_bams.sh
# Stage 1 without realignment: G1 (mirror BAM -> 12x / 1x / 0.1x analysis BAMs) -> 13 task lists.
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/gse71378/full/p3_env.sh
cd "$S2"; rm -f "$S2_META/STAGE1_COMPLETE"
N=$(grep -vcE '^#|^\s*$' samples.txt)
g=$(sbatch --parsable --array=0-$((N-1)) G1_bam.slurm)
t=$(sbatch --parsable --dependency=afterany:$g 13_tasks.slurm)
echo "bam=$g tasks=$t" | tee "$S2_LOG/stage1_bam_jobids.txt"
echo "Submitted $N samples. Check: bash 3_check_p3.sh"
