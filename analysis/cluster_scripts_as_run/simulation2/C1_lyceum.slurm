#!/bin/bash
#SBATCH --job-name=S2_lyceum
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/lyceum_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/lyceum_%A_%a.err
# C1  LYCEUM in its NATIVE setting (GRCh38 BAMs, its own hg38 exon targets, released pipeline),
#     keeping the per-exon probabilities. One task per mixture (104 tasks).
#     Simulation 1 needed a GRCh37 adaptation; this run has none.
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/simulation2/full/s2_env.sh
module load SAMtools/1.16.1-GCC-11.3.0
unset PYTHONPATH; export OMP_NUM_THREADS=8
IFS=$'\t' read -r S BAM D X < <(task_line "$S2_META/tasks_lyceum.tsv" "$SLURM_ARRAY_TASK_ID")
W="$S2_OUT/LYCEUM/$S"; mkdir -p "$W"/{read_depths,processed,stats,calls}
[[ -f "$W/.done" ]] && { echo "SKIP $S"; exit 0; }
T="$LYC_REPO/hg38_hglft_genome_64dc_dcbaa0_unique.bed"; cd "$W"
[[ -s "read_depths/$S.bam.txt" ]] || sambamba depth base -t 8 -L "$T" "$BAM" > "read_depths/$S.bam.txt"
[[ -s "processed/${S}_labeled_data.npy" ]] || "$PYL" "$LYC_REPO/scripts/preprocess_sample.py" --readdepth read_depths --output processed --target "$T"
[[ -s stats/s2_stats_lookup.npy ]] || "$PYL" "$LYC_REPO/scripts/mean_std_calculator.py" --processed_samples processed --setname s2 --output stats
"$PYL" "$LYC_REPO/scripts/LYCEUM_call_keepconf.py" --model lyceum --input processed --output calls \
   --cnv exonlevel --batch_size 16 --stats_lookup stats/s2_stats_lookup.npy
test -s "calls/${S}_exonlevel.csv"; rm -f "read_depths/$S.bam.txt"
touch "$W/.done"; echo "LYCEUM_DONE $S"
