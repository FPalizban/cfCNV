#!/bin/bash
#SBATCH --job-name=cfCNV_actualdepth
#SBATCH --cpus-per-task=20
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/actualdepth_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/actualdepth_%A_%a.err

set -euo pipefail

source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/config.sh
module load SAMtools/1.16.1-GCC-11.3.0

MAN=$SIM/healthy_cfDNA_spikein/manifests/analysis_manifest.tsv

LINE=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+2)) 'NR==n' "$MAN")

IFS=$'\t' read -r ID BACKGROUND CLASS REP TARGET BAM TRUTH <<< "$LINE"

OUTDIR=$QC/simulation1_depths
mkdir -p "$OUTDIR"

COV=$OUTDIR/${ID}.coverage.tsv
DEPTH=$OUTDIR/${ID}.mean_depth.txt

echo "ID=$ID"
echo "Target=${TARGET}x"
echo "BAM=$BAM"

samtools coverage "$BAM" > "$COV"

awk '
BEGIN {FS="\t"}

NR>1 {
    chr=$1
    gsub(/^chr/,"",chr)

    if(chr ~ /^([1-9]|1[0-9]|2[0-2])$/) {
        len=$3-$2+1
        weighted += $7*len
        total += len
    }
}

END {
    if(total>0)
        printf "%.6f\n",weighted/total
    else
        print "NA"
}
' "$COV" > "$DEPTH"

echo "Measured=$(cat "$DEPTH")x"
