#!/bin/bash
#SBATCH --job-name=cfCNV_depth
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/depth_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/depth_%j.err

set -euo pipefail

source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/config.sh

# Load ONLY samtools for this job
module purge
module load SAMtools/1.16.1-GCC-11.3.0

BAM=$1
SAMPLE=$2

mkdir -p $QC/depth

echo "=========================================="
echo "Sample: $SAMPLE"
echo "BAM: $BAM"
echo "Started: $(date)"
echo "=========================================="

# Confirm BAM exists
if [ ! -f "$BAM" ]; then
    echo "ERROR: BAM does not exist:"
    echo "$BAM"
    exit 1
fi

# Basic integrity check
samtools quickcheck -v "$BAM"

# Save per-chromosome coverage statistics
samtools coverage "$BAM" \
    > $QC/depth/${SAMPLE}.coverage.tsv

# Calculate weighted mean autosomal sequencing depth.
# samtools coverage columns:
# 1=rname
# 2=startpos
# 3=endpos
# 4=numreads
# 5=covbases
# 6=coverage
# 7=meandepth
# 8=meanbaseq
# 9=meanmapq

awk '
BEGIN {
    FS="\t"
}

NR > 1 {

    chr=$1
    gsub(/^chr/, "", chr)

    if (chr ~ /^([1-9]|1[0-9]|2[0-2])$/) {

        len=$3-$2+1

        weighted_depth += $7 * len
        total_length += len
    }
}

END {

    if (total_length > 0)
        printf "%.6f\n", weighted_depth/total_length
    else
        print "NA"
}
' $QC/depth/${SAMPLE}.coverage.tsv \
> $QC/depth/${SAMPLE}.mean_autosomal_depth.txt

echo
echo "Estimated autosomal mean depth:"
cat $QC/depth/${SAMPLE}.mean_autosomal_depth.txt

echo
echo "Finished: $(date)"
