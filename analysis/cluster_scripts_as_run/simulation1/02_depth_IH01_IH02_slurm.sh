#!/bin/bash
#SBATCH --job-name=cfCNV_depth
#SBATCH --array=0-1
#SBATCH --cpus-per-task=25
#SBATCH --mem=200G
#SBATCH --time=3-00:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/depth_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/depth_%A_%a.err

set -euo pipefail

# ------------------------------------------------------------
# Project configuration
# ------------------------------------------------------------
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/config.sh

# Load only the module needed for this task
module load SAMtools/1.16.1-GCC-11.3.0

mkdir -p "$QC/depth"
mkdir -p "$LOGS"

# ------------------------------------------------------------
# Samples
# Array task 0 = IH01
# Array task 1 = IH02
# ------------------------------------------------------------
SAMPLES=(IH01 IH02)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
BAM="$SHENDURE_DIR/${SAMPLE}.bam"

COVERAGE="$QC/depth/${SAMPLE}.coverage.tsv"
DEPTH_OUT="$QC/depth/${SAMPLE}.mean_autosomal_depth.txt"

echo "======================================================"
echo "cfCNV sequencing-depth calculation"
echo "======================================================"
echo "Date:              $(date)"
echo "Host:              $(hostname)"
echo "User:              $USER"
echo "Slurm job ID:      ${SLURM_JOB_ID:-NA}"
echo "Array task ID:     ${SLURM_ARRAY_TASK_ID:-NA}"
echo "Sample:            $SAMPLE"
echo "BAM:               $BAM"
echo "Coverage output:   $COVERAGE"
echo "Depth output:      $DEPTH_OUT"
echo "======================================================"

# ------------------------------------------------------------
# Confirm BAM exists
# ------------------------------------------------------------
if [[ ! -s "$BAM" ]]; then
    echo "ERROR: BAM missing or empty:"
    echo "$BAM"
    exit 1
fi

echo
echo "BAM size:"
ls -lh "$BAM"

# ------------------------------------------------------------
# BAM integrity check
# ------------------------------------------------------------
echo
echo "Running samtools quickcheck..."

if samtools quickcheck -v "$BAM"; then
    echo "BAM quickcheck: PASS"
else
    echo "BAM quickcheck: FAIL"
    exit 1
fi

# ------------------------------------------------------------
# Calculate chromosome-level coverage
#
# IMPORTANT:
# samtools coverage in this installed version does not use -@
# ------------------------------------------------------------
echo
echo "Starting samtools coverage at $(date)..."

/usr/bin/time -v \
samtools coverage "$BAM" \
> "$COVERAGE" \
2> "$QC/depth/${SAMPLE}.coverage.time.txt"

echo "samtools coverage finished at $(date)"

# ------------------------------------------------------------
# Verify coverage output
# ------------------------------------------------------------
if [[ ! -s "$COVERAGE" ]]; then
    echo "ERROR: coverage file is empty."
    exit 1
fi

echo
echo "First rows of coverage:"
head -5 "$COVERAGE"

# ------------------------------------------------------------
# Calculate weighted mean autosomal depth
#
# samtools coverage columns:
# 1 rname
# 2 startpos
# 3 endpos
# 4 numreads
# 5 covbases
# 6 coverage
# 7 meandepth
# 8 meanbaseq
# 9 meanmapq
# ------------------------------------------------------------
awk '
BEGIN {
    FS="\t"
}

NR > 1 {

    chr=$1
    gsub(/^chr/, "", chr)

    if (chr ~ /^([1-9]|1[0-9]|2[0-2])$/) {

        len=$3-$2+1

        weighted += $7 * len
        total += len
    }
}

END {

    if (total > 0)
        printf "%.6f\n", weighted/total
    else
        print "NA"
}
' "$COVERAGE" > "$DEPTH_OUT"

# ------------------------------------------------------------
# Validate final value
# ------------------------------------------------------------
VALUE=$(cat "$DEPTH_OUT")

echo
echo "======================================================"
echo "RESULT"
echo "======================================================"
echo "$SAMPLE mean autosomal depth = ${VALUE}x"

if echo "$VALUE" | grep -Eq '^[0-9]+([.][0-9]+)?$'; then
    echo "Depth calculation: PASS"
else
    echo "Depth calculation: FAIL"
    exit 1
fi

echo
echo "Finished at: $(date)"
echo "======================================================"
