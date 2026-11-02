#!/bin/bash
set -euo pipefail

source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/config.sh

SIM1=$SIM/healthy_cfDNA_spikein
MAN=$SIM1/manifests

if ! command -v sbatch >/dev/null 2>&1
then
    echo "ERROR: sbatch is not available on this host."
    echo "Run this script from your Slurm submission/login node."
    exit 1
fi

NPARENT=$(($(wc -l < "$MAN/parent_tasks.tsv") - 1))
NSPIKE=$(($(wc -l < "$MAN/spike_tasks.tsv") - 1))
NDEPTH=$(($(wc -l < "$MAN/depth_tasks.tsv") - 1))
NNEG=$(($(wc -l < "$MAN/negative_tasks.tsv") - 1))

echo "======================================================"
echo "FULL cfCNV SIMULATION 1"
echo "======================================================"
echo "Parent jobs:          $NPARENT"
echo "Spike jobs:           $NSPIKE"
echo "Positive depth jobs:  $NDEPTH"
echo "Negative jobs:        $NNEG"
echo "======================================================"

# ------------------------------------------------------------
# 3 parent BAMs
# ------------------------------------------------------------

PARENT_JOB=$(sbatch \
    --parsable \
    --array=0-$((NPARENT-1))%3 \
    "$PROJECT/00_metadata/07_make_parent_array.sh")

echo "Parent job ID: $PARENT_JOB"

# ------------------------------------------------------------
# 15 spike-in parent BAMs
# Starts only after ALL parent jobs succeed.
# ------------------------------------------------------------

SPIKE_JOB=$(sbatch \
    --parsable \
    --dependency=afterok:$PARENT_JOB \
    --array=0-$((NSPIKE-1))%3 \
    "$PROJECT/00_metadata/08_spike_array.sh")

echo "Spike job ID: $SPIKE_JOB"

# ------------------------------------------------------------
# 19 negative-control depth BAMs
# Only need parent BAMs.
# ------------------------------------------------------------

NEG_JOB=$(sbatch \
    --parsable \
    --dependency=afterok:$PARENT_JOB \
    --array=0-$((NNEG-1))%8 \
    "$PROJECT/00_metadata/09_downsample_array.sh" \
    "$MAN/negative_tasks.tsv")

echo "Negative job ID: $NEG_JOB"

# ------------------------------------------------------------
# 95 positive simulated depth BAMs
# Requires spike-in BAMs first.
# ------------------------------------------------------------

DEPTH_JOB=$(sbatch \
    --parsable \
    --dependency=afterok:$SPIKE_JOB \
    --array=0-$((NDEPTH-1))%10 \
    "$PROJECT/00_metadata/09_downsample_array.sh" \
    "$MAN/depth_tasks.tsv")

echo "Positive-depth job ID: $DEPTH_JOB"

cat > "$MAN/submitted_jobs.txt" <<JOBS
parent_job=$PARENT_JOB
spike_job=$SPIKE_JOB
negative_job=$NEG_JOB
positive_depth_job=$DEPTH_JOB
JOBS

echo
echo "======================================================"
echo "ALL SIMULATION-1 JOBS SUBMITTED"
echo "======================================================"
cat "$MAN/submitted_jobs.txt"
