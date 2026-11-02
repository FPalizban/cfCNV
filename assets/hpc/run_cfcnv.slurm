#!/bin/bash
#SBATCH --job-name=cfCNV_head
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --output=cfCNV_head_%j.log
# The "head" job: Nextflow itself runs here and submits every tool step as its own Slurm job.
# Edit the paths below, then:  sbatch run_cfcnv.slurm
set -euo pipefail
module load singularity/4.4.1 2>/dev/null || true
module load Java/17 2>/dev/null || true                      # any Java 17 or newer (see: module spider Java)
export NXF_SINGULARITY_CACHEDIR=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/singularity_cache
export NXF_OPTS='-Xms1g -Xmx4g'
NF=${HOME}/bin/nextflow
$NF run FPalizban/cfCNV -r v1.0.0 -profile chop \
    --input samples.csv \
    --fasta /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/02_reference/GRCh37/fasta/human_g1k_v37.fasta \
    --genome GRCh37 \
    --bicseq2_resources /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/06_tools/BICseq2_resources/GRCh37_101mer \
    --bicseq2_map_prefix human_g1k_v37.101mer \
    --outdir results -resume
