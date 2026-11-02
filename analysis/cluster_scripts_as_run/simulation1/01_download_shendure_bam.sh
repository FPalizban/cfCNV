#!/bin/bash
#SBATCH --job-name=ShendureBAM
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=3-00:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/ShendureBAM.%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/ShendureBAM.%A_%a.err

set -euo pipefail

source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/config.sh

SAMPLE=$1

BASE=https://kircherlab.bihealth.org/download/cfDNA

mkdir -p $SHENDURE_DIR
cd $SHENDURE_DIR

echo "Downloading ${SAMPLE}.bam"

wget --continue \
     --tries=20 \
     --timeout=60 \
     ${BASE}/${SAMPLE}.bam

wget --continue \
     --tries=20 \
     --timeout=60 \
     ${BASE}/${SAMPLE}.bam.bai

echo "Finished $SAMPLE"
