#!/bin/bash
#SBATCH --job-name=P3_bwa_index
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/phase3/bwa_index_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/phase3/bwa_index_%j.err
# X0  BWA index of the GRCh38 GIABv3 reference (only if it does not exist yet; ~2-3 h).
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/phase3/full/p3_env.sh
aln_modules
[[ -s "$REF38.bwt" ]] && { echo "index exists"; exit 0; }
bwa index -a bwtsw "$REF38"
echo BWA_INDEX_DONE
