#!/bin/bash
#SBATCH --job-name=bicseq2_preflight
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/bicseq2_preflight_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/bicseq2_preflight_%j.err
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/00_env.sh
"$TOOLS/bin/run_bicseq2" bash -lc 'find / -type f \( -name "NBICseq-norm.pl" -o -name "NBICseq-seg.pl" -o -name "prep_mappability.sh" \) 2>/dev/null | sort'
echo "BIC-seq2 installation OK."
echo "Before production BIC-seq2 runs, provide GRCh37 per-chromosome FASTA + read-length-specific mappability files."
