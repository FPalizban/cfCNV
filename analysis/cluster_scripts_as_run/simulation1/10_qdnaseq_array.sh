#!/bin/bash
#SBATCH --job-name=qdnaseq_s1
#SBATCH --cpus-per-task=25
#SBATCH --mem=150G
#SBATCH --time=08:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/qdnaseq_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/qdnaseq_%A_%a.err
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/00_env.sh
line=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+2)) 'NR==n{print}' "$SIM1_MAN/qdnaseq_tasks.tsv")
IFS=$'\t' read -r sample bam bg depth rep cls bin <<< "$line"
module load R/4.2.2-foss-2022b
module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
export R_LIBS_USER="$TOOLS/Rlib"
od="$SIM1_OUT/QDNAseq/$sample/$bin"; mkdir -p "$od"
Rscript - "$sample" "$bam" "$bin" "$od" <<'RS'
a=commandArgs(TRUE); s=a[1]; bam=a[2]; b=as.integer(a[3]); od=a[4]
.libPaths(c(Sys.getenv("R_LIBS_USER"),.libPaths()))
suppressPackageStartupMessages({library(QDNAseq);library(QDNAseq.hg19)})
x=getBinAnnotations(binSize=as.integer(b/1000),genome="hg19")
x=binReadCounts(x,bamfiles=bam)
x=applyFilters(x,residual=TRUE,blacklist=TRUE)
x=estimateCorrection(x); x=correctBins(x); x=normalizeBins(x)
x=smoothOutlierBins(x); x=segmentBins(x,transformFun="sqrt")
x=normalizeSegmentedBins(x); x=callBins(x)
# Export QDNAseq CNV calls using the official QDNAseq interface.
# SEG contains non-neutral called segments.
QDNAseq::exportBins(
    x,
    file=file.path(od, paste0(s, ".calls.tsv")),
    format="seg",
    type="calls"
)

# VCF additionally preserves explicit DEL/DUP call direction,
# QDNAseq call score, number of bins, and segment log2 ratio.
QDNAseq::exportBins(
    x,
    file=file.path(od, paste0(s, ".calls.vcf")),
    format="vcf",
    type="calls"
)
RS
