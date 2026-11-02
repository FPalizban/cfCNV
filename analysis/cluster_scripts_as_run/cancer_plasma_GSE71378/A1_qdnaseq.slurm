#!/bin/bash
#SBATCH --job-name=G7_qdnaseq
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/gse71378/qdnaseq_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/gse71378/qdnaseq_%A_%a.err
# A1  QDNAseq, same workflow as Simulation 1, QDNAseq.hg38 bins, 6 bin sizes (624 tasks).
#     Also saves the segmented object (.rds) that ACE (A5) re-uses.
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/gse71378/full/p3_env.sh
IFS=$'\t' read -r S BAM D B < <(task_line "$S2_META/tasks_qdnaseq.tsv" "$SLURM_ARRAY_TASK_ID")
OD="$S2_OUT/QDNAseq/$S/$B"; mkdir -p "$OD"
[[ -f "$OD/.done" && -s "$OD/$S.qdnaseq_seg.rds" ]] && { echo "SKIP $S $B"; exit 0; }
r_modules
/usr/bin/time -v Rscript - "$S" "$BAM" "$B" "$OD" <<'RS'
a <- commandArgs(TRUE); s <- a[1]; bam <- a[2]; b <- as.integer(a[3]); od <- a[4]
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages({library(QDNAseq); library(QDNAseq.hg19)})
x <- getBinAnnotations(binSize=b/1000, genome="hg19")
x <- binReadCounts(x, bamfiles=bam, bamnames=s)
n <- sum(Biobase::assayDataElement(x, "counts")); cat("reads counted in bins:", n, "\n")
if (n < 1000) stop("almost no reads counted - chromosome naming mismatch")
x <- applyFilters(x, residual=TRUE, blacklist=TRUE)
x <- estimateCorrection(x); x <- correctBins(x); x <- normalizeBins(x)
x <- smoothOutlierBins(x); x <- segmentBins(x, transformFun="sqrt")
x <- normalizeSegmentedBins(x); x <- callBins(x)
saveRDS(x, file.path(od, paste0(s, ".qdnaseq_seg.rds")))
exportBins(x, file=file.path(od, paste0(s, ".calls.vcf")), format="vcf", type="calls")
exportBins(x, file=file.path(od, paste0(s, ".segments.igv")), format="igv", type="segments", logTransform=TRUE)
RS
touch "$OD/.done"; echo "QDNASEQ_DONE $S $B"
