#!/usr/bin/env Rscript
# QDNAseq, same steps and settings as the cfCNV paper.
# usage: run_qdnaseq.R <sample_id> <bam> <bin_kb> <hg19|hg38> <prefix>
# writes <prefix>.calls.vcf, <prefix>.segments.igv, <prefix>.qdnaseq_seg.rds (re-used by ACE)
a <- commandArgs(TRUE); s <- a[1]; bam <- a[2]; b <- as.integer(a[3]); g <- a[4]; pre <- a[5]
suppressPackageStartupMessages(library(QDNAseq))
if (g == "hg19") suppressPackageStartupMessages(library(QDNAseq.hg19)) else suppressPackageStartupMessages(library(QDNAseq.hg38))
x <- getBinAnnotations(binSize = b, genome = g)
x <- binReadCounts(x, bamfiles = bam, bamnames = s)
n <- sum(Biobase::assayDataElement(x, "counts")); cat("reads counted in bins:", n, "\n")
if (n < 1000) stop("almost no reads counted - check chromosome naming of the BAM against the bin annotation")
x <- applyFilters(x, residual = TRUE, blacklist = TRUE)
x <- estimateCorrection(x); x <- correctBins(x); x <- normalizeBins(x)
x <- smoothOutlierBins(x); x <- segmentBins(x, transformFun = "sqrt")
x <- normalizeSegmentedBins(x); x <- callBins(x)
saveRDS(x, paste0(pre, ".qdnaseq_seg.rds"))
exportBins(x, file = paste0(pre, ".calls.vcf"), format = "vcf", type = "calls")
exportBins(x, file = paste0(pre, ".segments.igv"), format = "igv", type = "segments", logTransform = TRUE)
