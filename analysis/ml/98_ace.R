#!/usr/bin/env Rscript
# =====================================================================
# cfCNV - QDNAseq segmentation + ACE tumour-cellularity estimate, 1 BAM.
#
# Rscript 98_ace.R <bam> <sample_id> <build hg19|hg38> <bin_kb> <outdir> [ploidy=2]
#
# Writes <outdir>/<sample>.ace.tsv           (best and last fit, n minima)
#        <outdir>/<sample>.ace_segments.tsv  (absolute-CN segments at best fit)
#        <outdir>/<sample>.qdnaseq_seg.rds
# =====================================================================
suppressPackageStartupMessages({
  library(QDNAseq); library(ACE); library(Biobase)
})

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 5) stop("usage: 98_ace.R bam sample build bin_kb outdir [ploidy]")
bam <- a[1]; sample <- a[2]; build <- a[3]; bin_kb <- as.integer(a[4])
outdir <- a[5]; ploidy <- if (length(a) >= 6) as.integer(a[6]) else 2L
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (build == "hg38") {
  suppressPackageStartupMessages(library(QDNAseq.hg38))
} else if (build == "hg19") {
  suppressPackageStartupMessages(library(QDNAseq.hg19))
} else stop("build must be hg19 or hg38")

bins <- getBinAnnotations(binSize = bin_kb, genome = build)
rc   <- binReadCounts(bins, bamfiles = bam, bamnames = sample)

tot <- sum(assayDataElement(rc, "counts"))
cat("reads counted in bins:", tot, "\n")
if (tot == 0) stop("0 reads counted: chromosome names in the BAM (chr1 vs 1) ",
                   "do not match the QDNAseq bins for ", build)

rc  <- applyFilters(rc, residual = TRUE, blacklist = TRUE)
rc  <- estimateCorrection(rc)
cn  <- correctBins(rc)
cn  <- normalizeBins(cn)
cn  <- smoothOutlierBins(cn)
seg <- segmentBins(cn, transformFun = "sqrt")
seg <- normalizeSegmentedBins(seg)
saveRDS(seg, file.path(outdir, paste0(sample, ".qdnaseq_seg.rds")))

model <- singlemodel(seg, QDNAseqobjectsample = 1, ploidy = ploidy)
if (!all(c("minima", "rerror") %in% names(model)))
  stop("ACE singlemodel() fields changed; found: ", paste(names(model), collapse = ","))

best <- model$minima[which.min(model$rerror)]
last <- tail(model$minima, 1)

write.table(data.frame(sample_id = sample, build = build, bin_kb = bin_kb,
                       ploidy = ploidy, ace_bestfit = best, ace_lastfit = last,
                       n_minima = length(model$minima), reads_in_bins = tot),
            file.path(outdir, paste0(sample, ".ace.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

segs <- tryCatch(
  getadjustedsegments(seg, QDNAseqobjectsample = 1, cellularity = best, ploidy = ploidy),
  error = function(e) { message("segment export skipped: ", conditionMessage(e)); NULL })
if (!is.null(segs))
  write.table(segs, file.path(outdir, paste0(sample, ".ace_segments.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

file.create(file.path(outdir, ".done"))
cat("ACE DONE", sample, "best cellularity =", best, "\n")
