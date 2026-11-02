#!/usr/bin/env Rscript
# ACE cellularity (tumour fraction) from a QDNAseq segmented object, as in the cfCNV paper:
# ploidy fixed at 2 (singlemodel) and free ploidy (squaremodel).
# usage: run_ace.R <sample_id> <qdnaseq_seg.rds> <out.tsv>
a <- commandArgs(TRUE); s <- a[1]; rds <- a[2]; out <- a[3]
suppressPackageStartupMessages({library(QDNAseq); library(ACE)})
x <- readRDS(rds)
m <- singlemodel(x, QDNAseqobjectsample = 1, ploidy = 2)
best2 <- m$minima[which.min(m$rerror)]
sq <- tryCatch(squaremodel(x, QDNAseqobjectsample = 1), error = function(e) NULL)
bestc <- NA; bestp <- NA
if (!is.null(sq) && !is.null(sq$minimadf)) {
  d <- sq$minimadf; i <- which.min(d$error); bestc <- d$cellularity[i]; bestp <- d$ploidy[i] }
write.table(data.frame(sample_id = s, ace_cellularity_ploidy2 = best2, n_minima = length(m$minima),
                       ace_cellularity_free = bestc, ace_ploidy_free = bestp),
            out, sep = "\t", quote = FALSE, row.names = FALSE)
