#!/bin/bash
#SBATCH --job-name=P3_ace
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/phase3/ace_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/phase3/ace_%A_%a.err
# A5  ACE tumour-fraction (cellularity) estimate from the QDNAseq segmentation of A1
#     (100 kb, 500 kb, 1 Mb; 312 tasks). Fits ploidy 2 (singlemodel) and free ploidy (squaremodel).
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/phase3/full/p3_env.sh
IFS=$'\t' read -r S BAM D B < <(task_line "$S2_META/tasks_ace.tsv" "$SLURM_ARRAY_TASK_ID")
RDS="$S2_OUT/QDNAseq/$S/$B/$S.qdnaseq_seg.rds"; test -s "$RDS"
OD="$S2_OUT/ACE/$S/$B"; mkdir -p "$OD"
[[ -f "$OD/.done" ]] && { echo "SKIP $S $B"; exit 0; }
r_modules
Rscript - "$S" "$RDS" "$OD" <<'RS'
a <- commandArgs(TRUE); s <- a[1]; rds <- a[2]; od <- a[3]
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
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
            file.path(od, paste0(s, ".ace.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
cat("ACE", s, "cellularity(ploidy 2) =", best2, " free:", bestc, bestp, "\n")
RS
touch "$OD/.done"
