#!/bin/bash
# =====================================================================
# cfCNV - one-time install of ACE (+ QDNAseq.hg38 bins) and Fragle.
# Run on a login node (needs internet).
#   R_MODULE  = the R module you used for QDNAseq (edit below)
# =====================================================================
set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project
TOOLS="$PROJECT/06_tools"
R_MODULE="${R_MODULE:-R}"          # e.g. R/4.3.2  <- set to your QDNAseq module

# ---------------- ACE + QDNAseq.hg38 ----------------
module load "$R_MODULE"
Rscript -e '
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
  for (p in c("QDNAseq","QDNAseq.hg19","ACE","Biobase")) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, update=FALSE, ask=FALSE)
  if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", repos="https://cloud.r-project.org")
  if (!requireNamespace("QDNAseq.hg38", quietly=TRUE)) remotes::install_github("asntech/QDNAseq.hg38@main", upgrade="never")
  for (p in c("QDNAseq","QDNAseq.hg19","QDNAseq.hg38","ACE")) cat(p, as.character(packageVersion(p)), "\n")
' | tee "$TOOLS/ACE_versions.txt"

# ---------------- Fragle (authors' conda-pack env) ----------------
mkdir -p "$TOOLS/ML_based"
cd "$TOOLS/ML_based"
if [[ ! -d FRAGLE ]]; then
    git clone https://github.com/skandlab/FRAGLE.git
fi
cd FRAGLE
git rev-parse HEAD > "$TOOLS/Fragle_commit.txt"
if [[ ! -d env ]]; then
    wget -O Fragle.tar.gz "https://zenodo.org/records/13968359/files/Fragle.tar.gz?download=1"
    mkdir -p env && tar -xzf Fragle.tar.gz -C env && rm Fragle.tar.gz
fi
source env/bin/activate
python -c "import torch, sklearn, pysam; print('Fragle env OK', torch.__version__, sklearn.__version__)"
python main.py -h | head -20
source env/bin/deactivate
echo "INSTALL DONE"
