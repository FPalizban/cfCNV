#!/bin/bash
# Run once on the cluster (login node is fine) to record the tool versions behind the paper's results:
#   bash record_tool_versions.sh > tool_versions_used.txt
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/config.sh
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/00_env.sh
echo "== R / Bioconductor modules"; module load R/4.2.2-foss-2022b R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2 2>/dev/null
Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); for (p in c("QDNAseq","QDNAseq.hg19","QDNAseq.hg38","ACE","ichorCNA","HMMcopy","DNAcopy")) cat(p, tryCatch(as.character(packageVersion(p)), error=function(e) "not found"), "\n")'
echo "== ichorCNA source"; git -C "$ICHOR_DIR" describe --tags --always 2>/dev/null; git -C "$ICHOR_DIR" log -1 --format='%H %cd' 2>/dev/null
echo "== CNVkit"; module load CNVkit/0.9.10-foss-2022b-R-4.2.2 2>/dev/null && cnvkit.py version
echo "== Spectre"; "$SPECTRE" --version 2>&1 | head -3; head -20 "$SPECTRE"
echo "== mosdepth"; "$MOSDEPTH" --version
echo "== readCounter (HMMcopy utils)"; ls -l "$READCOUNTER"
echo "== samtools / bedtools"; module load SAMtools/1.16.1-GCC-11.3.0 BEDTools/2.30.0-GCC-11.3.0 2>/dev/null; samtools --version | head -1; bedtools --version
echo "== BIC-seq2 image"; ls -l "$TOOLS"/containers/*.sif 2>/dev/null; cat "$TOOLS/bin/run_bicseq2" 2>/dev/null | head -20
echo "== Fragle"; ls "$TOOLS"/FRAGLE-main 2>/dev/null | head; "$TOOLS"/Fragle_env/bin/python -c "import torch, sklearn; print('torch', torch.__version__, 'sklearn', sklearn.__version__)" 2>/dev/null
echo "== Singularity"; singularity --version 2>/dev/null || /cm/shared/apps_chop/singularity/4.4.1/bin/singularity --version
