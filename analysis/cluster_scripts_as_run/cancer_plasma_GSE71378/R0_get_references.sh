#!/bin/bash
# R0_get_references.sh - run ONCE on the LOGIN node (needs internet), ~10-20 min:  bash R0_get_references.sh
# Public GRCh37 references for annotating and evaluating the real-cancer calls; then condenses
# them into small tables (R1_prepare_refs.py). Must end with REFS_OK.
#   TCGA PanCanAtlas arm calls (Taylor et al. 2018) + SNP6 segments  -> per-cancer-type arm and 1-Mb
#                                                                       gain/loss frequencies
#   GENCODE v19 genes, OncoKB cancer genes (oncogene / tumour suppressor)
#   gnomAD-SV v2.1 common CNVs (AF >= 1%), UCSC segmental duplications, ENCODE blacklist v2,
#   UCSC cytobands (arm boundaries), MSigDB Hallmark gene sets
set -uo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/gse71378/full/p3_env.sh
R="$P/02_reference/annotation_GRCh37"; RAW="$R/raw"; mkdir -p "$RAW" "$R/prepared"; cd "$RAW"; bad=0
get() {  # get <file> <url>
  if [[ -s "$1" ]]; then echo "present  $1"; return; fi
  wget -q -c --tries=5 -O "$1" "$2" && [[ -s "$1" ]] && echo "OK       $1 ($(du -h "$1" | cut -f1))" || { echo "FAILED   $1  <- $2"; rm -f "$1"; bad=1; }
}
get tcga_arm_calls.txt        https://api.gdc.cancer.gov/data/4c35f34f-b0f3-4891-8794-4840dd748aad
get tcga_snp6_whitelisted.seg https://api.gdc.cancer.gov/data/00a32f7a-c85f-4f86-850d-be53973cbc4d
get gencode.v19.gtf.gz        https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
get oncokb_cancer_genes.json  https://www.oncokb.org/api/v1/utils/cancerGeneList
get gnomad_v2.1_sv.bed.gz     https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.bed.gz
get hg19_segdups.txt.gz       https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz
get hg19_blacklist.v2.bed.gz  https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz
get hg19_cytoband.txt.gz      https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
get hallmark.gmt              https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt
echo; echo "== condensing =="
"$SYSPY" "$S2/R1_prepare_refs.py" "$RAW" "$R/prepared" || bad=1
ls -lh "$R/prepared"
echo; [[ $bad -eq 0 ]] && echo REFS_OK || echo "REFS: some files FAILED - send me this output (the rest can still be used)"
