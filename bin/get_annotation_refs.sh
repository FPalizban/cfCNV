#!/usr/bin/env bash
# Download the public references used to annotate real-data calls (GRCh37), then condense them into the
# small tables read by annotate_calls.py and consensus.py (--tcga). Needs internet; ~10-20 min, ~1.5 GB raw.
#   TCGA PanCanAtlas arm calls (Taylor et al. 2018) + SNP6 segments, GENCODE v19 genes, OncoKB cancer genes,
#   gnomAD-SV v2.1 (AF >= 1%), UCSC segmental duplications, ENCODE blacklist v2, UCSC cytobands
#   + the UCSC hg38ToHg19 chain (to annotate GRCh38 calls)
# usage: bash get_annotation_refs.sh <out_dir>        -> <out_dir>/prepared  (pass this as --annotation_refs)
set -uo pipefail
OUT=${1:?usage: get_annotation_refs.sh <out_dir>}; HERE=$(cd "$(dirname "$0")" && pwd)
RAW="$OUT/raw"; mkdir -p "$RAW" "$OUT/prepared"; cd "$RAW"; bad=0
get() {
  if [[ -s "$1" ]]; then echo "present  $1"; return; fi
  if command -v wget >/dev/null; then wget -q -c --tries=5 -O "$1" "$2"; else curl -fsSL --retry 5 -o "$1" "$2"; fi
  [[ -s "$1" ]] && echo "OK       $1" || { echo "FAILED   $1  <- $2"; rm -f "$1"; bad=1; }
}
get tcga_arm_calls.txt        https://api.gdc.cancer.gov/data/4c35f34f-b0f3-4891-8794-4840dd748aad
get tcga_snp6_whitelisted.seg https://api.gdc.cancer.gov/data/00a32f7a-c85f-4f86-850d-be53973cbc4d
get gencode.v19.gtf.gz        https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
get oncokb_cancer_genes.json  https://www.oncokb.org/api/v1/utils/cancerGeneList
get gnomad_v2.1_sv.bed.gz     https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.bed.gz
get hg19_segdups.txt.gz       https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz
get hg19_blacklist.v2.bed.gz  https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz
get hg19_cytoband.txt.gz      https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
get hg38ToHg19.over.chain.gz  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
python3 "$HERE/prepare_annotation_refs.py" "$RAW" "$OUT/prepared" || bad=1
cp -f "$RAW/hg38ToHg19.over.chain.gz" "$OUT/prepared/" 2>/dev/null || true
ls -lh "$OUT/prepared"
[[ $bad -eq 0 ]] && echo REFS_OK || echo "some downloads FAILED - re-run to resume"
