#!/bin/bash
#SBATCH --job-name=P3_bic_res
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/phase3/bic_resources_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/phase3/bic_resources_%j.err
# =====================================================================
# B0  BIC-seq2 resources for GRCh38 (GIABv3, autosomes) and 150-bp reads - run ONCE (long).
# Same procedure as the Simulation 1 GRCh37 101-mer resources (61_prepare + 61c fix):
# per-chromosome FASTA, GEM 150-mer mappability (2 mismatches), repaired WIG directives,
# bigWig -> bedGraph -> per-chromosome lists of uniquely mappable positions.
# =====================================================================
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/phase3/full/p3_env.sh
module load SAMtools/1.16.1-GCC-11.3.0
WRAP="$TOOLS/bin/run_bicseq2"; R="$BICRES38"; L=$BIC_READLEN
MER="GRCh38_GIABv3.${L}mer"; MAPROOT="$R/mappability"; MAPDIR="$MAPROOT/$MER"
mkdir -p "$R/reference" "$MAPDIR"
[[ -f "$R/RESOURCES_READY" ]] && { echo "already ready"; exit 0; }
printf 'chr%s\n' {1..22} > "$R/chromosomes.txt"
while read -r C; do [[ -s "$R/reference/$C.fa" ]] || samtools faidx "$REF38" "$C" > "$R/reference/$C.fa"; done < "$R/chromosomes.txt"
cut -f1,2 "$REF38_AUTO.fai" > "$MAPROOT/$MER.sizes"
if [[ ! -s "$MAPROOT/$MER.wig" ]]; then
  "$WRAP" bash -s -- "$REF38_AUTO" "$MAPROOT" "$MER" "$L" <<'INNER'
set -euo pipefail
REF="$1"; MAPROOT="$2"; MER="$3"; L="$4"; cd "$MAPROOT"
for X in gem-indexer gem-mappability gem-2-wig wigToBigWig bigWigToBedGraph; do command -v "$X"; done
[[ -s GRCh38_auto.gem ]] || gem-indexer -i "$REF" -o GRCh38_auto -T 24
gem-mappability -m 2 -I GRCh38_auto.gem -l "$L" -o "$MER" -T 24 2> "$MER.mappability.log"
gem-2-wig -I GRCh38_auto.gem -i "$MER.mappability" -o "$MER"
INNER
fi
# repair the WIG header lines exactly as in Simulation 1 (61c): keep only key=value fields
awk '/^(variableStep|fixedStep)[[:space:]]/{o=$1; for(i=2;i<=NF;i++) if(index($i,"=")>0) o=o" "$i; print o; next}{print}' \
    "$MAPROOT/$MER.wig" > "$MAPROOT/$MER.fixed.wig"
"$WRAP" wigToBigWig "$MAPROOT/$MER.fixed.wig" "$MAPROOT/$MER.sizes" "$MAPROOT/$MER.bw"
"$WRAP" bigWigToBedGraph "$MAPROOT/$MER.bw" "$MAPROOT/$MER.bedGraph"
rm -f "$MAPDIR"/*.txt
awk -v out="$MAPDIR" -v mer="$MER" 'BEGIN{OFS="\t"} $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ && $4==1 {print $2,$3 >> (out "/" mer "." $1 ".txt")}' \
    "$MAPROOT/$MER.bedGraph"
n=$(ls "$MAPDIR"/*.txt | wc -l); echo "mappability files: $n"; [[ $n -eq 22 ]]
rm -f "$MAPROOT/$MER.fixed.wig"
touch "$R/RESOURCES_READY"; echo "BIC_RESOURCES_READY $R"
