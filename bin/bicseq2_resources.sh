#!/usr/bin/env bash
# BIC-seq2 resources, as built for the cfCNV paper (run once per reference and read length; ~1-2 days):
# per-chromosome FASTA, GEM mappability (2 mismatches), repaired WIG directives, bigWig -> bedGraph ->
# per-chromosome lists of uniquely mappable positions.
# usage: bicseq2_resources.sh <ref.fa> <read_length> <out_dir> <threads>
#   -> <out_dir>/reference/<chr>.fa  and  <out_dir>/mappability/<prefix>/<prefix>.<chr>.txt  (prefix = <out_dir basename>)
set -euo pipefail
REF=$1; L=$2; OUT=$3; T=${4:-8}
MER=$(basename "$OUT"); MAPDIR="$OUT/mappability/$MER"; mkdir -p "$OUT/reference" "$MAPDIR"
CHRS=$(awk '{c=$1; sub(/^chr/,"",c)} c ~ /^([1-9]|1[0-9]|2[0-2])$/ {print $1}' "$REF.fai")
for C in $CHRS; do samtools faidx "$REF" "$C" > "$OUT/reference/$C.fa"; done
samtools faidx "$REF" $CHRS > auto.fa; samtools faidx auto.fa; cut -f1,2 auto.fa.fai > "$MER.sizes"
gem-indexer -i auto.fa -o auto_index -T "$T"
gem-mappability -m 2 -I auto_index.gem -l "$L" -o "$MER" -T "$T"
gem-2-wig -I auto_index.gem -i "$MER.mappability" -o "$MER"
awk '/^(variableStep|fixedStep)[[:space:]]/{o=$1; for(i=2;i<=NF;i++) if(index($i,"=")>0) o=o" "$i; print o; next}{print}' \
    "$MER.wig" > "$MER.fixed.wig"
wigToBigWig "$MER.fixed.wig" "$MER.sizes" "$MER.bw"
bigWigToBedGraph "$MER.bw" "$MER.bedGraph"
awk -v out="$MAPDIR" -v mer="$MER" 'BEGIN{OFS="\t"} $4==1 {print $2,$3 >> (out "/" mer "." $1 ".txt")}' "$MER.bedGraph"
n=$(ls "$MAPDIR"/*.txt | wc -l); echo "mappability files: $n"; [[ $n -eq 22 ]]
rm -f auto.fa* auto_index* "$MER".*wig "$MER.mappability" "$MER.bw" "$MER.bedGraph"
touch "$OUT/RESOURCES_READY"
