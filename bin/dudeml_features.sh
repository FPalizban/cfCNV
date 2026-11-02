#!/usr/bin/env bash
# dudeML window statistics and feature vectors for one BAM, as in the cfCNV paper (5 flanking windows).
# usage: dudeml_features.sh <bam> <sample> <bin_bp> <truth.tsv or NONE>
#   -> <sample>.sample.bed (for prediction) and, with truth, <sample>.train.bed (labelled, for training)
set -euo pipefail
BAM=$1; S=$2; B=$3; TRUTH=${4:-NONE}; D=${DUDEML_DIR:-/opt/dudeML}
samtools idxstats "$BAM" | awk '{c=$1; sub(/^chr/,"",c)} c ~ /^([1-9]|1[0-9]|2[0-2])$/ {print $1}' > chroms.txt
dudeml_winstat.py --bam "$BAM" --sample "$S" --outdir . --chroms chroms.txt --bins "$B"
mkdir -p fv
python "$D/dudeML.py" fvecSample -i "$S.$B.winStat.bed" -o "$S.sample.bed" -id "$S" -d fv -w "$B" -s "$B" -windows 5
test -s "$S.sample.bed"
if [[ "$TRUTH" != "NONE" && -s "$TRUTH" ]]; then
  P=$(head -n1 chroms.txt | grep -q '^chr' && echo chr || true)
  awk -v p="$P" 'BEGIN{OFS="\t"} !/^#/ {c=$1; sub(/^chr/,"",c); if(tolower($4)=="gain") print p c,$2,$3,"dup",0.5,2}' "$TRUTH" > dups.bed
  awk -v p="$P" 'BEGIN{OFS="\t"} !/^#/ {c=$1; sub(/^chr/,"",c); if(tolower($4)=="loss") print p c,$2,$3,"del",0.5,0}' "$TRUTH" > dels.bed
  python "$D/dudeML.py" fvecTrain -i "$S.$B.winStat.bed" -o "$S.train.bed" -d fv -w "$B" -s "$B" -windows 5 -dups dups.bed -dels dels.bed
  test -s "$S.train.bed"
fi
rm -rf fv "$S.$B.winStat.bed"
