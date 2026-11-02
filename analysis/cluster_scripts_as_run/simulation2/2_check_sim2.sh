#!/bin/bash
# 2_check_sim2.sh - progress + first errors (login node, any time)
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/simulation2/full/s2_env.sh
c() { find "$S2_OUT/$1" -name "$2" 2>/dev/null | wc -l; }
printf "%-24s %5s / %s\n" "QDNAseq" "$(c QDNAseq .done)" 624 "CNVkit" "$(($(c CNVkit .done)))" 624 \
  "ichorCNA" "$(c ichorCNA .done)" 312 "Spectre" "$(c Spectre .done)" 624 "ACE" "$(c ACE .done)" 312 \
  "dudeML features" "$(c dudeML/features .features_done)" 104 "dudeML predictions" "$(c dudeML/predictions .done)" 624 \
  "CNV-MEANN features" "$(c CNV-MEANN/features_v2 .done)" 624 "CNV-MEANN applied" "$(c CNV-MEANN/calls .done)" 6 \
  "BIC-seq2 resources" "$([[ -f $BICRES38/RESOURCES_READY ]] && echo 1 || echo 0)" 1 \
  "BIC-seq2 unique" "$(find "$S2_OUT/BICseq2/unique" -name .done 2>/dev/null | wc -l)" 104 \
  "BIC-seq2 segmented" "$(find "$S2_OUT/BICseq2" -path '*bin_*' -name .done 2>/dev/null | wc -l)" 624 \
  "LYCEUM" "$(c LYCEUM .done)" 104
echo "scores: $([[ -s $S2_MET/SCORING_SUMMARY.txt ]] && echo ready || echo not yet)"
echo; for j in $(tr ' ' '\n' < "$S2_LOG/sim2_jobids.txt" | cut -d= -f2); do
  sacct -j "$j" -X -n -P -o JobName,State 2>/dev/null; done | sort | uniq -c
echo; for f in $(ls -t "$S2_LOG"/*.err 2>/dev/null | xargs grep -l -iE "error|traceback|killed|cannot" 2>/dev/null | head -4); do
  echo "--- $f"; grep -iE "error|traceback|killed|cannot" "$f" | grep -v "^\s*$" | head -3; done
