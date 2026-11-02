#!/bin/bash
# 3_check_p3.sh - progress and first errors (login node, any time)
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/gse71378/full/p3_env.sh
n() { [[ -s "$S2_META/tasks_$1.tsv" ]] && echo $(( $(wc -l < "$S2_META/tasks_$1.tsv") - 1 )) || echo "?"; }
c() { find "$S2_OUT/$1" -name "$2" 2>/dev/null | wc -l; }
printf "%-22s %5s / %s\n" \
  "BAMs ready" "$(find "$P3_BAM" -name .align_done | wc -l)" "$(grep -vcE '^#|^\s*$' "$S2/samples.txt")"
echo "stage 1: $([[ -f $S2_META/STAGE1_COMPLETE ]] && echo complete || echo running)"
if [[ -f "$S2_META/STAGE1_COMPLETE" ]]; then
  printf "%-22s %5s / %s\n" QDNAseq "$(c QDNAseq .done)" "$(n qdnaseq)" CNVkit "$(c CNVkit .done)" "$(n cnvkit)" \
    ichorCNA "$(c ichorCNA .done)" "$(n ichor)" Spectre "$(c Spectre .done)" "$(n spectre)" ACE "$(c ACE .done)" "$(n ace)" \
    "BIC-seq2 resources" "$([[ -f $BICRES38/RESOURCES_READY ]] && echo 1 || echo 0)" 1 \
    "BIC-seq2 unique" "$(find "$S2_OUT/BICseq2/unique" -name .done 2>/dev/null | wc -l)" "$(n bic_unique)" \
    "BIC-seq2 segmented" "$(find "$S2_OUT/BICseq2" -path '*bin_*' -name .done 2>/dev/null | wc -l)" "$(n bic_norm)" \
    Fragle "$(c Fragle .done)" "$(n fragle)"
fi
echo; squeue -u "$USER" -h -o "%.10i %.14j %.3t %.10M %R" | head -15
echo; for f in $(ls -t "$S2_LOG"/*.err 2>/dev/null | xargs grep -l -iE "error|traceback|killed|cannot|failed" 2>/dev/null | head -4); do
  echo "--- $f"; grep -iE "error|traceback|killed|cannot|failed" "$f" | head -3; done
