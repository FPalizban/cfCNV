#!/bin/bash
# 2_check.sh - progress of the ML dataset jobs (read only):   bash 2_check.sh
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/ml_dataset/ml_env.sh
echo "== job ids:  $(cat "$MLLOG/ml_jobids.txt" 2>/dev/null)"
echo "== still in the queue:"; squeue -u "$USER" -h -o "%.10i %.16j %.3t %.10M %R" | grep -E "ML_" || echo "   none"
N=$(( $(wc -l < "$MLO/manifest.tsv") - 1 ))
done_n=$(find "$MLO/counts_10kb" -name "*.done" 2>/dev/null | wc -l)
echo "== counts finished: $done_n / $N BAMs"
echo "== reference annotation: $(ls "$MLO"/reference/*.tsv.gz 2>/dev/null | wc -l) / 2 builds"
[[ -f "$MLO/hf_dataset/README.md" ]] && echo "== dataset: DONE ($MLO/hf_dataset)" || echo "== dataset: not yet"
[[ -f "$MLO/ml_results_for_claude.tar.gz" ]] && echo "== baseline: DONE -> $MLO/ml_results_for_claude.tar.gz" || echo "== baseline: not yet"
echo "== errors in logs (last lines of any log that has 'Error' or 'Traceback'):"
found=0
for f in $(grep -l -E "Traceback|Error|error:|CANCELLED|oom" "$MLLOG"/*.err 2>/dev/null | head -5); do
  found=1; echo "--- $f"; tail -5 "$f"
done
[[ $found == 0 ]] && echo "   none"
