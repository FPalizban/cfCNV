#!/bin/bash
# 3_collect_sim2.sh - after the jobs finished (login node).  -> ~/cfcnv_sim2_for_claude.tar.gz
# Scores, per-mixture summaries, ichorCNA/ACE estimates and a few logs. No BAMs, no large calls.
set -uo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/simulation2/full/s2_env.sh
cd "$P"; L=$(mktemp)
find 09_metrics/simulation2 -type f > "$L"
ls 00_metadata/simulation2/HG008_mixtures_realised.tsv 00_metadata/simulation2/pool_stats.tsv logs/simulation2/sim2_jobids.txt >> "$L" 2>/dev/null
find 07_tool_outputs/simulation2/ichorCNA -name '*.params.txt' >> "$L"
find 07_tool_outputs/simulation2/ACE -name '*.ace.tsv' >> "$L"
find 07_tool_outputs/simulation2/CNV-MEANN/calls -name 'prediction_log.tsv' >> "$L"
ls -t logs/simulation2/score_*.out 2>/dev/null | head -2 >> "$L"
for p in qdnaseq cnvkit ichor spectre ace dude_feat dude_pred meann_feat meann_apply bic_resources bic_unique bic_normseg lyceum; do
  ls -t logs/simulation2/${p}_*.err 2>/dev/null | head -1 >> "$L"; ls -t logs/simulation2/${p}_*.out 2>/dev/null | head -1 >> "$L"
done
tar -czf "$HOME/cfcnv_sim2_for_claude.tar.gz" -T "$L" 2>/dev/null
echo "$(wc -l < "$L") files -> $HOME/cfcnv_sim2_for_claude.tar.gz ($(du -h "$HOME/cfcnv_sim2_for_claude.tar.gz" | cut -f1))"
rm -f "$L"
