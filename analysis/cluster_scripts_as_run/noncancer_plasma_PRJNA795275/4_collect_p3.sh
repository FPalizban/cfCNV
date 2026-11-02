#!/bin/bash
# 4_collect_p3.sh - after stage 2 (login node).  -> ~/cfcnv_phase3_for_claude.tar.gz
# Inventory, QC, fragment-length histograms and the callers' segment/call files (no BAMs).
set -uo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/phase3/full/p3_env.sh
cd "$P"; L=$(mktemp); Q="$S2_MET/qc"; mkdir -p "$Q"
for d in "$P3_BAM"/*/; do s=$(basename "$d")
  grep -P '^IS\t' "$d/$s.stats.txt" 2>/dev/null | cut -f2-3 > "$Q/$s.insert_sizes.tsv"
  cp "$d/$s.mosdepth.summary.txt" "$d/$s.markdup.txt" "$d/derived.tsv" "$Q/" 2>/dev/null
  [[ -f "$Q/derived.tsv" ]] && mv "$Q/derived.tsv" "$Q/$s.derived.tsv"
done
ls 00_metadata/phase3/*.tsv 00_metadata/phase3/*.txt 2>/dev/null >> "$L"
find 09_metrics/phase3 -type f >> "$L"
O=07_tool_outputs/phase3
find $O/QDNAseq -name '*.segments.igv' -o -name '*.calls.vcf' >> "$L"
find $O/CNVkit -name '*.cns' >> "$L"
find $O/ichorCNA -name '*.seg' -o -name '*.seg.txt' -o -name '*.params.txt' -o -name '*.cna.seg' >> "$L"
find $O/Spectre -path '*/calls/*' \( -name '*.vcf*' -o -name '*_cnv.bed.gz' \) >> "$L"
find $O/BICseq2 -name '*.cnv' >> "$L"
find $O/ACE -name '*.ace.tsv' >> "$L"
find $O/Fragle -name '*.csv' >> "$L"
for p in download align derive tasks qdnaseq cnvkit ichor spectre ace bic_resources bic_unique bic_normseg fragle; do
  ls -t logs/phase3/${p}_*.err 2>/dev/null | head -1 >> "$L"; ls -t logs/phase3/${p}_*.out 2>/dev/null | head -1 >> "$L"
done
tar -czf "$HOME/cfcnv_phase3_for_claude.tar.gz" -T "$L" 2>/dev/null
echo "$(wc -l < "$L") files -> $HOME/cfcnv_phase3_for_claude.tar.gz ($(du -h "$HOME/cfcnv_phase3_for_claude.tar.gz" | cut -f1))"
rm -f "$L"
