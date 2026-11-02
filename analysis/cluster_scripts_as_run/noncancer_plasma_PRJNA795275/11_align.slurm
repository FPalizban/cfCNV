#!/bin/bash
#SBATCH --job-name=P3_align
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/phase3/align_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/phase3/align_%A_%a.err
# 11  One sample: reads (subsampled to ~12x for deep samples, mates kept together) -> BWA-MEM
#     (GRCh38 GIABv3) -> fixmate -> sort -> duplicates removed -> index -> samtools stats
#     (fragment-length histogram) -> mosdepth depth summary.
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/phase3/full/p3_env.sh
IFS=$'\t' read -r S FILES LAYOUT FRAC < <(task_line "$S2_META/tasks_align.tsv" "$SLURM_ARRAY_TASK_ID")
OD="$P3_BAM/$S"; mkdir -p "$OD"; BAM="$OD/$S.bam"
[[ -f "$OD/.align_done" ]] && { echo "SKIP $S"; exit 0; }
IFS=';' read -ra F <<< "$FILES"
for f in "${F[@]}"; do [[ -f "$f.ok" ]] || { echo "not downloaded yet: $f"; exit 3; }; done
aln_modules
TMP="$OD/tmp"; mkdir -p "$TMP"
RG="@RG\tID:$S\tSM:$S\tPL:ILLUMINA"
if [[ "$LAYOUT" == "PAIRED" ]]; then
  "$SYSPY" "$S2/fq_stream.py" "$FRAC" 11 PAIRED "${F[@]}" \
   | bwa mem -t 16 -p -Y -K 100000000 -R "$RG" "$REF38" - 2> "$OD/$S.bwa.log" \
   | samtools fixmate -m -u -@ 4 - - \
   | samtools sort -u -@ 8 -m 1500M -T "$TMP/s" - \
   | samtools markdup -r -@ 8 -f "$OD/$S.markdup.txt" - "$BAM"
else
  "$SYSPY" "$S2/fq_stream.py" "$FRAC" 11 SINGLE "${F[@]}" \
   | bwa mem -t 16 -Y -K 100000000 -R "$RG" "$REF38" - 2> "$OD/$S.bwa.log" \
   | samtools sort -u -@ 8 -m 1500M -T "$TMP/s" - \
   | samtools markdup -r -@ 8 -f "$OD/$S.markdup.txt" - "$BAM"
fi
samtools index -@ 8 "$BAM"
samtools stats -@ 8 "$BAM" > "$OD/$S.stats.txt"
"$MOSDEPTH" -t 4 -n --fast-mode "$OD/$S" "$BAM"
rm -rf "$TMP"; touch "$OD/.align_done"; echo "ALIGN_DONE $S"
