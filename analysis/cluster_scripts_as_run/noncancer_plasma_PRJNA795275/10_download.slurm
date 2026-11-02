#!/bin/bash
#SBATCH --job-name=P3_download
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/phase3/download_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/phase3/download_%A_%a.err
# 10  Download one FASTQ file from ENA (resumes partial files, checks md5, 5 attempts).
#     Runs without an ENA copy fall back to the SRA toolkit (prefetch + fasterq-dump) if available.
set -uo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/phase3/full/p3_env.sh
IFS=$'\t' read -r S RUN URL MD5 LOCAL < <(task_line "$S2_META/tasks_download.tsv" "$SLURM_ARRAY_TASK_ID")
mkdir -p "$(dirname "$LOCAL")"
[[ -f "$LOCAL.ok" ]] && { echo "SKIP $LOCAL"; exit 0; }
if [[ "$URL" == "SRA" ]]; then
  [[ -s "$S2_META/p3_sra_module.sh" ]] && source "$S2_META/p3_sra_module.sh"
  command -v fasterq-dump >/dev/null || { echo "no ENA copy and no SRA toolkit for $RUN"; exit 3; }
  cd "$(dirname "$LOCAL")" && prefetch "$RUN" && fasterq-dump --split-files -e 1 "$RUN" && gzip -f "${RUN}"*.fastq \
    && touch "$LOCAL.ok" && exit 0; exit 4
fi
# Interrupted transfers are resumed (-c, unlimited retries) until the server reports the file
# complete; only a COMPLETE file with a wrong md5 is deleted and fetched again (3 times).
for i in 1 2 3; do
  until wget -c -q --tries=0 --waitretry=30 --read-timeout=120 --timeout=120 -P "$(dirname "$LOCAL")" "$URL"; do
    echo "transfer interrupted at $(du -h "$LOCAL" 2>/dev/null | cut -f1), resuming"; sleep 30; done
  if [[ -z "$MD5" ]] || [[ "$(md5sum "$LOCAL" | cut -d' ' -f1)" == "$MD5" ]]; then
    touch "$LOCAL.ok"; echo "DOWNLOADED $LOCAL ($(du -h "$LOCAL" | cut -f1))"; exit 0
  fi
  echo "complete file but md5 mismatch (attempt $i) - fetching again"; rm -f "$LOCAL"; sleep 60
done
echo "FAILED $URL"; exit 2
