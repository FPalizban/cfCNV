#!/bin/bash
# 5_resume_stage1.sh - push stage 1 forward with whatever is ready (login node), any time:
#   bash 5_resume_stage1.sh
# * aligns every sample whose FASTQ files are all downloaded and that is not aligned yet
#   (no waiting for the slowest download)
# * if no download job is in the queue any more, resubmits unfinished downloads
#   (wget resumes partial files)
# * then derive + task lists (finished samples are skipped in seconds)
# Safe to run as often as you like; run it again after the last download finishes.
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/phase3/full/p3_env.sh
cd "$S2"
join() { local IFS=,; echo "$*"; }
queued() { squeue -u "$USER" -h -o "%j" | grep -cx "$1" || true; }
if [[ $(queued P3_download) -eq 0 ]]; then
  dl=(); i=0
  while IFS=$'\t' read -r S RUN URL MD5 LOCAL; do [[ -f "$LOCAL.ok" ]] || dl+=("$i"); i=$((i+1)); done \
    < <(tail -n +2 "$S2_META/tasks_download.tsv")
  [[ ${#dl[@]} -gt 0 ]] && echo "downloads resubmitted (${dl[*]}): $(sbatch --parsable --array=$(join "${dl[@]}")%4 10_download.slurm)"
else echo "downloads still in the queue: $(queued P3_download) job(s)"; fi
[[ $(queued P3_align) -gt 0 ]] && { echo "alignment jobs already queued - run this again when they finish"; exit 0; }
al=(); wait_dl=0; i=0
while IFS=$'\t' read -r S FILES LAYOUT FRAC; do
  if [[ ! -f "$P3_BAM/$S/.align_done" ]]; then
    ok=1; IFS=';' read -ra F <<< "$FILES"; for f in "${F[@]}"; do [[ -f "$f.ok" ]] || ok=0; done
    [[ $ok -eq 1 ]] && al+=("$i") || wait_dl=$((wait_dl+1))
  fi; i=$((i+1))
done < <(tail -n +2 "$S2_META/tasks_align.tsv")
echo "samples ready to align now: ${#al[@]}   still waiting for downloads: $wait_dl"
dep=""
if [[ ${#al[@]} -gt 0 ]]; then a=$(sbatch --parsable --array=$(join "${al[@]}")%6 11_align.slurm); dep="--dependency=afterany:$a"; echo "align -> $a"; fi
n=$(( $(wc -l < "$S2_META/tasks_align.tsv") - 2 ))
v=$(sbatch --parsable $dep --array=0-$n 12_derive.slurm); echo "derive=$v"
# the tool task lists are rebuilt only when no caller job is queued: rebuilding them while stage 2
# runs would shift the rows that the queued jobs read
s2=$(squeue -u "$USER" -h -o "%j" | grep -cE '^P3_(qdnaseq|cnvkit|ichor|spectre|ace|bic_|fragle)' || true)
if [[ $s2 -eq 0 ]]; then echo "task lists: $(sbatch --parsable --dependency=afterany:$v 13_tasks.slurm)"
else echo "stage 2 is still running ($s2 jobs): task lists NOT rebuilt - run this script again after stage 2 has finished"; fi
[[ $wait_dl -gt 0 ]] && echo "NOTE: $wait_dl sample(s) still downloading - run bash 5_resume_stage1.sh again after the downloads finish."
