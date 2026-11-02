#!/bin/bash
# 0_slurm_usage.sh - run time, peak memory and state of every job since 1 Sep 2026 (Simulations 1
# and 2, rescue and checks), for Results 3.6. Login node, ~1 min:   bash 0_slurm_usage.sh
# -> ~/cfcnv_slurm_usage.tsv.gz  (one line per job step; MaxRSS is on the .batch step)
sacct -u "$USER" -S 2026-09-01 -E now -P \
  -o JobID,JobName%60,State,ExitCode,Elapsed,ElapsedRaw,TotalCPU,AllocCPUS,ReqMem,MaxRSS,MaxVMSize,Start,End,NodeList \
  | gzip > "$HOME/cfcnv_slurm_usage.tsv.gz"
echo "$(zcat "$HOME/cfcnv_slurm_usage.tsv.gz" | wc -l) lines -> $HOME/cfcnv_slurm_usage.tsv.gz"
