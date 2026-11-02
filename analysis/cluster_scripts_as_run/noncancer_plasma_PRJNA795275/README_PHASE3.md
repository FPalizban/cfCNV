# cfCNV - remaining results: run times/memory (3.6) and Phase 3 real cfDNA (3.7)

All commands: copy-paste one block at a time. Lines starting with # are comments.

## A. Run times and memory (3.6) - 1 minute, any time
```bash
cd /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/phase3/full
bash 0_slurm_usage.sh            # -> ~/cfcnv_slurm_usage.tsv.gz  (send it to Claude)
```

## B. Phase 3 (PRJNA795275)
1. Setup (login node, internet, ~10-20 min). Must end with SETUP_P3_OK. Paste the inventory table it prints to Claude.
   `bash 0_setup_p3.sh`
2. Stage 1 - download, align, depth series (1-2 days):  `bash 1_submit_stage1.sh`
3. Check any time:  `bash 3_check_p3.sh`   (stage 1 is done when it prints "stage 1: complete")
4. Stage 2 - callers:  `bash 2_submit_stage2.sh`
5. When 3_check_p3.sh shows all tools complete:  `bash 4_collect_p3.sh`  -> ~/cfcnv_phase3_for_claude.tar.gz

Re-submitting any step is safe: finished work is skipped.
Disk: FASTQ under 01_raw_data/PRJNA795275 and BAMs under 04_aligned/PRJNA795275.
The FASTQ files can be deleted after stage 1 if space is short (they can be downloaded again).
