#!/bin/bash
# =====================================================================
# 2_submit_stage2.sh - Phase 3 stage 2, the callers (login node), AFTER stage 1 has finished
# (3_check_p3.sh says "stage 1: complete"):   bash 2_submit_stage2.sh
#   A1 QDNAseq | A2 CNVkit | A3 ichorCNA | A4 Spectre | A5 ACE (after A1)
#   B0 BIC-seq2 mappability for this read length (only if missing; up to ~2 days) -> B1 -> B2
#   F1 Fragle (paired-end BAMs)
# =====================================================================
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/phase3/full/p3_env.sh
cd "$S2"
[[ -f "$S2_META/STAGE1_COMPLETE" ]] || { echo "Stage 1 is not complete yet - run bash 3_check_p3.sh"; exit 1; }
last() { echo $(( $(wc -l < "$S2_META/tasks_$1.tsv") - 2 )); }
sub() { sbatch --parsable "$@"; }
echo "read length $BIC_READLEN -> BIC-seq2 resources $BICRES38"
a1=$(sub --array=0-$(last qdnaseq)%40 A1_qdnaseq.slurm)
a2=$(sub --array=0-$(last cnvkit)%40 A2_cnvkit.slurm)
a3=$(sub --array=0-$(last ichor)%10 A3_ichorcna.slurm)
a4=$(sub --array=0-$(last spectre)%30 A4_spectre.slurm)
a5=$(sub --dependency=afterany:$a1 --array=0-$(last ace)%30 A5_ace.slurm)
b1=$(sub --array=0-$(last bic_unique)%20 B1_bicseq2_unique.slurm)
if [[ -f "$BICRES38/RESOURCES_READY" ]]; then
  b0=ready; b2=$(sub --dependency=afterany:$b1 --array=0-$(last bic_norm)%40 B2_bicseq2_normseg.slurm)
else
  b0=$(sub B0_bicseq2_resources.slurm)
  b2=$(sub --dependency=afterok:$b0,afterany:$b1 --array=0-$(last bic_norm)%40 B2_bicseq2_normseg.slurm)
fi
f1=none
if [[ $(last fragle) -ge 0 && -s "$FRAGLE_SRC/main.py" ]]; then f1=$(sub --array=0-$(last fragle)%10 F1_fragle.slurm); fi
echo "qdnaseq=$a1 cnvkit=$a2 ichor=$a3 spectre=$a4 ace=$a5 bic_res=$b0 bic_unique=$b1 bic_seg=$b2 fragle=$f1" | tee "$S2_LOG/stage2_jobids.txt"
echo "Submitted stage 2. Check with: bash 3_check_p3.sh"
