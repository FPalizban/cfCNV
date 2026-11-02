#!/bin/bash
#SBATCH --job-name=G7_bic_seg
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/gse71378/bic_normseg_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/gse71378/bic_normseg_%A_%a.err
# B2  BIC-seq2 normalisation + segmentation, same parameters as Simulation 1 (lambda 2),
#     read length 150, 6 bin sizes (624 tasks).
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/gse71378/full/p3_env.sh
WRAP="$TOOLS/bin/run_bicseq2"; R="$BICRES38"; MER="human_g1k_v37.${BIC_READLEN}mer"
IFS=$'\t' read -r S BAM D B < <(task_line "$S2_META/tasks_bic_norm.tsv" "$SLURM_ARRAY_TASK_ID")
U="$S2_OUT/BICseq2/unique/$S"; OD="$S2_OUT/BICseq2/$S/bin_$B"; mkdir -p "$OD"/{norm,seg/tmp,tmp}
[[ -f "$OD/.done" ]] && { echo "SKIP $S $B"; exit 0; }
[[ -f "$R/RESOURCES_READY" && -f "$U/.done" && -s "$U/fragment_size.txt" ]] || { echo "inputs not ready"; exit 3; }
NCFG="$OD/$S.norm-config.tsv"; SCFG="$OD/$S.seg-config.tsv"
printf "chromName\tfaFile\tMapFile\treadPosFile\tbinFileNorm\n" > "$NCFG"; printf "chromName\tbinFileNorm\n" > "$SCFG"
for C in $(seq 1 22); do
  printf "%s\t%s\t%s\t%s\t%s\n" "$C" "$R/reference/$C.fa" "$R/mappability/$MER/$MER.$C.txt" \
         "$U/$S.chr$C.seq" "$OD/norm/$S.chr$C.norm.bin" >> "$NCFG"
  printf "%s\t%s\n" "$C" "$OD/norm/$S.chr$C.norm.bin" >> "$SCFG"
done
"$WRAP" perl /NBICseq-norm_v0.2.4/NBICseq-norm.pl --tmp="$OD/tmp" -l "$BIC_READLEN" -s "$(cat "$U/fragment_size.txt")" \
   -b "$B" --fig "$OD/$S.GC.pdf" "$NCFG" "$OD/$S.norm_parameters.txt"
"$WRAP" perl /NBICseq-seg_v0.7.2/NBICseq-seg.pl --lambda=2 --tmp="$OD/seg/tmp" --fig="$OD/seg/$S.png" "$SCFG" "$OD/seg/$S.cnv"
test -f "$OD/seg/$S.cnv"; rm -rf "$OD/tmp" "$OD/seg/tmp"
touch "$OD/.done"; echo "BIC_SEG_DONE $S $B"
