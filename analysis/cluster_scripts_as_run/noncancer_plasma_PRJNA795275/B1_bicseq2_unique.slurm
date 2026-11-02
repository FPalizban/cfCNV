#!/bin/bash
#SBATCH --job-name=P3_bic_uniq
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/phase3/bic_unique_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/phase3/bic_unique_%A_%a.err
# B1  BIC-seq2 uniquely-mapped read positions per chromosome + median fragment size,
#     identical to Simulation 1 (62_bicseq2_unique). One task per mixture (104 tasks).
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/phase3/full/p3_env.sh
WRAP="$TOOLS/bin/run_bicseq2"
IFS=$'\t' read -r S BAM D X < <(task_line "$S2_META/tasks_bic_unique.tsv" "$SLURM_ARRAY_TASK_ID")
OUT="$S2_OUT/BICseq2/unique/$S"; mkdir -p "$OUT"
[[ -f "$OUT/.done" ]] && { echo "SKIP $S"; exit 0; }
"$WRAP" bash -c '
set -euo pipefail
BAM="$1"; S="$2"; OUT="$3"; H=/samtools-0.1.7a_getUnique-0.1.3/misc/samtools.pl; [[ -f "$H" ]]
for C in $(seq 1 22); do
  ( set -o pipefail; samtools view "$BAM" "chr$C" | perl "$H" unique - | cut -f4 > "$OUT/$S.chr$C.seq" ) &
  (( C % 4 == 0 )) && wait
done; wait' _ "$BAM" "$S" "$OUT"
for C in $(seq 1 22); do [[ -f "$OUT/$S.chr$C.seq" ]] || { echo "missing chr$C"; exit 20; }; done
"$WRAP" bash -c 'set +o pipefail; samtools view -f 2 -F 3852 "$1" | awk '"'"'$9>0 && $9<2000 {print $9; n++; if(n>=100000) exit}'"'"'' _ "$BAM" \
  | sort -n | awk '{a[NR]=$1} END{ if(NR<100) exit 1; print (NR%2)? a[(NR+1)/2] : int((a[NR/2]+a[NR/2+1])/2+0.5)}' > "$OUT/fragment_size.txt"
echo "fragment size: $(cat "$OUT/fragment_size.txt")"
touch "$OUT/.done"; echo "BIC_UNIQUE_DONE $S"
