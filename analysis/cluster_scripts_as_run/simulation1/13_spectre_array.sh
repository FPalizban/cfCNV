#!/bin/bash
#SBATCH --job-name=spectre_s1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=06:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/spectre_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/spectre_%A_%a.err
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/00_env.sh
line=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+2)) 'NR==n{print}' "$SIM1_MAN/spectre_tasks.tsv")
IFS=$'\t' read -r sample bam bg depth rep cls <<< "$line"
od="$SIM1_OUT/Spectre/$sample";mkdir -p "$od/mosdepth" "$od/calls"
p="$od/mosdepth/$sample"
"$MOSDEPTH" -t 8 -x -b 1000 -Q 20 "$p" "$bam"
cov="$p.regions.bed.gz";[[ -s "$cov" ]]||exit 4
/usr/bin/time -v "$SPECTRE" CNVCaller --coverage "$cov" --sample-id "$sample" --output-dir "$od/calls" --reference "$REF37_FASTA" --min-cnv-len 50000 --cancer --threads 8
bed=$(find "$od/calls" -type f -name '*.bed'|head -1);[[ -s "$bed" ]]||exit 5
python3 - "$sample" "$bed" "$od/$sample.calls.tsv" <<'PY'
import csv,sys
s,i,o=sys.argv[1:]
with open(o,"w") as g:
 w=csv.writer(g,delimiter="\t",lineterminator="\n");w.writerow(["sample_id","chrom","start","end","cnv_type","copy_number","segment_value"])
 for l in open(i):
  if not l.strip() or l.startswith("#"):continue
  p=l.rstrip().split("\t");txt="\t".join(p[3:]).upper()
  t="LOSS" if ("DEL" in txt or "LOSS" in txt) else ("GAIN" if ("DUP" in txt or "GAIN" in txt) else None)
  if t:w.writerow([s,p[0],p[1],p[2],t,"",""])
PY
