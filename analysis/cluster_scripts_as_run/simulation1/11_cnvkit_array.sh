#!/bin/bash
#SBATCH --job-name=cnvkit_s1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=06:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/cnvkit_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/cnvkit_%A_%a.err
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/00_env.sh
line=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+2)) 'NR==n{print}' "$SIM1_MAN/cnvkit_tasks.tsv")
IFS=$'\t' read -r sample bam bg depth rep cls bin <<< "$line"
module load CNVkit/0.9.10-foss-2022b-R-4.2.2
od="$SIM1_OUT/CNVkit/$sample/$bin"; mkdir -p "$od"
/usr/bin/time -v cnvkit.py batch "$bam" -m wgs -n -f "$REF37_FASTA" --target-avg-size "$bin" -d "$od" -p 4
cns=$(find "$od" -maxdepth 1 -name '*.cns'|head -1); [[ -s "$cns" ]] || exit 3
python3 - "$sample" "$cns" "$od/$sample.calls.tsv" <<'PY'
import csv,sys
s,i,o=sys.argv[1:]
with open(i) as f,open(o,"w") as g:
 r=csv.DictReader(f,delimiter="\t");w=csv.writer(g,delimiter="\t",lineterminator="\n")
 w.writerow(["sample_id","chrom","start","end","cnv_type","copy_number","segment_value"])
 for x in r:
  v=float(x["log2"]);t="LOSS" if v<-.2 else ("GAIN" if v>.2 else "NEUTRAL")
  if t!="NEUTRAL":w.writerow([s,x["chromosome"],x["start"],x["end"],t,"",x["log2"]])
PY
