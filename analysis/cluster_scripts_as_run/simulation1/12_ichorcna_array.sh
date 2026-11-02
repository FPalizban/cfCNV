#!/bin/bash
#SBATCH --job-name=ichor_s1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=06:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/ichor_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/ichor_%A_%a.err
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/00_env.sh
line=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+2)) 'NR==n{print}' "$SIM1_MAN/ichor_tasks.tsv")
IFS=$'\t' read -r sample bam bg depth rep cls bin <<< "$line"
module load R/4.2.2-foss-2022b
module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
export R_LIBS_USER="$TOOLS/Rlib"
od="$SIM1_OUT/ichorCNA/$sample/1000000"; mkdir -p "$od"
gc="$ICHOR_EXT/gc_hg19_1000kb.wig"; map="$ICHOR_EXT/map_hg19_1000kb.wig"
cent="$ICHOR_EXT/GRCh37.p13_centromere_UCSC-gapTable.txt"; pon="$ICHOR_EXT/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds"
for f in "$gc" "$map" "$cent"; do [[ -s "$f" ]]||{ echo "missing $f";exit 4;};done
wig="$od/$sample.wig"
"$READCOUNTER" --window 1000000 --quality 20 --chromosome "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y" "$bam" > "$wig"
cmd=(Rscript "$ICHOR_DIR/scripts/runIchorCNA.R" --id "$sample" --WIG "$wig" --gcWig "$gc" --mapWig "$map" --centromere "$cent"
 --ploidy 'c(2,3)' --normal 'c(0.5,0.6,0.7,0.8,0.9,0.95,0.99)' --maxCN 5 --includeHOMD False
 --chrs 'c(1:22,"X")' --chrTrain 'c(1:22)' --estimateNormal True --estimatePloidy True --estimateScPrevalence True
 --scStates 'c(1,3)' --txnE 0.9999 --txnStrength 10000 --plotFileType png --plotYLim 'c(-2,4)' --outDir "$od")
[[ -s "$pon" ]] && cmd+=(--normalPanel "$pon")
/usr/bin/time -v "${cmd[@]}"
seg=$(find "$od" -maxdepth 1 -type f \( -name '*seg*.txt' -o -name '*segments*.txt' \)|head -1)
[[ -s "$seg" ]] || { find "$od" -maxdepth 1 -type f -print;exit 5; }
python3 - "$sample" "$seg" "$od/$sample.calls.tsv" <<'PY'
import csv,sys,re
s,i,o=sys.argv[1:]
with open(i) as f:h=f.readline()
d="\t" if "\t" in h else None
with open(i) as f:
 if d: rows=list(csv.DictReader(f,delimiter="\t"))
 else:
  hdr=re.split(r"\s+",f.readline().strip());rows=[dict(zip(hdr,re.split(r"\s+",l.strip()))) for l in f if l.strip()]
def key(x,pat):
 for k in x:
  if re.search(pat,k,re.I):return k
with open(o,"w") as g:
 w=csv.writer(g,delimiter="\t",lineterminator="\n");w.writerow(["sample_id","chrom","start","end","cnv_type","copy_number","segment_value"])
 for x in rows:
  ck=key(x,"chr|chrom");sk=key(x,"start");ek=key(x,"end");nk=key(x,"copy|cn")
  if not all([ck,sk,ek,nk]):continue
  try:n=float(x[nk])
  except:continue
  t="LOSS" if n<2 else ("GAIN" if n>2 else "NEUTRAL")
  if t!="NEUTRAL":w.writerow([s,x[ck],x[sk],x[ek],t,n,""])
PY
