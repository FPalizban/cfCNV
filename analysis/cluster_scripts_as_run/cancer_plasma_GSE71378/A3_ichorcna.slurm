#!/bin/bash
#SBATCH --job-name=G7_ichor
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/gse71378/ichor_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/gse71378/ichor_%A_%a.err
# A3  ichorCNA on GRCh38, three configurations (312 tasks):
#   1000kb_default = tumour defaults (normal 0.5-0.99, ploidy 2-3, subclonal states) - the usual setting
#   1000kb_lowTF   = ichorCNA's recommended low-tumour-fraction setting (normal 0.95-0.999)
#   500kb_default  = second resolution
# Run at most 10 at a time (submit script does this): ichorCNA downloads chromosome info from UCSC.
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/gse71378/full/p3_env.sh
IFS=$'\t' read -r S BAM D CFG < <(task_line "$S2_META/tasks_ichor.tsv" "$SLURM_ARRAY_TASK_ID")
W=${CFG%%kb_*}; MODE=${CFG#*kb_}; WIN=$((W*1000))
OD="$S2_OUT/ichorCNA/$S/$CFG"; mkdir -p "$OD"
[[ -f "$OD/.done" ]] && { echo "SKIP $S $CFG"; exit 0; }
r_modules
E="$ICHOR_EXT38"; gc="$E/gc_hg19_${W}kb.wig"; map="$E/map_hg19_${W}kb.wig"
cent="$E/GRCh37.p13_centromere_UCSC-gapTable.txt"
pon="$E/HD_ULP_PoN_$( [[ $W == 1000 ]] && echo 1Mb || echo 500kb )_median_normAutosome_mapScoreFiltered_median.rds"
for f in "$gc" "$map" "$cent" "$pon"; do [[ -s "$f" ]] || { echo "missing $f"; exit 4; }; done
wig="$OD/$S.wig"
[[ -s "$wig" ]] || "$READCOUNTER" --window "$WIN" --quality 20 --chromosome "$(printf '%s,' {1..22} X | sed 's/,$//')" "$BAM" > "$wig"
if [[ "$MODE" == "lowTF" ]]; then
  PA=(--ploidy 'c(2)' --normal 'c(0.95,0.99,0.995,0.999)' --maxCN 3 --estimateScPrevalence False --scStates 'c()')
else
  PA=(--ploidy 'c(2,3)' --normal 'c(0.5,0.6,0.7,0.8,0.9,0.95,0.99)' --maxCN 5 --estimateScPrevalence True --scStates 'c(1,3)')
fi
for try in 1 2 3; do            # retry: the UCSC download can fail when many jobs start together
  /usr/bin/time -v Rscript "$ICHOR_DIR/scripts/runIchorCNA.R" --id "$S" --WIG "$wig" \
    --gcWig "$gc" --mapWig "$map" --centromere "$cent" --normalPanel "$pon" "${PA[@]}" \
    --includeHOMD False --chrs 'c(1:22, "X")' --chrTrain 'c(1:22)' \
    --estimateNormal True --estimatePloidy True --txnE 0.9999 --txnStrength 10000 \
    --genomeBuild hg19 --genomeStyle NCBI --plotFileType png --outDir "$OD" && break
  echo "ichorCNA attempt $try failed - retrying in 60 s"; sleep 60
done
test -s "$OD/$S.params.txt"
touch "$OD/.done"; echo "ICHOR_DONE $S $CFG"
