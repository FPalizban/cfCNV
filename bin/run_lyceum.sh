#!/usr/bin/env bash
# Run the released, pretrained LYCEUM model on one BAM (as in the cfCNV paper; model and weights unchanged).
# usage: run_lyceum.sh <bam> <GRCh37|GRCh38> <chain or NO_FILE> <threads> <out_segments.tsv>
set -euo pipefail
BAM=$1; BUILD=$2; CHAIN=$3; T=${4:-4}; OUT=$5
L=${LYCEUM_DIR:-/opt/LYCEUM}; TGT="$L/hg38_hglft_genome_64dc_dcbaa0_unique.bed"
mkdir -p read_depths processed stats calls
# sambamba requires the header to declare SO:coordinate. An indexed BAM is coordinate-sorted (samtools refuses to
# index anything else), so when only the header line is missing, a local copy with a corrected header is used.
samtools view -H "$BAM" > header.orig.sam
if ! grep -q $'^@HD.*\tSO:coordinate' header.orig.sam; then
  awk 'BEGIN{OFS="\t"} NR==1 && !/^@HD/ {print "@HD\tVN:1.6\tSO:coordinate"} /^@HD/ {gsub(/\tSO:[^\t]*/,""); print $0"\tSO:coordinate"; next} {print}' header.orig.sam > header.sam
  samtools reheader header.sam "$BAM" > input.sorted_header.bam
  samtools index input.sorted_header.bam
  BAM=input.sorted_header.bam
  echo "BAM header did not declare SO:coordinate; using a copy with a corrected header" >&2
fi
samtools idxstats "$BAM" | cut -f1 > bam_chroms.txt
if [[ "$BUILD" == "GRCh37" ]]; then
  [[ -s "$CHAIN" ]] || { echo "GRCh37 BAMs need the hg38ToHg19 chain (--liftover_chain)"; exit 2; }
  lyceum_tools.py targets --hg38-bed "$TGT" --bam-chroms bam_chroms.txt --out-bed targets.bam.bed --out-map map.tsv --chain "$CHAIN"
else
  lyceum_tools.py targets --hg38-bed "$TGT" --bam-chroms bam_chroms.txt --out-bed targets.bam.bed --out-map map.tsv
fi
sambamba depth base -t "$T" -L targets.bam.bed "$BAM" > depth.bam.txt
lyceum_tools.py depth --in depth.bam.txt --map map.tsv --out read_depths/LYC.bam.txt --out-bed targets.hg38.bed
rm -f depth.bam.txt
python "$L/scripts/preprocess_sample.py" --readdepth read_depths --output processed --target targets.hg38.bed
python "$L/scripts/mean_std_calculator.py" --processed_samples processed --setname cfcnv --output stats
python "$L/scripts/LYCEUM_call.py" --model lyceum --input processed --output calls --cnv exonlevel \
    --batch_size 16 --stats_lookup stats/cfcnv_stats_lookup.npy
lyceum_tools.py calls --in calls/LYC_exonlevel.csv --map map.tsv --out "$OUT"
cp calls/LYC_exonlevel.csv "${OUT%.tsv}.exonlevel.csv"
rm -rf read_depths processed input.sorted_header.bam input.sorted_header.bam.bai header.sam header.orig.sam
