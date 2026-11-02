#!/usr/bin/env bash
# Run the released, pretrained LYCEUM model on one BAM (as in the cfCNV paper; model and weights unchanged).
# usage: run_lyceum.sh <bam> <GRCh37|GRCh38> <chain or NO_FILE> <threads> <out_segments.tsv>
set -euo pipefail
BAM=$1; BUILD=$2; CHAIN=$3; T=${4:-4}; OUT=$5
L=${LYCEUM_DIR:-/opt/LYCEUM}; TGT="$L/hg38_hglft_genome_64dc_dcbaa0_unique.bed"
mkdir -p read_depths processed stats calls
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
rm -rf read_depths processed
