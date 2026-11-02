#!/bin/bash
#SBATCH --job-name=cfCNV_dl_pilot
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/download_pilot.%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/download_pilot.%j.err

set -euo pipefail

source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/config.sh
source $PROJECT/00_metadata/load_modules.sh

mkdir -p $RAW/shendure_cfDNA/sra $RAW/shendure_cfDNA/fastq
mkdir -p $RAW/PRJNA795275/sra $RAW/PRJNA795275/fastq

# Download 2 Shendure/Snyder runs
head -n 2 $PROJECT/00_metadata/shendure_sra_accessions.txt | while read SRR; do
  echo "Downloading Shendure/Snyder: $SRR"
  prefetch --max-size 200G $SRR -O $RAW/shendure_cfDNA/sra

  fasterq-dump $RAW/shendure_cfDNA/sra/$SRR \
    -O $RAW/shendure_cfDNA/fastq \
    --split-files \
    --threads $THREADS

  gzip -f $RAW/shendure_cfDNA/fastq/${SRR}*.fastq
done

# Download 2 PRJNA795275 runs
head -n 2 $PROJECT/00_metadata/PRJNA795275_sra_accessions.txt | while read SRR; do
  echo "Downloading PRJNA795275: $SRR"
  prefetch --max-size 200G $SRR -O $RAW/PRJNA795275/sra

  fasterq-dump $RAW/PRJNA795275/sra/$SRR \
    -O $RAW/PRJNA795275/fastq \
    --split-files \
    --threads $THREADS

  gzip -f $RAW/PRJNA795275/fastq/${SRR}*.fastq
done

echo "Pilot download finished."
