#!/bin/bash
#SBATCH --job-name=HG008_reference
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/reference_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/reference_%j.err

set -euo pipefail
module load SAMtools/1.16.1-GCC-11.3.0

P=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project
D="$P/02_reference/hg38_GIABv3"
mkdir -p "$D"

REF="$D/GRCh38_GIABv3.fa"
GZ="$D/GRCh38_GIABv3.fasta.gz"

URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz"

if [[ ! -s "$REF" ]]; then
    wget -c -O "$GZ" "$URL"
    gzip -t "$GZ"
    gzip -dc "$GZ" > "$REF.pending"
    mv "$REF.pending" "$REF"
fi

if [[ ! -s "$REF.fai" ]]; then
    samtools faidx "$REF"
fi

cut -f1,2 "$REF.fai" > "$D/chrom.sizes"

echo "REFERENCE_READY=$REF"
ls -lh "$REF" "$REF.fai"
