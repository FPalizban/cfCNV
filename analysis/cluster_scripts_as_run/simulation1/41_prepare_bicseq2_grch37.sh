#!/bin/bash
#SBATCH --job-name=bicseq2_map37
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/bicseq2_map_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/bicseq2_map_%j.err

set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project

source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

module load SAMtools/1.16.1-GCC-11.3.0
module load singularity/4.4.1

IMAGE="$TOOLS/containers/bicseq2_20210625.sif"

BICROOT="$REF37_DIR/BICseq2"
BICREF="$BICROOT/per_chrom"
BICMAP="$BICROOT/mappability"
CHRLIST="$BICROOT/chromosomes.txt"

mkdir -p \
"$BICROOT" \
"$BICREF" \
"$BICMAP" \
"$LOGS/simulation1"

[[ -s "$REF37_FASTA" ]] || {
    echo "ERROR: GRCh37 FASTA missing: $REF37_FASTA"
    exit 2
}

[[ -s "$IMAGE" ]] || {
    echo "ERROR: BIC-seq2 container missing: $IMAGE"
    exit 3
}

printf "%s\n" {1..22} > "$CHRLIST"

echo "Splitting GRCh37 into per-chromosome FASTAs..."

while read -r chr
do
    if [[ ! -s "$BICREF/${chr}.fa" ]]; then
        samtools faidx "$REF37_FASTA" "$chr" \
        > "$BICREF/${chr}.fa"
    fi
done < "$CHRLIST"


FIRST_BAM=$(awk -F'\t' 'NR==2 {print $2}' "$SIM1_MAN/master.tsv")

set +o pipefail
READ_LENGTH=$(
    samtools view "$FIRST_BAM" |
    awk 'NR==1 {print length($10); exit}'
)
set -o pipefail

[[ "$READ_LENGTH" =~ ^[0-9]+$ ]] || {
    echo "ERROR: could not determine read length"
    exit 4
}

echo "$READ_LENGTH" > "$BICROOT/read_length.txt"

echo "Detected read length: $READ_LENGTH"

MER="GRCh37.${READ_LENGTH}mer"

echo "Generating BIC-seq2 mappability files..."
echo "MER=$MER"

singularity exec \
-B "$PROJECT:$PROJECT" \
"$IMAGE" \
bash /BICSEQ2/src/prep_mappability.sh \
-l "$READ_LENGTH" \
-r GRCh37 \
-m 8 \
-n 32 \
"$REF37_FASTA" \
"$BICMAP" \
"$CHRLIST"


# Correct chromosome-specific splitting for numeric GRCh37 contigs
BEDGRAPH="$BICMAP/${MER}.bedGraph"

mkdir -p "$BICMAP/$MER"

if [[ -s "$BEDGRAPH" ]]; then

    while read -r chr
    do
        awk -v C="$chr" \
        '$1==C && $4==1 {print $2, $3}' \
        "$BEDGRAPH" \
        > "$BICMAP/$MER/${MER}.${chr}.txt"
    done < "$CHRLIST"

fi


echo "Checking outputs..."

missing=0

while read -r chr
do
    f="$BICMAP/$MER/${MER}.${chr}.txt"

    if [[ ! -s "$f" ]]; then
        echo "MISSING: $f"
        missing=$((missing+1))
    fi
done < "$CHRLIST"

if [[ "$missing" -ne 0 ]]; then
    echo "ERROR: $missing chromosome mappability files missing"
    exit 5
fi

cat > "$BICROOT/bicseq2_reference.env" <<ENV
export BICSEQ2_READ_LENGTH=$READ_LENGTH
export BICSEQ2_MER=$MER
export BICSEQ2_REFD=$BICREF
export BICSEQ2_MAPD=$BICMAP
export BICSEQ2_CHRLIST=$CHRLIST
ENV

echo "======================================"
echo "BIC-seq2 GRCh37 resources READY"
echo "READ_LENGTH=$READ_LENGTH"
echo "MER=$MER"
echo "REFD=$BICREF"
echo "MAPD=$BICMAP"
echo "======================================"
