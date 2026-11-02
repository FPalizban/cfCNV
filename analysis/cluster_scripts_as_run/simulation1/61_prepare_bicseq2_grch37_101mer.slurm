#!/bin/bash
#SBATCH --job-name=bic_map101
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/bic_map101_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/bic_map101_%j.err

set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project
source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

module load SAMtools/1.16.1-GCC-11.3.0

TOOLS="$PROJECT/06_tools"
WRAP="$TOOLS/bin/run_bicseq2"

BICRES="$TOOLS/BICseq2_resources/GRCh37_101mer"
CHRFILE="$BICRES/chromosomes.txt"
REFDIR="$BICRES/ref"
MAPROOT="$BICRES/mappability"

REF="${REF37_FASTA:-$PROJECT/02_reference/GRCh37/fasta/human_g1k_v37.fasta}"

REFBASE="human_g1k_v37"
READ_LENGTH=101
MER="${REFBASE}.${READ_LENGTH}mer"
MAPDIR="$MAPROOT/$MER"

mkdir -p "$BICRES" "$REFDIR" "$MAPROOT" "$MAPDIR"

seq 1 22 > "$CHRFILE"

echo "========================================"
echo "BIC-seq2 GRCh37 101-mer resources"
echo "========================================"
echo "Reference:   $REF"
echo "Read length: $READ_LENGTH"
echo "Chromosomes: 1-22"
echo

[[ -s "$REF" ]] || {
    echo "ERROR: reference missing: $REF"
    exit 2
}

if [[ ! -s "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi

# Ensure expected GRCh37 autosomes exist.
while read -r chr
do
    grep -q -P "^${chr}\t" "${REF}.fai" || {
        echo "ERROR: chromosome $chr absent from reference"
        exit 3
    }
done < "$CHRFILE"


# --------------------------------------------------
# Per-chromosome FASTA files required by NBICseq-norm
# --------------------------------------------------

echo "Preparing per-chromosome reference FASTAs..."

while read -r chr
do
    OUTFA="$REFDIR/${REFBASE}.${chr}.fa"

    if [[ ! -s "$OUTFA" ]]; then
        samtools faidx "$REF" "$chr" > "$OUTFA"
    fi

done < "$CHRFILE"


# --------------------------------------------------
# Determine whether mappability is already complete
# --------------------------------------------------

COMPLETE=1

while read -r chr
do
    MAP="$MAPDIR/${MER}.${chr}.txt"

    if [[ ! -s "$MAP" ]]; then
        COMPLETE=0
    fi

done < "$CHRFILE"


if [[ "$COMPLETE" -eq 1 ]]; then

    echo
    echo "All 22 mappability files already exist."
    echo "Skipping GEM regeneration."

else

    echo
    echo "Generating whole-genome 101-mer mappability..."
    echo "This is the long step."

    "$WRAP" bash -s -- \
        "$REF" \
        "$MAPROOT" \
        "$REFBASE" \
        "$MER" \
        "${REF}.fai" <<'INNER'

set -euo pipefail

REF="$1"
MAPROOT="$2"
REFBASE="$3"
MER="$4"
FAI="$5"

cd "$MAPROOT"

echo "Checking required programs..."

for X in \
    gem-indexer \
    gem-mappability \
    gem-2-wig \
    wigToBigWig \
    bigWigToBedGraph
do
    command -v "$X" || {
        echo "ERROR: $X not available in BIC-seq2 container"
        exit 10
    }
done

echo
echo "===== GEM INDEX ====="

gem-indexer \
    -i "$REF" \
    -o "$MAPROOT/$REFBASE" \
    -T 8


echo
echo "===== GEM MAPPABILITY ====="

gem-mappability \
    -m 2 \
    -I "$MAPROOT/${REFBASE}.gem" \
    -l 101 \
    -o "$MAPROOT/$MER" \
    -T 24 \
    2> "$MAPROOT/${MER}.mappability.log"


echo
echo "===== GEM TO WIG ====="

gem-2-wig \
    -I "$MAPROOT/${REFBASE}.gem" \
    -i "$MAPROOT/${MER}.mappability" \
    -o "$MER"


# Use the actual FASTA index as chromosome sizes.
cut -f1,2 "$FAI" > "$MAPROOT/${MER}.sizes.fai"


echo
echo "===== WIG TO BIGWIG ====="

wigToBigWig \
    "$MAPROOT/${MER}.wig" \
    "$MAPROOT/${MER}.sizes.fai" \
    "$MAPROOT/${MER}.bw"


echo
echo "===== BIGWIG TO BEDGRAPH ====="

bigWigToBedGraph \
    "$MAPROOT/${MER}.bw" \
    "$MAPROOT/${MER}.bedGraph"

INNER


    # --------------------------------------------------
    # Exact chromosome splitting.
    #
    # Do NOT use a simple grep for chromosome "1",
    # because that can also match 10,11,... or positions.
    # --------------------------------------------------

    mkdir -p "$MAPDIR"

    rm -f "$MAPDIR"/"${MER}".*.txt

    awk \
        -v out="$MAPDIR" \
        -v mer="$MER" \
        '
        BEGIN {
            OFS="\t"
            for (i=1; i<=22; i++)
                keep[i]=1
        }

        ($1 in keep) && $4==1 {
            fn=out "/" mer "." $1 ".txt"
            print $2,$3 >> fn
        }
        ' \
        "$MAPROOT/${MER}.bedGraph"

fi


# --------------------------------------------------
# Final validation
# --------------------------------------------------

echo
echo "===== VALIDATING RESOURCES ====="

NFA=0
NMAP=0

while read -r chr
do
    FA="$REFDIR/${REFBASE}.${chr}.fa"
    MAP="$MAPDIR/${MER}.${chr}.txt"

    [[ -s "$FA" ]] || {
        echo "ERROR: missing FASTA $FA"
        exit 20
    }

    [[ -s "$MAP" ]] || {
        echo "ERROR: missing mappability $MAP"
        exit 21
    }

    NFA=$((NFA+1))
    NMAP=$((NMAP+1))

done < "$CHRFILE"

echo "Per-chromosome FASTAs : $NFA"
echo "Mappability files     : $NMAP"

[[ "$NFA" -eq 22 ]]
[[ "$NMAP" -eq 22 ]]

touch "$BICRES/RESOURCES_READY"

echo
echo "========================================"
echo "BICSEQ2 RESOURCES READY"
echo "========================================"
