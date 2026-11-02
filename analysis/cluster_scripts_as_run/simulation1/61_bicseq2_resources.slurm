#!/bin/bash
#SBATCH --job-name=bic_map101
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
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

REF="$PROJECT/02_reference/GRCh37/fasta/human_g1k_v37.fasta"

REFSPLIT="$BICRES/reference"
MAPROOT="$BICRES/mappability"

REFBASE="human_g1k_v37"
READ_LENGTH=101
MER="${REFBASE}.${READ_LENGTH}mer"
MAPDIR="$MAPROOT/$MER"

mkdir -p \
"$REFSPLIT" \
"$MAPROOT" \
"$MAPDIR"

seq 1 22 > "$CHRFILE"

echo "============================================"
echo "BIC-seq2 GRCh37 101-mer resource generation"
echo "============================================"
echo "Reference   : $REF"
echo "Read length : $READ_LENGTH"
echo

[[ -s "$REF" ]] || {
    echo "ERROR: reference not found"
    exit 2
}

if [[ ! -s "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi


# -------------------------------------------------
# SPLIT GRCh37 INTO AUTOSOMAL FASTA FILES
# -------------------------------------------------

echo "Preparing chromosome FASTAs..."

while read -r CHR
do

    OUT="$REFSPLIT/${CHR}.fa"

    if [[ ! -s "$OUT" ]]; then
        samtools faidx "$REF" "$CHR" > "$OUT"
    fi

done < "$CHRFILE"


# -------------------------------------------------
# CHECK IF MAPPABILITY ALREADY EXISTS
# -------------------------------------------------

READY=1

while read -r CHR
do
    F="$MAPDIR/${MER}.${CHR}.txt"

    if [[ ! -s "$F" ]]; then
        READY=0
        break
    fi

done < "$CHRFILE"


if [[ "$READY" -eq 0 ]]; then

    echo
    echo "Checking programs inside BIC-seq2 container..."

    "$WRAP" bash -lc '
        command -v gem-indexer
        command -v gem-mappability
        command -v gem-2-wig
        command -v wigToBigWig
        command -v bigWigToBedGraph
    '

    echo
    echo "Running official BIC-seq2 mappability workflow..."

    "$WRAP" \
    bash /BICSEQ2/src/prep_mappability.sh \
        -l 101 \
        -m 8 \
        -n 24 \
        -r "$REFBASE" \
        "$REF" \
        "$MAPROOT" \
        "$CHRFILE"


    # -------------------------------------------------
    # IMPORTANT:
    # official older script uses grep CHR; for chromosome
    # names 1,2,... this can match 10/11/etc.
    # Re-split the generated BedGraph EXACTLY by column 1.
    # -------------------------------------------------

    echo
    echo "Creating exact autosomal mappability files..."

    BEDGRAPH="$MAPROOT/${MER}.bedGraph"

    [[ -s "$BEDGRAPH" ]] || {
        echo "ERROR: BedGraph not produced"
        exit 10
    }

    rm -f "$MAPDIR"/"${MER}".*.txt

    awk \
    -v out="$MAPDIR" \
    -v mer="$MER" \
    '
    BEGIN {
        OFS="\t"
        for(i=1;i<=22;i++)
            keep[i]=1
    }

    ($1 in keep) && $4==1 {
        fn=out "/" mer "." $1 ".txt"
        print $2,$3 > fn
    }
    ' "$BEDGRAPH"

else

    echo "Existing mappability resources detected; skipping generation."

fi


# -------------------------------------------------
# VALIDATE
# -------------------------------------------------

echo
echo "Validating all chromosomes..."

NREF=0
NMAP=0

while read -r CHR
do

    RF="$REFSPLIT/${CHR}.fa"
    MF="$MAPDIR/${MER}.${CHR}.txt"

    [[ -s "$RF" ]] || {
        echo "ERROR missing reference: $RF"
        exit 20
    }

    [[ -s "$MF" ]] || {
        echo "ERROR missing mappability: $MF"
        exit 21
    }

    NREF=$((NREF+1))
    NMAP=$((NMAP+1))

done < "$CHRFILE"

echo "Reference chromosomes   = $NREF"
echo "Mappability chromosomes = $NMAP"

[[ "$NREF" -eq 22 ]]
[[ "$NMAP" -eq 22 ]]

touch "$BICRES/RESOURCES_READY"

echo
echo "============================================"
echo "BICSEQ2 RESOURCES READY"
echo "============================================"
