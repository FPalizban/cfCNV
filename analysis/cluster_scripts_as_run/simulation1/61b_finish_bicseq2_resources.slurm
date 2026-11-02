#!/bin/bash
#SBATCH --job-name=bic_mapfix
#SBATCH --cpus-per-task=25
#SBATCH --mem=120G
#SBATCH --time=08:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/bic_mapfix_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/bic_mapfix_%j.err

set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project

source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

TOOLS="$PROJECT/06_tools"
WRAP="$TOOLS/bin/run_bicseq2"

BICRES="$TOOLS/BICseq2_resources/GRCh37_101mer"
MAPROOT="$BICRES/mappability"

REF="$PROJECT/02_reference/GRCh37/fasta/human_g1k_v37.fasta"

MER="human_g1k_v37.101mer"

WIG="$MAPROOT/${MER}.wig"
SIZES="$MAPROOT/${MER}.sizes.fixed"
BW="$MAPROOT/${MER}.bw"
BEDGRAPH="$MAPROOT/${MER}.bedGraph"

MAPDIR="$MAPROOT/$MER"
CHRFILE="$BICRES/chromosomes.txt"

mkdir -p "$MAPDIR"

echo "=========================================="
echo "Resume BIC-seq2 mappability preparation"
echo "=========================================="

[[ -s "$WIG" ]] || {
    echo "ERROR: existing WIG missing: $WIG"
    exit 2
}

[[ -s "${REF}.fai" ]] || {
    echo "ERROR: FASTA index missing: ${REF}.fai"
    exit 3
}

echo "Using existing WIG:"
ls -lh "$WIG"

echo
echo "Creating corrected chromosome sizes..."

cut -f1,2 "${REF}.fai" \
| awk '$2 ~ /^[0-9]+$/ {print $1 "\t" $2}' \
> "$SIZES"

echo "First chromosome sizes:"
head "$SIZES"

grep -q $'^1\t249250621$' "$SIZES" || {
    echo "ERROR: corrected sizes file does not contain expected GRCh37 chr1"
    exit 4
}

echo
echo "Converting WIG -> BigWig..."

"$WRAP" \
wigToBigWig \
"$WIG" \
"$SIZES" \
"$BW"

[[ -s "$BW" ]] || {
    echo "ERROR: BigWig not created"
    exit 5
}

echo
echo "Converting BigWig -> BedGraph..."

"$WRAP" \
bigWigToBedGraph \
"$BW" \
"$BEDGRAPH"

[[ -s "$BEDGRAPH" ]] || {
    echo "ERROR: BedGraph not created"
    exit 6
}

echo
echo "Creating exact autosomal mappability files..."

rm -f "$MAPDIR"/"${MER}".*.txt

awk \
-v out="$MAPDIR" \
-v mer="$MER" \
'
BEGIN {
    OFS="\t"
}
$1 ~ /^([1-9]|1[0-9]|2[0-2])$/ && $4==1 {
    fn=out "/" mer "." $1 ".txt"
    print $2,$3 >> fn
}
' "$BEDGRAPH"

echo
echo "Validating chromosomes..."

N=0

while read -r CHR
do
    F="$MAPDIR/${MER}.${CHR}.txt"

    [[ -s "$F" ]] || {
        echo "ERROR: missing mappability file $F"
        exit 10
    }

    N=$((N+1))
done < "$CHRFILE"

[[ "$N" -eq 22 ]] || {
    echo "ERROR: expected 22 chromosomes; found $N"
    exit 11
}

touch "$BICRES/RESOURCES_READY"

echo
echo "=========================================="
echo "BICSEQ2 RESOURCES READY"
echo "Chromosomes = $N"
echo "=========================================="
