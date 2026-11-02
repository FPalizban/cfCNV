#!/bin/bash
#SBATCH --job-name=bic_mapfix
#SBATCH --cpus-per-task=20
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
FIXED_WIG="$MAPROOT/${MER}.fixed.wig"
SIZES="$MAPROOT/${MER}.sizes.fixed"
BW="$MAPROOT/${MER}.bw"
BEDGRAPH="$MAPROOT/${MER}.bedGraph"

MAPDIR="$MAPROOT/$MER"
CHRFILE="$BICRES/chromosomes.txt"

mkdir -p "$MAPDIR"

echo "=========================================="
echo "BIC-seq2 mappability repair"
echo "=========================================="

[[ -s "$WIG" ]] || {
    echo "ERROR: WIG missing: $WIG"
    exit 2
}

[[ -s "${REF}.fai" ]] || {
    echo "ERROR: reference FAI missing: ${REF}.fai"
    exit 3
}

echo
echo "Original WIG first lines:"
head -5 "$WIG"

echo
echo "Creating corrected chromosome sizes..."

cut -f1,2 "${REF}.fai" \
| awk '$2 ~ /^[0-9]+$/ {print $1 "\t" $2}' \
> "$SIZES"

grep -q $'^1\t249250621$' "$SIZES" || {
    echo "ERROR: GRCh37 chr1 size check failed"
    exit 4
}

echo
echo "Sanitizing GEM 1.423 WIG directives..."

awk '
/^(variableStep|fixedStep)[[:space:]]/ {
    out=$1
    for(i=2;i<=NF;i++) {
        if(index($i,"=")>0)
            out=out " " $i
    }
    print out
    next
}
{
    print
}
' "$WIG" > "$FIXED_WIG"

[[ -s "$FIXED_WIG" ]] || {
    echo "ERROR: fixed WIG was not created"
    exit 5
}

echo
echo "Corrected WIG first lines:"
head -5 "$FIXED_WIG"

echo
echo "Checking corrected directives..."

awk '
/^(variableStep|fixedStep)[[:space:]]/ {
    for(i=2;i<=NF;i++) {
        if(index($i,"=")==0) {
            print "ERROR malformed directive:", $0 > "/dev/stderr"
            bad=1
            exit
        }
    }
}
END {
    if(bad)
        exit 1
}
' "$FIXED_WIG"

rm -f "$BW" "$BEDGRAPH"

echo
echo "Converting corrected WIG -> BigWig..."

"$WRAP" \
wigToBigWig \
"$FIXED_WIG" \
"$SIZES" \
"$BW"

[[ -s "$BW" ]] || {
    echo "ERROR: BigWig not created"
    exit 6
}

echo
echo "BigWig:"
ls -lh "$BW"

echo
echo "Converting BigWig -> BedGraph..."

"$WRAP" \
bigWigToBedGraph \
"$BW" \
"$BEDGRAPH"

[[ -s "$BEDGRAPH" ]] || {
    echo "ERROR: BedGraph not created"
    exit 7
}

echo
echo "BedGraph:"
ls -lh "$BEDGRAPH"

echo
echo "Creating exact chromosome 1-22 mappability files..."

rm -f "$MAPDIR"/"${MER}".*.txt

awk \
-v out="$MAPDIR" \
-v mer="$MER" '
BEGIN {
    OFS="\t"
}
$1 ~ /^([1-9]|1[0-9]|2[0-2])$/ && $4 == 1 {
    fn=out "/" mer "." $1 ".txt"
    print $2,$3 >> fn
}
' "$BEDGRAPH"

echo
echo "Validating 22 autosomes..."

N=0

while read -r CHR
do
    F="$MAPDIR/${MER}.${CHR}.txt"

    [[ -s "$F" ]] || {
        echo "ERROR: missing/empty mappability file: $F"
        exit 10
    }

    N=$((N+1))
done < "$CHRFILE"

[[ "$N" -eq 22 ]] || {
    echo "ERROR: expected 22 map files; found $N"
    exit 11
}

touch "$BICRES/RESOURCES_READY"

echo
echo "Mappability files = $N"
echo "RESOURCE MARKER:"
ls -lh "$BICRES/RESOURCES_READY"

# Fixed WIG is only an intermediate copy.
# Original GEM WIG remains available.
rm -f "$FIXED_WIG"

echo
echo "=========================================="
echo "BICSEQ2 RESOURCES READY"
echo "=========================================="
