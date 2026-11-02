#!/bin/bash
#SBATCH --job-name=cfcnv_ref37
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/ref37_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/ref37_%j.err

set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project

source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

module load SAMtools/1.16.1-GCC-11.3.0

mkdir -p "$REF37_DIR"
mkdir -p "$PROJECT/logs"

FASTA="${REF37_FASTA:-$REF37_DIR/human_g1k_v37.fasta}"
GZ="${FASTA}.gz"

URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"

echo "=========================================="
echo "cfCNV GRCh37 REFERENCE PREPARATION"
echo "FASTA = $FASTA"
echo "=========================================="

# --------------------------------------------------
# Reuse an already downloaded FASTA if one exists
# --------------------------------------------------

if [[ ! -s "$FASTA" ]]; then

    ALT=$(find "$REF37_DIR" \
        -maxdepth 2 \
        -type f \
        -name 'human_g1k_v37.fasta' \
        -size +2G \
        2>/dev/null | head -1 || true)

    if [[ -n "$ALT" && "$ALT" != "$FASTA" ]]; then
        echo "Reusing existing FASTA:"
        echo "$ALT"

        mkdir -p "$(dirname "$FASTA")"
        ln -sf "$ALT" "$FASTA"
    fi
fi


# --------------------------------------------------
# Download compressed FASTA only when necessary
# --------------------------------------------------

if [[ ! -s "$FASTA" ]]; then

    if [[ ! -s "$GZ" ]]; then

        ALT_GZ=$(find "$REF37_DIR" \
            -maxdepth 2 \
            -type f \
            -name 'human_g1k_v37.fasta.gz' \
            -size +800M \
            2>/dev/null | head -1 || true)

        if [[ -n "$ALT_GZ" && "$ALT_GZ" != "$GZ" ]]; then

            echo "Reusing downloaded gzip file:"
            echo "$ALT_GZ"

            mkdir -p "$(dirname "$GZ")"
            ln -sf "$ALT_GZ" "$GZ"

        else

            echo "Downloading GRCh37..."

            mkdir -p "$(dirname "$GZ")"

            wget \
                --tries=5 \
                --timeout=60 \
                --continue \
                -O "$GZ" \
                "$URL"
        fi
    fi


    # --------------------------------------------------
    # Decompress.
    #
    # The 1000 Genomes file can produce:
    # "decompression OK, trailing garbage ignored"
    # and gzip return code 2 despite producing a valid
    # FASTA. Therefore validate the output explicitly.
    # --------------------------------------------------

    echo "Decompressing GRCh37..."

    TMP="${FASTA}.tmp.$$"

    rm -f "$TMP"

    set +e
    gzip -dc "$GZ" > "$TMP"
    RC=$?
    set -e

    if [[ "$RC" -ne 0 && "$RC" -ne 2 ]]; then
        echo "ERROR: gzip returned $RC"
        rm -f "$TMP"
        exit 2
    fi

    if [[ ! -s "$TMP" ]]; then
        echo "ERROR: decompressed FASTA is empty"
        exit 3
    fi

    if ! head -1 "$TMP" | grep -q '^>'; then
        echo "ERROR: decompressed file does not look like FASTA"
        rm -f "$TMP"
        exit 4
    fi

    mv "$TMP" "$FASTA"
fi


# --------------------------------------------------
# Basic FASTA validation
# --------------------------------------------------

echo
echo "Reference file:"
ls -lh "$FASTA"

[[ -s "$FASTA" ]] || {
    echo "ERROR: FASTA missing"
    exit 5
}


# --------------------------------------------------
# FASTA index
# --------------------------------------------------

if [[ ! -s "${FASTA}.fai" ]]; then
    echo "Building FASTA index..."
    samtools faidx "$FASTA"
fi


# --------------------------------------------------
# Validate GRCh37 chr1 length
# --------------------------------------------------

CHR1LEN=$(awk '$1=="1"{print $2}' "${FASTA}.fai")

echo "chr1 length = $CHR1LEN"

if [[ "$CHR1LEN" != "249250621" ]]; then
    echo "ERROR: reference is not expected human_g1k_v37 GRCh37"
    exit 6
fi


# --------------------------------------------------
# Check BAM/reference chromosome naming
# --------------------------------------------------

FIRST_BAM=$(awk -F'\t' 'NR==2 {print $2}' \
"$SIM1_MAN/master.tsv")

[[ -s "$FIRST_BAM" ]] || {
    echo "ERROR: example BAM missing"
    echo "$FIRST_BAM"
    exit 7
}

BAM_CONTIG=$(
    samtools view -H "$FIRST_BAM" |
    awk -F'\t' '
    /^@SQ/ {
        for(i=1;i<=NF;i++) {
            if($i ~ /^SN:/) {
                sub(/^SN:/,"",$i)
                print $i
                exit
            }
        }
    }'
)

REF_CONTIG=$(head -1 "${FASTA}.fai" | cut -f1)

echo "BAM first contig = $BAM_CONTIG"
echo "REF first contig = $REF_CONTIG"

if [[ "$BAM_CONTIG" != "$REF_CONTIG" ]]; then
    echo "ERROR: BAM/reference chromosome naming mismatch"
    exit 8
fi


# --------------------------------------------------
# Autosomal BED
# --------------------------------------------------

awk '
$1 ~ /^([1-9]|1[0-9]|2[0-2])$/ {
    print $1 "\t0\t" $2
}
' "${FASTA}.fai" > "$REF37_DIR/autosomes.bed"


echo
echo "=========================================="
echo "GRCh37 REFERENCE READY"
echo "=========================================="
echo "FASTA = $FASTA"
echo "INDEX = ${FASTA}.fai"
echo "=========================================="
