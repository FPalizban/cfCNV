#!/bin/bash

source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/config.sh
module load SAMtools/1.16.1-GCC-11.3.0

MAN=$SIM/healthy_cfDNA_spikein/manifests/analysis_manifest.tsv

TOTAL=0
READY=0
BAD=0

while IFS=$'\t' read -r \
ID BACKGROUND CLASS REP DEPTH BAM TRUTH
do

    [[ "$ID" == "sample_id" ]] && continue

    TOTAL=$((TOTAL+1))

    if [[ -s "$BAM" ]] && samtools quickcheck "$BAM" 2>/dev/null
    then
        READY=$((READY+1))
    else
        echo "MISSING/FAILED: $ID"
        echo "  $BAM"
        BAD=$((BAD+1))
    fi

done < "$MAN"

echo
echo "========================================"
echo "cfCNV SIMULATION 1 STATUS"
echo "========================================"
echo "Expected: $TOTAL"
echo "Ready:    $READY"
echo "Missing:  $BAD"
echo "========================================"
