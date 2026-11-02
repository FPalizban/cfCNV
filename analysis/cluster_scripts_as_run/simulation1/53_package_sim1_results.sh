#!/bin/bash

set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project

source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

EXPORT_DIR="$PROJECT/13_exports"
STAGE="$EXPORT_DIR/cfCNV_Sim1_results_for_local"
ZIP="$EXPORT_DIR/cfCNV_Sim1_results_for_local.zip"

mkdir -p "$EXPORT_DIR"

# Remove old staging directory/archive
rm -rf "$STAGE"
rm -f "$ZIP"

mkdir -p \
    "$STAGE/09_metrics" \
    "$STAGE/08_harmonized_calls" \
    "$STAGE/truth" \
    "$STAGE/manifests" \
    "$STAGE/scripts"

echo "=========================================="
echo "Packaging cfCNV Simulation-1 results"
echo "=========================================="

# --------------------------------------------------
# 1. ALL METRICS / ANALYSIS TABLES
# --------------------------------------------------

if [[ -d "$PROJECT/09_metrics/simulation1" ]]; then
    cp -a \
        "$PROJECT/09_metrics/simulation1" \
        "$STAGE/09_metrics/"
else
    echo "WARNING: 09_metrics/simulation1 not found"
fi


# --------------------------------------------------
# 2. HARMONIZED CALLS, IF STORED SEPARATELY
# --------------------------------------------------

if [[ -d "$PROJECT/08_harmonized_calls/simulation1" ]]; then
    cp -a \
        "$PROJECT/08_harmonized_calls/simulation1" \
        "$STAGE/08_harmonized_calls/"
fi


# --------------------------------------------------
# 3. SIMULATION GROUND TRUTH
# --------------------------------------------------

TRUTH_DIR="$SIM/healthy_cfDNA_spikein/truth"

if [[ -d "$TRUTH_DIR" ]]; then
    cp -a \
        "$TRUTH_DIR"/. \
        "$STAGE/truth/"
else
    echo "WARNING: truth directory not found"
fi


# --------------------------------------------------
# 4. MANIFESTS
# --------------------------------------------------

for F in \
    master.tsv \
    qdnaseq_tasks.tsv \
    cnvkit_tasks.tsv \
    ichor_tasks.tsv \
    spectre_tasks.tsv \
    submitted_jobs.tsv
do
    if [[ -f "$SIM1_MAN/$F" ]]; then
        cp "$SIM1_MAN/$F" "$STAGE/manifests/"
    fi
done


# --------------------------------------------------
# 5. SCORING / HARMONIZATION SCRIPTS
# --------------------------------------------------

for F in \
    50_score_sim1.py \
    51_score_sim1.slurm \
    52_collect_resources.sh \
    20_collect_calls.py \
    21_score_sim1.py \
    22_bin_metrics.py
do
    if [[ -f "$PROJECT/00_metadata/$F" ]]; then
        cp "$PROJECT/00_metadata/$F" \
           "$STAGE/scripts/"
    fi
done


# --------------------------------------------------
# 6. CREATE SMALL README
# --------------------------------------------------

cat > "$STAGE/README.txt" <<TXT
cfCNV Simulation-1 local analysis package
=========================================

Generated:
$(date)

Project:
cfCNV Benchmarking

This archive contains the analysis-ready outputs needed for
downstream statistics, regression analyses, and visualization.

Important directories:

09_metrics/simulation1/
    Main benchmark tables and calculated metrics.

truth/
    Simulation-1 ground-truth CNVs.

manifests/
    Sample/tool/bin/depth task manifests.

scripts/
    Exact harmonization/scoring scripts used to generate metrics.

Large BAM files and raw sequencing files are intentionally excluded.

Large raw caller-output directories are also excluded because the
harmonized call and metric tables are sufficient for downstream
statistical analysis and figure generation.
TXT


# --------------------------------------------------
# 7. FILE INVENTORY
# --------------------------------------------------

find "$STAGE" -type f \
    | sed "s|$STAGE/||" \
    | sort \
    > "$STAGE/FILE_MANIFEST.txt"


# --------------------------------------------------
# 8. CHECKSUMS
# --------------------------------------------------

(
    cd "$STAGE"

    find . -type f \
        ! -name 'SHA256SUMS.txt' \
        -print0 \
        | sort -z \
        | xargs -0 sha256sum \
        > SHA256SUMS.txt
)


# --------------------------------------------------
# 9. CREATE ZIP
# --------------------------------------------------

cd "$EXPORT_DIR"

zip -rq \
    "$ZIP" \
    "$(basename "$STAGE")"


# --------------------------------------------------
# 10. FINAL REPORT
# --------------------------------------------------

echo
echo "=========================================="
echo "PACKAGE COMPLETE"
echo "=========================================="
echo
echo "Archive:"
echo "$ZIP"
echo

ls -lh "$ZIP"

echo
echo "Files included:"
find "$STAGE" -type f | wc -l

echo
echo "Main metrics:"
find "$STAGE/09_metrics" \
    -type f \
    -printf '%f\n' \
    2>/dev/null \
    | sort

echo
echo "=========================================="
echo "READY FOR SCP"
echo "=========================================="
