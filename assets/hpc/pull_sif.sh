#!/bin/bash
# Download the ready-made cfCNV Singularity images (.sif) into one folder. Each image is a single file:
# nothing is unpacked and no temporary or cache files are written outside that folder.
# usage: bash pull_sif.sh <folder> [tag]
set -euo pipefail
D=${1:?usage: pull_sif.sh <folder> [tag]}; TAG=${2:-1.0.0}; REG=ghcr.io/fpalizban
mkdir -p "$D"
export SINGULARITY_TMPDIR="$D" APPTAINER_TMPDIR="$D" SINGULARITY_CACHEDIR="$D/.cache" APPTAINER_CACHEDIR="$D/.cache"
SING=$(command -v singularity || command -v apptainer)
for n in cfcnv-tools cfcnv-r cfcnv-cnvkit cfcnv-spectre cfcnv-fragle cfcnv-ml cfcnv-lyceum; do
  f="$D/$n-$TAG.sif"
  if [[ -s "$f" ]]; then echo "present  $f"; continue; fi
  "$SING" pull --disable-cache "$f" "oras://$REG/$n:$TAG-sif" && echo "OK       $f" || { echo "FAILED   $n"; rm -f "$f"; }
done
rm -rf "$D/.cache"
ls -lh "$D"/*.sif
