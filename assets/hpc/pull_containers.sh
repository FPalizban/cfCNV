#!/bin/bash
# Download the six cfCNV images once into the Singularity cache that Nextflow uses (run on a node with internet).
#   export NXF_SINGULARITY_CACHEDIR=/path/with/space/singularity_cache
#   bash assets/hpc/pull_containers.sh 1.0.0
set -euo pipefail
TAG=${1:-1.0.0}; REG=${2:-ghcr.io/fpalizban}
D=${NXF_SINGULARITY_CACHEDIR:?set NXF_SINGULARITY_CACHEDIR first}; mkdir -p "$D"
SING=$(command -v singularity || command -v apptainer)
for n in cfcnv-tools cfcnv-r cfcnv-cnvkit cfcnv-spectre cfcnv-fragle cfcnv-bicseq2; do
  img="$D/$(echo "${REG}/${n}:${TAG}" | sed 's#[/:]#-#g').img"     # the file name Nextflow looks for
  if [[ -s "$img" ]]; then echo "present  $img"; continue; fi
  "$SING" pull "$img" "docker://${REG}/${n}:${TAG}" && echo "OK       $img"
done
ls -lh "$D"
