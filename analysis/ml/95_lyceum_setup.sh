#!/bin/bash
# =====================================================================
# cfCNV - one-time LYCEUM installation (run on a login node).
# Released environment file + repo; no changes to model code.
# =====================================================================
set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project
TOOLS="$PROJECT/06_tools/ML_based"
mkdir -p "$TOOLS"
cd "$TOOLS"

if [[ ! -d LYCEUM ]]; then
    git clone https://github.com/ciceklab/LYCEUM.git
fi
cd LYCEUM
git rev-parse HEAD > "$PROJECT/06_tools/LYCEUM_commit.txt"
echo "LYCEUM commit: $(cat "$PROJECT/06_tools/LYCEUM_commit.txt")"

# Conda env from the authors' file (Linux). Load your cluster's conda first,
# e.g.  module load Miniforge3   (adjust to the module name on your cluster)
if ! conda env list | grep -q '^lyceum_env '; then
    conda env create --name lyceum_env -f LYCEUM_environment.yml
fi

# sambamba is needed for read depth over the exon BED
conda run -n lyceum_env bash -c 'command -v sambamba || conda install -y -c bioconda sambamba'

# Record what the released target BED looks like (chromosome naming!)
echo
echo "Target BED (first lines):"
head -3 hg38_hglft_genome_64dc_dcbaa0_unique.bed
cut -f1 hg38_hglft_genome_64dc_dcbaa0_unique.bed | sort -u | head -30 | tr '\n' ' '
echo
echo
echo "Smoke test on the authors' toy data (CPU):"
conda run -n lyceum_env bash -c '
  set -e
  mkdir -p read_depths stats
  for f in ./example_data/*.bam; do
      samtools index "$f" 2>/dev/null || sambamba index "$f"
      sambamba depth base -L hg38_hglft_genome_64dc_dcbaa0_unique.bed "$f" > read_depths/$(basename "$f").txt
  done
  python ./scripts/preprocess_sample.py --readdepth ./read_depths --output ./processed_samples --target hg38_hglft_genome_64dc_dcbaa0_unique.bed
  python ./scripts/mean_std_calculator.py --processed_samples ./processed_samples --setname toySet --output ./stats/
  python ./scripts/LYCEUM_call.py --model lyceum --input ./processed_samples --output ./lyceum_calls_output --cnv exonlevel --batch_size 16 --stats_lookup ./stats/toySet_stats_lookup.npy
'
echo
echo "Toy output (first lines) - note the label values in the last column:"
head -5 lyceum_calls_output/* 2>/dev/null || ls lyceum_calls_output
