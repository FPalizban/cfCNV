# =====================================================================
# cfCNV ML dataset + baseline model - shared settings. Every script loads this file.
# Scripts on the cluster:  $P/00_metadata/ml_dataset
# Everything this package WRITES goes to:  $P/13_ml_dataset  and  $P/logs/ml_dataset
# Everything else (BAMs, truth files, references, manifests) is only READ, never changed.
# =====================================================================
export P=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project
export ML="$P/00_metadata/ml_dataset"           # these scripts
export MLO="$P/13_ml_dataset"                   # all outputs of this package
export MLLOG="$P/logs/ml_dataset"               # Slurm logs of this package

# python with numpy / pandas / scikit-learn / pysam (the one used for dudeML and Simulation 2)
export PYD="$P/06_tools/envs/dudeml/bin/python"

# references (read only)
export REF37="$P/02_reference/GRCh37/fasta/human_g1k_v37.fasta"          # Sim1 + GSE71378
export REF38="$P/02_reference/hg38_GIABv3/GRCh38_GIABv3.fa"              # HG008 + PRJNA795275
export BLACKLIST38="$P/02_reference/hg38_GIABv3/hg38-blacklist.v2.bed"

# inputs of the four arms (read only)
export SIM1_MANIFEST="$P/05_simulation/healthy_cfDNA_spikein/manifests/analysis_manifest.tsv"
export SIM1_TRUTH_DIR="$P/05_simulation/healthy_cfDNA_spikein/truth"
export SIM2_MANIFEST="$P/00_metadata/simulation2/HG008_mixtures_realised.tsv"
export TRUTH38="$P/01_raw_data/NIST_GIAB_HG008/benchmark/GRCh38_HG008-T-V0.5_somatic-CNV_PASS.draftbenchmark.calls.bed"
export PRJNA_MANIFEST="$P/00_metadata/phase3/p3_bams.tsv"
export GSE_MANIFEST="$P/00_metadata/gse71378/p3_bams.tsv"
export GSE_LABELS="$ML/data/cancer_consensus_window_states.tsv"         # consensus of the callers (from the repo)
export GSE_TF="$ML/data/cancer_tumour_fraction.tsv"                     # ichorCNA estimates (from the repo)

mkdir -p "$MLO" "$MLLOG"

task_line() { awk -F'\t' -v n=$(( $2 + 2 )) 'NR==n' "$1"; }             # row N (0-based) of a table with header
pyd() { module load Python/3.10.8-GCCcore-12.2.0 >/dev/null 2>&1; "$PYD" "$@"; }
