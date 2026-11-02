# =====================================================================
# cfCNV Simulation 2 (HG008, GRCh38) - shared settings. Every script loads this file.
# Folder of these scripts on the cluster:  $P/00_metadata/simulation2/full
# =====================================================================
export P=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project
source "$P/00_metadata/config.sh"
source "$P/00_metadata/00_env.sh"                 # TOOLS, READCOUNTER, SPECTRE, MOSDEPTH, ICHOR_DIR ...
export S2="$P/00_metadata/simulation2/full"       # these scripts
export S2_META="$P/00_metadata/simulation2"       # manifests (mixtures_v2, tasks)
export S2_MIX="$P/05_simulation/simulation2_HG008/mixtures_v2"
export S2_OUT="$P/07_tool_outputs/simulation2"
export S2_MET="$P/09_metrics/simulation2"
export S2_LOG="$P/logs/simulation2"
export REF38="$P/02_reference/hg38_GIABv3/GRCh38_GIABv3.fa"
export REF38_AUTO="$P/02_reference/hg38_GIABv3/GRCh38_GIABv3.autosomes.fa"
export BLACKLIST38="$P/02_reference/hg38_GIABv3/hg38-blacklist.v2.bed"
export TRUTH38="$P/01_raw_data/NIST_GIAB_HG008/benchmark/GRCh38_HG008-T-V0.5_somatic-CNV_PASS.draftbenchmark.calls.bed"
export ICHOR_EXT38="$ICHOR_DIR/inst/extdata"
export R_LIBS_USER="$TOOLS/Rlib"
export PYD="$P/06_tools/envs/dudeml/bin/python"     # numpy / pandas / scikit-learn / pysam
export PYL="$P/06_tools/envs/lyceum/bin/python"     # LYCEUM (torch)
export LYC_REPO="$P/06_tools/ML_based/LYCEUM"
export BICRES38="$TOOLS/BICseq2_resources/GRCh38_GIABv3_150mer"
export BIC_READLEN=150                             # Element AVITI 2 x 150 bp (mean 147 bp)
export SIM1_DUDE_MODELS="$P/07_tool_outputs/simulation1/dudeML/models"
export SIM1_MEANN_V2="$P/07_tool_outputs/simulation1_rescue/CNV-MEANN_v2/features"
export BINS="10000 30000 50000 100000 500000 1000000"
export PATH="$TOOLS/bin:$PATH"
mkdir -p "$S2_OUT" "$S2_MET" "$S2_LOG"

task_line() { awk -F'\t' -v n=$(( $2 + 2 )) 'NR==n' "$1"; }        # row N (0-based) of a task file
pyd() { module load Python/3.10.8-GCCcore-12.2.0 >/dev/null 2>&1; "$PYD" "$@"; }   # dudeML python
r_modules() { module load R/4.2.2-foss-2022b; module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2; }
