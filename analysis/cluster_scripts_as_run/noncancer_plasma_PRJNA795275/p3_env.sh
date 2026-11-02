# =====================================================================
# cfCNV Phase 3 (real cfDNA, PRJNA795275, GRCh38) - shared settings. Every script loads this file.
# Scripts on the cluster:  $P/00_metadata/phase3/full
# The tool scripts (A*, B*) are the Simulation 2 scripts pointed at Phase 3 folders, so the
# variable names S2_* are kept on purpose.
# =====================================================================
export P=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project
source "$P/00_metadata/config.sh"
source "$P/00_metadata/00_env.sh"                 # TOOLS, READCOUNTER, SPECTRE, MOSDEPTH, ICHOR_DIR ...
export S2="$P/00_metadata/phase3/full"            # these scripts
export S2_META="$P/00_metadata/phase3"            # inventory, manifests, task lists
export S2_OUT="$P/07_tool_outputs/phase3"
export S2_MET="$P/09_metrics/phase3"
export S2_LOG="$P/logs/phase3"
export P3_FASTQ="$P/01_raw_data/PRJNA795275/fastq"
export P3_BAM="$P/04_aligned/PRJNA795275"
export REF38="$P/02_reference/hg38_GIABv3/GRCh38_GIABv3.fa"
export REF38_AUTO="$P/02_reference/hg38_GIABv3/GRCh38_GIABv3.autosomes.fa"
export BLACKLIST38="$P/02_reference/hg38_GIABv3/hg38-blacklist.v2.bed"
export ICHOR_EXT38="$ICHOR_DIR/inst/extdata"
export R_LIBS_USER="$TOOLS/Rlib"
# read length of these data (written by 13_tasks after alignment; 150 until then)
export BIC_READLEN=$(cat "$S2_META/p3_readlen.txt" 2>/dev/null || echo 150)
export BICRES38="$TOOLS/BICseq2_resources/GRCh38_GIABv3_${BIC_READLEN}mer"
export FRAGLE_ENV="$TOOLS/Fragle_env"               # conda-pack environment (Zenodo 13968359)
export FRAGLE_SRC="$TOOLS/FRAGLE-main"              # code (github.com/skandlab/FRAGLE)
export TARGET_MAX_X=12                              # deep samples are subsampled to ~12x before alignment
export PATH="$TOOLS/bin:$PATH"
mkdir -p "$S2_OUT" "$S2_MET" "$S2_LOG" "$P3_FASTQ" "$P3_BAM"

task_line() { awk -F'\t' -v n=$(( $2 + 2 )) 'NR==n' "$1"; }        # row N (0-based) of a task file
r_modules() { module load R/4.2.2-foss-2022b; module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2; }
aln_modules() { source "$S2_META/p3_modules.sh"; }                  # BWA + SAMtools, chosen by setup
SYSPY=$(command -v python3)                                          # system python (before modules)
export SYSPY
