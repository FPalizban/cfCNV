# =====================================================================
# cfCNV real cancer cfDNA, GSE71378 / PRJNA291063 (Snyder et al. 2016), GRCh38 - shared settings.
# The authors' own GRCh37 BAMs (kircherlab mirror) are used directly - NO realignment - with the
# Simulation 1 GRCh37 resources (same BAM source as IH01-IH03); tool settings as in Phase 3.
# Scripts on the cluster:  $P/00_metadata/gse71378/full
# The tool scripts (A*, B*) are the Simulation 2 scripts pointed at GSE71378 folders, so the
# variable names S2_* are kept on purpose.
# =====================================================================
export P=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project
source "$P/00_metadata/config.sh"
source "$P/00_metadata/00_env.sh"                 # TOOLS, READCOUNTER, SPECTRE, MOSDEPTH, ICHOR_DIR ...
export S2="$P/00_metadata/gse71378/full"            # these scripts
export S2_META="$P/00_metadata/gse71378"            # inventory, manifests, task lists
export S2_OUT="$P/07_tool_outputs/gse71378"
export S2_MET="$P/09_metrics/gse71378"
export S2_LOG="$P/logs/gse71378"
export P3_FASTQ="$P/01_raw_data/GSE71378/fastq"
export P3_BAM="$P/04_aligned/GSE71378"
export REF38="$P/02_reference/hg38_GIABv3/GRCh38_GIABv3.fa"
export REF38_AUTO="$P/02_reference/hg38_GIABv3/GRCh38_GIABv3.autosomes.fa"
export BLACKLIST38="$P/02_reference/hg38_GIABv3/hg38-blacklist.v2.bed"
export ICHOR_EXT38="$ICHOR_DIR/inst/extdata"
export REF37="$P/02_reference/GRCh37/fasta/human_g1k_v37.fasta"
export MIRROR=https://kircherlab.bihealth.org/download/cfDNA
export TARGET_FULL_X=12                             # deep BAMs are down-sampled to ~12x (as Phase 3)
export R_LIBS_USER="$TOOLS/Rlib"
# read length of these data (written by 13_tasks after alignment; 150 until then)
export BIC_READLEN=101                              # GRCh37 101-mer mappability, as in Simulation 1
export BICRES38="$TOOLS/BICseq2_resources/GRCh37_101mer"   # (variable name kept; GRCh37 here)
# Building a new mappability track takes up to ~2 days. If none exists for this read length, the
# nearest finished one is used instead (read length is reported; note it in Methods).
if false; then
  _best=""; _bd=9999
  for _d in "$TOOLS"/BICseq2_resources/GRCh38_GIABv3_*mer; do
    [[ -f "$_d/RESOURCES_READY" ]] || continue
    _l=${_d##*_}; _l=${_l%mer}; _x=$(( _l > BIC_READLEN ? _l - BIC_READLEN : BIC_READLEN - _l ))
    (( _x < _bd )) && { _bd=$_x; _best=$_d; }
  done
  [[ -n "$_best" ]] && export BICRES38="$_best"
fi
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
