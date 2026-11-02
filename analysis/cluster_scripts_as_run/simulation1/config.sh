export PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project

export RAW=$PROJECT/01_raw_data
export REF_DIR=$PROJECT/02_reference
export QC=$PROJECT/03_qc
export BAM_DIR=$PROJECT/04_processed_bam
export SIM=$PROJECT/05_simulation
export TOOLS=$PROJECT/06_tools
export OUT=$PROJECT/07_tool_outputs
export HARM=$PROJECT/08_harmonized_calls
export METRICS=$PROJECT/09_metrics
export BIO=$PROJECT/10_biological_validation
export FIG=$PROJECT/11_figures
export TABLES=$PROJECT/12_tables
export LOGS=$PROJECT/logs

export THREADS=8

export BINSIZES="100000 500000 1000000"
export DEPTHS="0.5 1 5 10"
export TUMOR_FRACTIONS="0.01 0.05 0.10 0.20 0.30"

# Reference paths will be filled later
export REF=$REF_DIR/hg38/hg38.fa
export CHROMSIZES=$REF_DIR/hg38/hg38.chrom.sizes

# ============================================================
# cfCNV benchmarking design
# ============================================================

# Shendure/Snyder BAMs are aligned to GRCh37
export REF37_DIR=$PROJECT/02_reference/GRCh37

# NIST HG008 benchmark is GRCh38
export REF38_DIR=$PROJECT/02_reference/GRCh38

# Benchmark bin-size grid
export BIN_SIZES="10000 25000 50000 100000 250000 500000 1000000"

# Controlled sequencing-depth grid
export DEPTH_LEVELS="0.1 0.25 0.5 1 2 5 10"

# NIST tumor-fraction grid
export TUMOR_FRACTIONS="0 0.01 0.02 0.05 0.10 0.20 0.30"

# Simulation replicates
export N_SIM_REPS=5

export SHENDURE_DIR=$RAW/shendure_cfDNA/original_BAM
export PATH=$TOOLS/bin:$PATH
