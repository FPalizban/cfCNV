#!/bin/bash
#SBATCH --job-name=Spectre-cnv_analysis
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=cnv_analysis.out
#SBATCH --error=cnv_analysis.err
#SBATCH --array=1-10%5  # Adjust the range and step size based on the number of BAM files

# Load necessary modules
module load python/3.9
#module load mosdepth/0.3.4

# Define paths and parameters
BAM_DIR="/work/FAC/FBM/DBC/cdessim2/default/fpalizba/cfDNA-CNV/validation/simulated-samples"
REFERENCE_GENOME="/work/FAC/FBM/DBC/cdessim2/default/fpalizba/cfDNA-CNV/tools/Spectre/hs37d5.fa"
OUTPUT_BASE_DIR="/work/FAC/FBM/DBC/cdessim2/default/fpalizba/cfDNA-CNV/validation/simulated-samples/run-results/Spectre-1000000-IH03"

# Get the BAM file corresponding to the array task index
bam_file=$(ls -1 ${BAM_DIR}/*.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)
sample_name=$(basename "${bam_file%.*}")

# Step 1: Run mosdepth
#mkdir -p "${OUTPUT_BASE_DIR}/${sample_name}_mosdepth"
./mosdepth -n --fast-mode --by 1000000 "${OUTPUT_BASE_DIR}/${sample_name}_mosdepth/${sample_name}" "${bam_file}"

# Step 2: Move mosdepth output to a specific directory
mkdir -p "${OUTPUT_BASE_DIR}/${sample_name}_mosdepth"
mv "${OUTPUT_BASE_DIR}/${sample_name}_mosdepth/${sample_name}".* "${OUTPUT_BASE_DIR}/${sample_name}_mosdepth/"

# Step 3: Run spectre.py CNVCaller
python3 spectre.py CNVCaller --bin-size 1000000 --coverage "${OUTPUT_BASE_DIR}/${sample_name}_mosdepth/" --sample-id "${sample_name}" --output-dir "${OUTPUT_BASE_DIR}/${sample_name}_res/" --reference "${REFERENCE_GENOME}"

echo "Analysis completed for ${sample_name}"
