#!/bin/bash
#SBATCH --job-name=sort_index_bam
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=sort_index_bam.out
#SBATCH --error=sort_index_bam.err

# Load necessary modules
module load gcc samtools

# Define input and output file paths
input_bam="/work/FAC/FBM/DBC/cdessim2/default/fpalizba/cfDNA-CNV/validation/IH01_simulated.bam"
sorted_bam="/work/FAC/FBM/DBC/cdessim2/default/fpalizba/cfDNA-CNV/validation/IH01_simulated_sorted.bam"
#indexed_bam="/path/to/your/sorted_output.bam.bai"

# Step 1: Sort BAM file
samtools sort -@ 4 -o "$sorted_bam" "$input_bam"

# Step 2: Index sorted BAM file
samtools index "$sorted_bam" 

echo "Sorting and indexing completed successfully."
