#!/bin/bash
#SBATCH --job-name=ichorCNA_array_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=array_job_%A_%a.out
#SBATCH --error=array_job_%A_%a.err
#SBATCH --array=1-10%5  # Adjust the range and step size based on the number of BAM files

# Load necessary modules
module load gcc samtools
module load gcc r

# Define paths and parameters
BAM_DIR="/work/FAC/FBM/DBC/cdessim2/default/fpalizba/cfDNA-CNV/validation/simulated-samples"
OUTPUT_DIR="/work/FAC/FBM/DBC/cdessim2/default/fpalizba/cfDNA-CNV/validation/simulated-samples/run-results/ichorCNA-1000000-IH03"

# Get the BAM file corresponding to the array task index
bam_file=$(ls -1 ${BAM_DIR}/*.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)
sample_name=$(basename "${bam_file%.*}")

# Step 1: Run readCounter
./hmmcopy_utils/bin/readCounter --window 1000000 --quality 10 --chromosome "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y" "${bam_file}" > "${OUTPUT_DIR}/${sample_name}.wig"

# Step 2: Run IchorCNA
Rscript scripts/runIchorCNA.R --id "${sample_name}" --WIG "${OUTPUT_DIR}/${sample_name}.wig" --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 --gcWig inst/extdata/gc_hg19_1000kb.wig --mapWig inst/extdata/map_hg19_1000kb.wig --centromere inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt --normalPanel inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" --estimateNormal True --estimatePloidy True --estimateScPrevalence True --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir "${OUTPUT_DIR}"

echo "Analysis completed for ${sample_name}"
