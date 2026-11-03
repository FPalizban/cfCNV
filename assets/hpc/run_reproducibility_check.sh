#!/bin/bash
# cfCNV reproducibility check on the CHOP cluster: re-runs the published containers on the paper's 0.1x GSE71378 BAMs
# (IC37, IH03) and compares calls and tumour fractions with the values reported in the paper.
# EVERYTHING is written inside RUN (your own project folder): images, Nextflow's files, temporary files, work, results.
# Inputs (BAMs, reference, BIC-seq2 image and resources) are only read.
#   bash run_reproducibility_check.sh            -> downloads images, submits one Slurm job, prints the job id
#   when the job has finished:  cat /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfcnv_run/check_result.txt
set -euo pipefail
MY=/mnt/isilon/cag_inf_op/palizbanf/cfCNV
RUN=$MY/cfcnv_run
P=$MY/cfCNV_project
B=$P/04_aligned/GSE71378
mkdir -p "$RUN"/{sif,nxf_home,tmp,work,results}
export NXF_HOME=$RUN/nxf_home NXF_TEMP=$RUN/tmp TMPDIR=$RUN/tmp NXF_SINGULARITY_CACHEDIR=$RUN/sif
export SINGULARITY_TMPDIR=$RUN/tmp APPTAINER_TMPDIR=$RUN/tmp SINGULARITY_CACHEDIR=$RUN/tmp APPTAINER_CACHEDIR=$RUN/tmp
module load singularity/4.4.1 2>/dev/null || true
command -v java >/dev/null || { echo "Java not found: load your Java module first (module load ...)"; exit 1; }
JH=$(dirname "$(dirname "$(readlink -f "$(command -v java)")")")
NF=$(command -v nextflow || echo "$HOME/bin/nextflow")

echo "== 1. inputs (read only)"
ls -l "$B/IC37/IC37_d0p1.bam" "$B/IC37/IC37_d0p1.bam.bai" "$B/IH03/IH03_d0p1.bam" "$B/IH03/IH03_d0p1.bam.bai"
ls -l "$P/02_reference/GRCh37/fasta/human_g1k_v37.fasta.fai" "$P/06_tools/containers/bicseq2_20210625.sif"

echo "== 2. container images (ready-made .sif files; nothing is unpacked)"
bash "$(dirname "$0")/pull_sif.sh" "$RUN/sif" 1.0.0
n=$(ls "$RUN"/sif/cfcnv-*-1.0.0.sif 2>/dev/null | wc -l)
[[ $n -eq 7 ]] || { echo "only $n of 7 images downloaded - check that all cfcnv packages are public, then run this again"; exit 2; }

echo "== 3. workflow code (pinned to the GitHub repository)"
"$NF" pull FPalizban/cfCNV -r main

printf "sample,bam,group,tcga_cohort\nIC37,%s,cancer,COAD+READ\nIH03,%s,control,\n" \
  "$B/IC37/IC37_d0p1.bam" "$B/IH03/IH03_d0p1.bam" > "$RUN/samples.csv"

cat > "$RUN/check.slurm" <<EOF
#!/bin/bash
#SBATCH --job-name=cfCNV_check
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=$RUN/check_%j.log
set -euo pipefail
cd $RUN
export NXF_HOME=$RUN/nxf_home NXF_TEMP=$RUN/tmp TMPDIR=$RUN/tmp NXF_SINGULARITY_CACHEDIR=$RUN/sif
export SINGULARITY_TMPDIR=$RUN/tmp APPTAINER_TMPDIR=$RUN/tmp SINGULARITY_CACHEDIR=$RUN/tmp APPTAINER_CACHEDIR=$RUN/tmp
export JAVA_HOME=$JH PATH=$JH/bin:\$PATH NXF_OPTS='-Xms1g -Xmx4g'
module load singularity/4.4.1 2>/dev/null || true
$NF run FPalizban/cfCNV -r main -profile chop -w $RUN/work \\
    --input $RUN/samples.csv --outdir $RUN/results --genome GRCh37 \\
    --fasta $P/02_reference/GRCh37/fasta/human_g1k_v37.fasta \\
    --bicseq2_resources $P/06_tools/BICseq2_resources/GRCh37_101mer --bicseq2_map_prefix human_g1k_v37.101mer
A=$RUN/nxf_home/assets/FPalizban/cfCNV/analysis
singularity exec --home $RUN $RUN/sif/cfcnv-tools-1.0.0.sif python3 \$A/compare_with_paper.py \\
    $RUN/results/summary/all_calls.tsv \$A/real_plasma/paper_calls_GSE71378_primary.tsv.gz \\
    --map IC37=IC37_d0p1 IH03=IH03_d0p1 \\
    --tf $RUN/results/summary/tumour_fraction.tsv --paper-tf \$A/real_plasma/paper_tf_GSE71378.tsv \\
    | tee $RUN/check_result.txt
EOF

echo "== 4. submitting the check (one small Slurm job; Nextflow submits each tool step as its own job)"
sbatch "$RUN/check.slurm"
echo
echo "Follow it with:   squeue -u $USER"
echo "When it is gone from squeue:   cat $RUN/check_result.txt    (if empty: tail -40 $RUN/check_*.log)"
