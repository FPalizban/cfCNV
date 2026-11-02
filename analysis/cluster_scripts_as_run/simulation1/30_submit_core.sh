#!/bin/bash
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/00_env.sh
python3 "$PROJECT/00_metadata/01_prepare_manifest.py"
NQ=$(($(wc -l<"$SIM1_MAN/qdnaseq_tasks.tsv")-1));NC=$(($(wc -l<"$SIM1_MAN/cnvkit_tasks.tsv")-1))
NI=$(($(wc -l<"$SIM1_MAN/ichor_tasks.tsv")-1));NS=$(($(wc -l<"$SIM1_MAN/spectre_tasks.tsv")-1))
JR=$(sbatch --parsable "$PROJECT/00_metadata/02_prepare_reference.sh")
JQ=$(sbatch --parsable --dependency=afterok:$JR --array=0-$((NQ-1))%18 "$PROJECT/00_metadata/10_qdnaseq_array.sh")
JC=$(sbatch --parsable --dependency=afterok:$JR --array=0-$((NC-1))%18 "$PROJECT/00_metadata/11_cnvkit_array.sh")
JI=$(sbatch --parsable --dependency=afterok:$JR --array=0-$((NI-1))%12 "$PROJECT/00_metadata/12_ichorcna_array.sh")
JS=$(sbatch --parsable --dependency=afterok:$JR --array=0-$((NS-1))%12 "$PROJECT/00_metadata/13_spectre_array.sh")
printf "reference\t%s\nQDNAseq\t%s\nCNVkit\t%s\nichorCNA\t%s\nSpectre\t%s\n" "$JR" "$JQ" "$JC" "$JI" "$JS" | tee "$SIM1_MAN/submitted_jobs.tsv"
