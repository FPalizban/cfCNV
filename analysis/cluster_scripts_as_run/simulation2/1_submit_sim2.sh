#!/bin/bash
# =====================================================================
# 1_submit_sim2.sh  - submits ALL Simulation 2 jobs (login node):   bash 1_submit_sim2.sh
# Dependent steps wait automatically. Job numbers -> logs/simulation2/sim2_jobids.txt
#   A1 QDNAseq 624 | A2 CNVkit 624 | A3 ichorCNA 312 | A4 Spectre 624
#   A5 ACE 312 (after A1) | A6 dudeML 104 -> 624 | A7 CNV-MEANN 624 -> 6
#   B  BIC-seq2: resources (1 long job) + 104 -> 624 | C LYCEUM 104
#   D1 scoring after A + C ; D2 scoring again after BIC-seq2
# =====================================================================
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/simulation2/full/s2_env.sh
cd "$S2"
last() { echo $(( $(wc -l < "$S2_META/tasks_$1.tsv") - 2 )); }     # last array index
sub() { sbatch --parsable "$@"; }
a1=$(sub --array=0-$(last qdnaseq)%40 A1_qdnaseq.slurm)
a2=$(sub --array=0-$(last cnvkit)%40 A2_cnvkit.slurm)
a3=$(sub --array=0-$(last ichor)%10 A3_ichorcna.slurm)
a4=$(sub --array=0-$(last spectre)%30 A4_spectre.slurm)
a5=$(sub --dependency=afterany:$a1 --array=0-$(last ace)%30 A5_ace.slurm)
a6a=$(sub --array=0-$(last dudeml_feat)%20 A6a_dudeml_features.slurm)
a6b=$(sub --dependency=afterany:$a6a --array=0-$(last dudeml_pred)%40 A6b_dudeml_predict.slurm)
a7a=$(sub --array=0-$(last meann_feat)%40 A7a_meann_features.slurm)
a7b=$(sub --dependency=afterany:$a7a A7b_meann_apply.slurm)
b0=$(sub B0_bicseq2_resources.slurm)
b1=$(sub --array=0-$(last bic_unique)%20 B1_bicseq2_unique.slurm)
b2=$(sub --dependency=afterok:$b0,afterany:$b1 --array=0-$(last bic_norm)%40 B2_bicseq2_normseg.slurm)
c1=$(sub --array=0-$(last lyceum)%10 C1_lyceum.slurm)
d1=$(sub --dependency=afterany:$a1:$a2:$a3:$a4:$a5:$a6b:$a7b:$c1 D_score_sim2.slurm)
d2=$(sub --dependency=afterany:$b2:$d1 D_score_sim2.slurm)
cat > "$S2_LOG/sim2_jobids.txt" <<EOF
qdnaseq=$a1 cnvkit=$a2 ichor=$a3 spectre=$a4 ace=$a5 dude_feat=$a6a dude_pred=$a6b meann_feat=$a7a meann_apply=$a7b bic_res=$b0 bic_unique=$b1 bic_seg=$b2 lyceum=$c1 score1=$d1 score2=$d2
EOF
cat "$S2_LOG/sim2_jobids.txt"; echo "Submitted. Progress:  bash 2_check_sim2.sh"
