#!/bin/bash
#SBATCH --job-name=cfcnv_score_s1
#SBATCH --cpus-per-task=20
#SBATCH --mem=150G
#SBATCH --time=08:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/scoring_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/scoring_%j.err
set -euo pipefail
PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project
source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"
mkdir -p "$PROJECT/logs/simulation1" "$PROJECT/09_metrics/simulation1"
FAI=""
for C in "${REF37_FASTA:-}.fai" "$REF37_DIR/human_g1k_v37.fasta.fai" "$REF37_DIR/fasta/human_g1k_v37.fasta.fai"; do
  if [[ -n "$C" && -s "$C" ]]; then FAI="$C"; break; fi
done
[[ -n "$FAI" ]] || { echo "ERROR: GRCh37 .fai not found" >&2; exit 2; }
python3 "$PROJECT/00_metadata/50_score_sim1.py" \
  --sim1-out "$SIM1_OUT" \
  --sim1-man "$SIM1_MAN" \
  --truth-dir "$SIM/healthy_cfDNA_spikein/truth" \
  --fai "$FAI" \
  --outdir "$PROJECT/09_metrics/simulation1"
echo "===== SCORING SUMMARY ====="
cat "$PROJECT/09_metrics/simulation1/SCORING_SUMMARY.txt"
