#!/bin/bash
#SBATCH --job-name=dudeml_check
#SBATCH --cpus-per-task=20
#SBATCH --mem=120G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/dudeml_preflight_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/dudeml_preflight_%j.err

set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project

source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

TOOLS="$PROJECT/06_tools"
DUDE="$TOOLS/ML_based/dudeML"
ENV="$TOOLS/envs/dudeml"

mkdir -p "$PROJECT/logs/simulation1"
mkdir -p "$PROJECT/07_tool_outputs/simulation1/dudeML"

source "$ENV/bin/activate"

echo "=========================================="
echo "dudeML PREFLIGHT"
echo "=========================================="

python --version

echo
echo "===== REPOSITORY ====="
echo "$DUDE"

cd "$DUDE"

git rev-parse HEAD 2>/dev/null || true

echo
echo "===== MAIN HELP ====="
python dudeML.py -h || true

for CMD in \
    winStat \
    fvecTrain \
    fvecSample \
    classify \
    predict \
    summarize \
    simCNV \
    recreateTotal
do
    echo
    echo "===== $CMD ====="
    python dudeML.py "$CMD" -h || true
done

echo
echo "===== BEDTOOLS ====="
command -v bedtools || true
bedtools --version || true

echo
echo "===== BENCHMARK MATRIX ====="
echo "samples = $(( $(wc -l < "$SIM1_MAN/master.tsv") - 1 ))"
echo "QDNAseq-style sample/bin tasks = $(( $(wc -l < "$SIM1_MAN/qdnaseq_tasks.tsv") - 1 ))"

echo
echo "===== BIN DISTRIBUTION ====="

awk -F'\t' '
NR==1 {
    for(i=1;i<=NF;i++)
        print i,$i
}
' "$SIM1_MAN/qdnaseq_tasks.tsv"

echo
echo "=========================================="
echo "dudeML PREFLIGHT COMPLETE"
echo "=========================================="
