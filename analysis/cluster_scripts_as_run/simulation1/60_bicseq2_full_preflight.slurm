#!/bin/bash
#SBATCH --job-name=bicseq2_check
#SBATCH --cpus-per-task=20
#SBATCH --mem=160G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/bicseq2_preflight_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/bicseq2_preflight_%j.err

set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project

source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

mkdir -p "$PROJECT/logs/simulation1"
mkdir -p "$PROJECT/07_tool_outputs/simulation1/BICseq2"
mkdir -p "$PROJECT/06_tools/BICseq2_resources/GRCh37"

TOOLS="$PROJECT/06_tools"

WRAPPER="$TOOLS/bin/run_bicseq2"
IMAGE="$TOOLS/containers/bicseq2_20210625.sif"
REPO="$TOOLS/classical/BICSEQ2"

SINGULARITY=/cm/shared/apps_chop/singularity/4.4.1/bin/singularity

echo "===================================================="
echo "cfCNV BIC-seq2 FULL PREFLIGHT"
echo "===================================================="
echo "DATE      : $(date)"
echo "HOST      : $(hostname)"
echo "PROJECT   : $PROJECT"
echo "WRAPPER   : $WRAPPER"
echo "IMAGE     : $IMAGE"
echo "REPO      : $REPO"
echo

echo "===== REQUIRED FILES ====="

for X in \
    "$WRAPPER" \
    "$IMAGE" \
    "$REPO"
do
    if [[ -e "$X" ]]; then
        echo "FOUND: $X"
    else
        echo "MISSING: $X"
    fi
done

echo
echo "===== WRAPPER ====="

if [[ -f "$WRAPPER" ]]; then
    sed -n '1,300p' "$WRAPPER"
fi

echo
echo "===== REPOSITORY FILES ====="

if [[ -d "$REPO" ]]; then
    find "$REPO" \
        -maxdepth 4 \
        -type f \
        \( \
          -iname '*bic*' \
          -o -iname '*case*' \
          -o -iname '*norm*' \
          -o -iname '*detect*' \
          -o -iname '*mapp*' \
          -o -iname '*.cwl' \
          -o -iname '*.yml' \
          -o -iname '*.yaml' \
          -o -iname '*.json' \
        \) \
        | sort \
        | head -300
fi

echo
echo "===== CONTAINER VERSION ====="

"$SINGULARITY" --version || true

echo
echo "===== CONTAINER CONTENTS ====="

"$SINGULARITY" exec "$IMAGE" \
    sh -lc '
      echo "--- BIC executables/scripts ---"
      find / \
        -maxdepth 5 \
        -type f \
        \( \
          -iname "NBICseq*.pl" \
          -o -iname "*bicseq*" \
          -o -iname "*mappability*" \
          -o -iname "get_unique*" \
          -o -iname "run_norm*" \
          -o -iname "run_detect*" \
        \) \
        2>/dev/null | sort | head -300

      echo
      echo "--- Important commands ---"
      command -v perl || true
      command -v samtools || true
      command -v gem-indexer || true
      command -v gem-mappability || true
      command -v gem-2-wig || true
    ' || true

echo
echo "===== LOCAL GRCh37 REFERENCE ====="

find "$REF37_DIR" \
    -maxdepth 3 \
    -type f \
    \( \
      -name '*.fa' \
      -o -name '*.fasta' \
      -o -name '*.fai' \
    \) \
    -printf '%s\t%p\n' \
    2>/dev/null \
    | sort -nr \
    | head -50

echo
echo "===== EXISTING BIC/MAPPABILITY RESOURCES ====="

find "$PROJECT" \
    -type f \
    \( \
      -iname '*mer*.txt' \
      -o -iname '*mappability*' \
      -o -iname '*.gem' \
    \) \
    2>/dev/null \
    | head -200

echo
echo "===== MASTER MANIFEST ====="

head -5 "$SIM1_MAN/master.tsv"

echo
echo "Samples:"
echo $(( $(wc -l < "$SIM1_MAN/master.tsv") - 1 ))

echo
echo "===== BAM CONTIG CHECK ====="

FIRST_BAM=$(awk -F'\t' 'NR==2 {print $2}' "$SIM1_MAN/master.tsv")

echo "Example BAM:"
echo "$FIRST_BAM"

"$TOOLS/bin/run_bicseq2" samtools view -H "$FIRST_BAM" \
    | grep '^@SQ' \
    | head -5 || true

echo
echo "===== READ LENGTH CHECK ====="

set +o pipefail

"$TOOLS/bin/run_bicseq2" samtools view "$FIRST_BAM" \
    | awk '
      NR <= 10000 {
          L[length($10)]++
      }
      NR == 10000 {
          exit
      }
      END {
          for (x in L)
              print x, L[x]
      }
    ' \
    | sort -k2,2nr \
    | head -10

set -o pipefail

echo
echo "===== AVAILABLE WORKFLOW ENGINES ====="

command -v cwltool || true
command -v cwl-runner || true
command -v python3 || true

if command -v cwltool >/dev/null 2>&1; then
    cwltool --version || true
fi

echo
echo "===================================================="
echo "BICSEQ2 PREFLIGHT COMPLETE"
echo "===================================================="
