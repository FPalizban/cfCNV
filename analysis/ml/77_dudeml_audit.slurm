#!/bin/bash
#SBATCH --job-name=dude_audit
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/dudeML/dude_audit_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/dudeML/dude_audit_%j.err

set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project

source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

DUDE="$PROJECT/07_tool_outputs/simulation1/dudeML"

echo "======================================"
echo "dudeML Simulation-1 audit"
echo "======================================"

FEATURES=$(find "$DUDE/features" \
    -type f \
    -name '.features_done' \
    | wc -l)

MODELS=$(find "$DUDE/models" \
    -type f \
    -name 'fold*.sav' \
    | wc -l)

MODELS2=$(find "$DUDE/models" \
    -type f \
    -name 'fold*.sav2' \
    | wc -l)

DONE=$(find "$DUDE/predictions" \
    -type f \
    -name '.done' \
    | wc -l)

SEG=$(find "$DUDE/predictions" \
    -type f \
    -name '*.segments.bed' \
    | wc -l)


echo "Feature-complete samples : $FEATURES / 114"
echo "Primary models           : $MODELS / 240"
echo "Ploidy models            : $MODELS2 / 240"
echo "Prediction configurations: $DONE / 684"
echo "Segment files            : $SEG / 684"


if [[ \
"$FEATURES" -eq 114 && \
"$MODELS" -eq 240 && \
"$MODELS2" -eq 240 && \
"$DONE" -eq 684 && \
"$SEG" -eq 684 \
]]
then

    echo
    echo "ALL DUDEML OUTPUTS COMPLETE"

else

    echo
    echo "DUDEML OUTPUTS INCOMPLETE"

fi
