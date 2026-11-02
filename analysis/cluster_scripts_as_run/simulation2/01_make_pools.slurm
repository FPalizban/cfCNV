#!/bin/bash
#SBATCH --job-name=HG008_pool
#SBATCH --cpus-per-task=12
#SBATCH --mem=160G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/pool_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/pool_%A_%a.err

set -euo pipefail
module load SAMtools/1.16.1-GCC-11.3.0

P=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project
D="$P/05_simulation/simulation2_HG008/pools"

BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/Element_AVITI_20240626"

if [[ "$SLURM_ARRAY_TASK_ID" == 0 ]]; then
    SAMPLE=HG008_T
    SOURCE_DEPTH=102
    TARGET_DEPTH=4
    SEED=51
    FILE=HG008-T_Element-StdInsert_102x_GRCh38-GIABv3.bam
else
    SAMPLE=HG008_N_D
    SOURCE_DEPTH=72
    TARGET_DEPTH=12
    SEED=52
    FILE=HG008-N-D_Element-StdInsert_72x_GRCh38-GIABv3.bam
fi

OUT="$D/${SAMPLE}_${TARGET_DEPTH}x.bam"

if [[ -s "$OUT" && -s "$OUT.bai" ]]; then
    samtools quickcheck -v "$OUT"
    echo "POOL_REUSED $OUT"
    exit 0
fi

PROB=$(python3 - "$TARGET_DEPTH" "$SOURCE_DEPTH" <<'PY'
import sys
p=float(sys.argv[1])/float(sys.argv[2])
assert 0 < p < 1
print(f"{p:.9f}")
PY
)

TMP="$OUT.pending.bam"
rm -f "$TMP" "$TMP.bai"

LOCAL=$(find "$P/01_raw_data/NIST_GIAB_HG008" \
    -type f -name "$FILE" -print -quit)

echo "SOURCE=$FILE"
echo "SUBSAMPLE_PROBABILITY=$PROB"

if [[ -n "$LOCAL" && -s "$LOCAL" ]]; then
    samtools view -@ 2 -b -s "${SEED}.${PROB#0.}" "$LOCAL" \
      | samtools sort -@ 4 -m 1G -o "$TMP" -
else
    curl -fL --retry 4 --retry-delay 15 "$BASE/$FILE" \
      | samtools view -@ 2 -b -s "${SEED}.${PROB#0.}" - \
      | samtools sort -@ 4 -m 1G -o "$TMP" -
fi

samtools quickcheck -v "$TMP"
samtools index -@ 4 "$TMP"

mv "$TMP" "$OUT"
mv "$TMP.bai" "$OUT.bai"

echo "POOL_READY $OUT"
