#!/bin/bash
#SBATCH --job-name=HG008_mix
#SBATCH --cpus-per-task=14
#SBATCH --mem=160G
#SBATCH --time=08:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/mix_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation2/mix_%A_%a.err

set -euo pipefail
module load SAMtools/1.16.1-GCC-11.3.0

P=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project
ROOT="$P/05_simulation/simulation2_HG008"
MAN="$P/00_metadata/simulation2/HG008_mixtures.tsv"

LINE=$(awk -F '\t' -v n="$((SLURM_ARRAY_TASK_ID+2))" \
  'NR==n {print $0}' "$MAN")

if [[ -z "$LINE" ]]; then
    echo "Missing mixture manifest row" >&2
    exit 2
fi

IFS=$'\t' read -r ID DEPTH F <<< "$LINE"

TPOOL="$ROOT/pools/HG008_T_4x.bam"
NPOOL="$ROOT/pools/HG008_N_D_12x.bam"
OUT="$ROOT/mixtures/$ID.bam"

if [[ -s "$OUT" && -s "$OUT.bai" ]]; then
    samtools quickcheck -v "$OUT"
    echo "MIXTURE_REUSED $ID"
    exit 0
fi

read -r PT PN < <(python3 - "$DEPTH" "$F" <<'PY'
import sys
d=float(sys.argv[1])
f=float(sys.argv[2])
pt=d*f/4
pn=d*(1-f)/12
assert 0 <= pt < 1 and 0 < pn < 1
print(f"{pt:.9f} {pn:.9f}")
PY
)

TMPDIR="$ROOT/tmp/$ID"
mkdir -p "$TMPDIR"

T="$TMPDIR/tumor.bam"
N="$TMPDIR/normal.bam"
COMBINED="$TMPDIR/mixture.bam"

SEED=$((1000+SLURM_ARRAY_TASK_ID))

samtools view -@ 2 -b \
  -s "${SEED}.${PN#0.}" \
  -o "$N" "$NPOOL"

if [[ "$PT" != "0.000000000" ]]; then
    samtools view -@ 2 -b \
      -s "$((SEED+100)).${PT#0.}" \
      -o "$T" "$TPOOL"

    samtools merge -@ 4 -f \
      "$COMBINED" "$T" "$N"
else
    mv "$N" "$COMBINED"
fi

samtools quickcheck -v "$COMBINED"
samtools index -@ 4 "$COMBINED"

mv "$COMBINED" "$OUT"
mv "$COMBINED.bai" "$OUT.bai"

rm -f "$T" "$N"
rmdir "$TMPDIR"

printf '%s\t%s\t%s\t%s\n' \
  "$ID" "$DEPTH" "$F" "$OUT" \
  > "$ROOT/mixtures/$ID.done"

echo "MIXTURE_READY $ID"
