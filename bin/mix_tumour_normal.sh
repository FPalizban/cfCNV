#!/usr/bin/env bash
# Simulation 2: in-silico tumour-normal mixture at a target depth and tumour fraction, identical to the
# cfCNV paper (HG008-T pooled to 4x, HG008-N to 12x): p_tumour = depth*f/Dt, p_normal = depth*(1-f)/Dn.
# usage: mix_tumour_normal.sh <tumour.bam> <Dt> <normal.bam> <Dn> <depth> <fraction> <seed> <out.bam> <threads>
set -euo pipefail
TB=$1; DT=$2; NB=$3; DN=$4; D=$5; F=$6; SEED=$7; OUT=$8; T=${9:-2}
read -r PT PN < <(python3 - "$D" "$F" "$DT" "$DN" <<'PY'
import sys
d, f, dt, dn = map(float, sys.argv[1:5])
pt, pn = d * f / dt, d * (1 - f) / dn
assert 0 <= pt < 1 and 0 < pn < 1, "requested depth exceeds the depth of a pool"
print(f"{pt:.9f} {pn:.9f}")
PY
)
samtools view -@ "$T" -b -s "${SEED}.${PN#0.}" -o normal.part.bam "$NB"
if [[ "$PT" != "0.000000000" ]]; then
  samtools view -@ "$T" -b -s "$((SEED+100)).${PT#0.}" -o tumour.part.bam "$TB"
  samtools merge -@ "$T" -f unsorted.bam tumour.part.bam normal.part.bam
  samtools sort -@ "$T" -o "$OUT" unsorted.bam
  rm -f tumour.part.bam unsorted.bam
else
  samtools sort -@ "$T" -o "$OUT" normal.part.bam
fi
rm -f normal.part.bam
samtools index -@ "$T" "$OUT"
samtools quickcheck -v "$OUT"
