#!/usr/bin/env bash
# Down-sample a BAM to a target mean autosomal depth.
# Depth = mapped autosomal reads x median read length / autosomal genome size.
# usage: downsample_bam.sh <in.bam> <target_x> <genome_size_bp|fai> <out.bam> <seed> <threads>
# writes <out.bam>.depth.txt : source_depth  target_depth  fraction
set -euo pipefail
BAM=$1; TARGET=$2; G=$3; OUT=$4; SEED=${5:-11}; T=${6:-2}
if [[ -f "$G" ]]; then
  G=$(awk '{c=$1; sub(/^chr/,"",c)} c ~ /^([1-9]|1[0-9]|2[0-2])$/ {s+=$2} END{print s}' "$G")
fi
READS=$(samtools idxstats "$BAM" | awk '{c=$1; sub(/^chr/,"",c)} c ~ /^([1-9]|1[0-9]|2[0-2])$/ {s+=$3} END{print s}')
RL=$( (set +o pipefail; samtools view -F 0x904 "$BAM" | head -n 20000 | awk '{print length($10)}' | sort -n \
      | awk '{a[NR]=$1} END{print (NR? a[int((NR+1)/2)] : 100)}') )
read -r DEPTH FRAC < <(python3 - "$READS" "$RL" "$G" "$TARGET" <<'PY'
import sys
r, rl, g, t = map(float, sys.argv[1:5])
d = r * rl / g
print(f"{d:.4f} {t / d if d else 1:.8f}")
PY
)
echo -e "source_depth\ttarget_depth\tfraction\n$DEPTH\t$TARGET\t$FRAC" > "$OUT.depth.txt"
if python3 -c "import sys; sys.exit(0 if float('$FRAC') >= 0.999 else 1)"; then
  echo "target ${TARGET}x >= source ${DEPTH}x: copying the input unchanged" >&2
  cp "$BAM" "$OUT"
else
  F=$(python3 -c "print(f'{float(\"$FRAC\"):.8f}'.split('.')[1])")
  samtools view -@ "$T" -b -s "${SEED}.${F}" -o "$OUT" "$BAM"
fi
samtools index -@ "$T" "$OUT"
