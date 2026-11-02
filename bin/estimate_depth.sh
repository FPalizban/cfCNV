#!/usr/bin/env bash
# Mean autosomal depth of a BAM: mapped autosomal reads x median read length / autosomal length (from the BAM header).
set -euo pipefail
BAM=$1
read -r READS LEN < <(samtools idxstats "$BAM" | awk '{c=$1; sub(/^chr/,"",c)} c ~ /^([1-9]|1[0-9]|2[0-2])$/ {r+=$3; l+=$2} END{print r, l}')
RL=$( (set +o pipefail; samtools view -F 0x904 "$BAM" | head -n 20000 | awk '{print length($10)}' | sort -n | awk '{a[NR]=$1} END{print (NR? a[int((NR+1)/2)] : 100)}') )
awk -v r="$READS" -v rl="$RL" -v l="$LEN" 'BEGIN{printf "%.4f\n", (l>0 ? r*rl/l : 0)}'
