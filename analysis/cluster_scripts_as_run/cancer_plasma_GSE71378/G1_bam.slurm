#!/bin/bash
#SBATCH --job-name=G7_bam
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/gse71378/bam_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/gse71378/bam_%A_%a.err
# G1  One sample, NO realignment: the authors' GRCh37 BAM (kircherlab mirror; IH03 re-used from
#     Simulation 1) -> depth (mosdepth) -> analysis BAMs: full (down-sampled to ~12x when deeper),
#     1x and 0.1x. Unmapped, secondary, supplementary and duplicate-flagged reads are dropped.
set -euo pipefail
trap 'echo "FAILED at line $LINENO: $BASH_COMMAND (exit $?)" >&2' ERR
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/gse71378/full/p3_env.sh
S=$(grep -vE '^#|^\s*$' "$S2/samples.txt" | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" | awk '{print $1}')
OD="$P3_BAM/$S"; mkdir -p "$OD"
[[ -f "$OD/.derive_done" ]] && { echo "SKIP $S"; exit 0; }
aln_modules
# 1. source BAM
SRC_IH="$P/01_raw_data/shendure_cfDNA/original_BAM/$S.bam"
if [[ -s "$SRC_IH" && -s "$SRC_IH.bai" ]]; then SRC="$SRC_IH"
else
  SRC="$P3_FASTQ/$S/$S.mirror.bam"; mkdir -p "$(dirname "$SRC")"
  until wget -c -q --tries=0 --waitretry=30 --read-timeout=120 -O "$SRC" "$MIRROR/$S.bam"; do echo "resuming download"; sleep 30; done
  wget -q -O "$SRC.bai" "$MIRROR/$S.bam.bai"
  [[ "$SRC.bai" -nt "$SRC" ]] || touch "$SRC.bai"
fi
samtools quickcheck "$SRC" || { echo "BAM incomplete: delete $SRC and resubmit"; exit 2; }
echo "source $SRC ($(du -h "$SRC" | cut -f1)); chromosome names: $(samtools view -H "$SRC" | awk '$1=="@SQ"{sub("SN:","",$2); print $2}' | head -3 | tr '\n' ' ')"
# 2. depth of the source BAM
[[ -s "$OD/${S}_source.mosdepth.summary.txt" ]] || "$MOSDEPTH" -t 4 -n --fast-mode "$OD/${S}_source" "$SRC"
DEP=$(awk '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2])$/ {L+=$2; B+=$3} END{printf "%.4f", B/L}' "$OD/${S}_source.mosdepth.summary.txt")
NP=$(set +o pipefail; samtools view -F 0x904 "$SRC" 2>/dev/null | head -200000 | awk 'int($2/1)%2==1{p++} END{print p+0}')   # head closes the pipe early: pipefail off here
LAYOUT=$( (( NP > 100000 )) && echo PAIRED || echo SINGLE )
echo "source depth ${DEP}x, layout $LAYOUT (paired reads in first 200k: $NP)"
# 3. analysis BAMs
mk() {  # mk <fraction> <out>
  if awk -v f="$1" 'BEGIN{exit !(f >= 0.999)}'; then samtools view -@ 8 -b -F 0xD04 -o "$2" "$SRC"
  else samtools view -@ 8 -b -F 0xD04 -s "11${1#0}" -o "$2" "$SRC"; fi
  samtools index -@ 8 "$2"; }
FULLX=$(awk -v d="$DEP" -v t="$TARGET_FULL_X" 'BEGIN{print (d > t) ? t : d}')
FR=$(awk -v d="$DEP" -v x="$FULLX" 'BEGIN{printf "%.6f", x/d}')
BAM="$OD/$S.bam"; [[ -s "$BAM.bai" ]] || mk "$FR" "$BAM"
samtools stats -@ 8 "$BAM" > "$OD/$S.stats.txt"
[[ -s "$OD/$S.mosdepth.summary.txt" ]] || "$MOSDEPTH" -t 4 -n --fast-mode "$OD/$S" "$BAM"
DF=$(awk '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2])$/ {L+=$2; B+=$3} END{printf "%.4f", B/L}' "$OD/$S.mosdepth.summary.txt")
{
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${S}_full" "$S" "$BAM" "$DF" full "$LAYOUT"
  for T in 1 0.1; do
    if awk -v d="$DF" -v t="$T" 'BEGIN{exit !(d > 1.3*t)}'; then
      L=$( [[ $T == 1 ]] && echo d1 || echo d0p1 ); X="$OD/${S}_${L}.bam"
      F=$(awk -v d="$DF" -v t="$T" 'BEGIN{printf "%.6f", t/d}')
      [[ -s "$X.bai" ]] || { samtools view -@ 4 -b -s "11${F#0}" -o "$X" "$BAM" && samtools index "$X"; }
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${S}_${L}" "$S" "$X" "$T" "$L" "$LAYOUT"
    fi
  done
} > "$OD/derived.tsv"
# the downloaded source BAM is kept until the run is checked (delete later to free space)
touch "$OD/.align_done" "$OD/.derive_done"; cat "$OD/derived.tsv"; echo "BAM_DONE $S"
