#!/bin/bash
#SBATCH --job-name=cfCNV_spikeQC
#SBATCH --cpus-per-task=20
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/spikeQC_%A_%a.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/spikeQC_%A_%a.err

set -euo pipefail

source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/config.sh

module load SAMtools/1.16.1-GCC-11.3.0
module load BEDTools/2.30.0-GCC-11.3.0

MAN=$SIM/healthy_cfDNA_spikein/manifests/positive_manifest.tsv

LINE=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+2)) 'NR==n' "$MAN")

IFS=$'\t' read -r ID BACKGROUND CLASS REP DEPTH BAM TRUTH <<< "$LINE"

NEG=$SIM/healthy_cfDNA_spikein/negative_controls/${BACKGROUND}/${BACKGROUND}_${DEPTH}x.noCNV.bam

OUTDIR=$QC/spike_validation
mkdir -p "$OUTDIR"

OUT=$OUTDIR/${ID}.signal.tsv

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

grep -v '^#' "$TRUTH" > "$TMP/truth.bed"

bedtools multicov \
-bed "$TMP/truth.bed" \
-bams "$BAM" \
> "$TMP/positive.tsv"

bedtools multicov \
-bed "$TMP/truth.bed" \
-bams "$NEG" \
> "$TMP/negative.tsv"

python3 - "$TMP/positive.tsv" "$TMP/negative.tsv" "$OUT" <<'PY'
import sys

pfile,nfile,outfile=sys.argv[1:]

with open(pfile) as p, open(nfile) as n, open(outfile,"w") as out:

    out.write(
        "chrom\tstart\tend\tcnv_type\tcopy_number\t"
        "size_class\tcnv_id\tcnv_size\t"
        "positive_reads\tnegative_reads\tread_ratio\n"
    )

    for lp,ln in zip(p,n):

        pfields=lp.rstrip().split("\t")
        nfields=ln.rstrip().split("\t")

        pos=int(pfields[-1])
        neg=int(nfields[-1])

        if neg > 0:
            ratio=pos/neg
            ratio_string=f"{ratio:.6f}"
        else:
            ratio_string="NA"

        out.write(
            "\t".join(
                pfields[:8] +
                [str(pos),str(neg),ratio_string]
            ) +
            "\n"
        )
PY

echo "$ID completed"
