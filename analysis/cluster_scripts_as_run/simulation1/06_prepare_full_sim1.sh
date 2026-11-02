#!/bin/bash
set -euo pipefail

source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/config.sh

module load SAMtools/1.16.1-GCC-11.3.0
module load BEDTools/2.30.0-GCC-11.3.0

SIM1=$SIM/healthy_cfDNA_spikein
MAN=$SIM1/manifests

mkdir -p \
    "$SIM1/parents" \
    "$SIM1/truth" \
    "$SIM1/spiked" \
    "$SIM1/depth_series" \
    "$SIM1/negative_controls" \
    "$SIM1/manifests" \
    "$SIM1/tmp" \
    "$REF37_DIR/bins" \
    "$LOGS"

echo "======================================================"
echo "PREPARING FULL SIMULATION 1"
echo "======================================================"

# ------------------------------------------------------------
# 1. GRCh37 chromosome sizes
# ------------------------------------------------------------

CHROMSIZES=$REF37_DIR/Shendure_GRCh37.chrom.sizes

if [[ ! -s "$CHROMSIZES" ]]
then

    echo "Creating GRCh37 chromosome-size file..."

    samtools view -H "$SHENDURE_DIR/IH01.bam" \
    | awk -F'\t' '
    $1=="@SQ" {
        sn="";
        ln="";

        for(i=1;i<=NF;i++){
            if($i ~ /^SN:/) sn=substr($i,4)
            if($i ~ /^LN:/) ln=substr($i,4)
        }

        if(sn ~ /^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$/)
            print sn"\t"ln
    }' > "$CHROMSIZES"

fi

echo
echo "Chromosome-size file:"
cat "$CHROMSIZES"

# ------------------------------------------------------------
# 2. Benchmark bin sizes
# ------------------------------------------------------------

echo
echo "Creating bin files..."

for BIN in \
10000 \
30000 \
50000 \
100000 \
500000 \
1000000
do

    bedtools makewindows \
        -g "$CHROMSIZES" \
        -w "$BIN" \
        > "$REF37_DIR/bins/GRCh37.${BIN}.bed"

done

cat > "$MAN/bin_sizes.tsv" <<'BINS'
bin_bp
10000
30000
50000
100000
500000
1000000
BINS

# ------------------------------------------------------------
# 3. Experimental design
# ------------------------------------------------------------

cat > "$MAN/sample_design.tsv" <<'DESIGN'
sample	native_depth	parent_depth	depths
IH01	97.656559	20	0.1,0.25,0.5,1,2,5,10,20
IH02	19.416836	10	0.1,0.25,0.5,1,2,5,10
IH03	1.085300	1	0.1,0.25,0.5,1
DESIGN

# ------------------------------------------------------------
# 4. Generate five truth replicates
# ------------------------------------------------------------

python3 <<'PY'
import os
import random
import csv

SIM = os.environ["SIM"]
REF37_DIR = os.environ["REF37_DIR"]
SHENDURE_DIR = os.environ["SHENDURE_DIR"]

root = os.path.join(SIM, "healthy_cfDNA_spikein")
truth_dir = os.path.join(root, "truth")
manifest_dir = os.path.join(root, "manifests")

os.makedirs(truth_dir, exist_ok=True)
os.makedirs(manifest_dir, exist_ok=True)

chromfile = os.path.join(
    REF37_DIR,
    "Shendure_GRCh37.chrom.sizes"
)

# ------------------------------------------------------------
# CNV design
#
# 30 events per replicate.
# Includes both deletions and duplications.
# ------------------------------------------------------------

cnv_design = {
    50_000: 4,
    100_000: 4,
    250_000: 4,
    500_000: 4,
    1_000_000: 4,
    5_000_000: 4,
    10_000_000: 2,
    25_000_000: 2,
    50_000_000: 2
}

chromosomes = []

with open(chromfile) as f:

    for line in f:

        chrom, length = line.rstrip().split()[:2]

        if chrom.isdigit() and 1 <= int(chrom) <= 22:
            chromosomes.append(
                (chrom, int(length))
            )

for rep in range(1, 6):

    random.seed(1000 + rep)

    events = []
    event_number = 1

    def overlaps(chrom, start, end):

        padding = 1_000_000

        for event in events:

            ec, es, ee = event[0], event[1], event[2]

            if ec != chrom:
                continue

            if start < ee + padding and end > es - padding:
                return True

        return False

    for size, number in cnv_design.items():

        for _ in range(number):

            placed = False

            for attempt in range(10000):

                chrom, chrom_length = random.choice(chromosomes)

                margin = max(
                    5_000_000,
                    int(size * 0.2)
                )

                if chrom_length <= size + 2 * margin:
                    continue

                start = random.randint(
                    margin,
                    chrom_length - size - margin
                )

                end = start + size

                if overlaps(chrom, start, end):
                    continue

                cnv_type = (
                    "loss"
                    if event_number % 2
                    else "gain"
                )

                copy_number = (
                    1 if cnv_type == "loss"
                    else 3
                )

                if size <= 250_000:
                    size_class = "small"

                elif size <= 5_000_000:
                    size_class = "medium"

                else:
                    size_class = "large"

                events.append(
                    (
                        chrom,
                        start,
                        end,
                        cnv_type,
                        copy_number,
                        size_class,
                        f"REP{rep}_CNV{event_number}",
                        size
                    )
                )

                event_number += 1
                placed = True
                break

            if not placed:
                raise RuntimeError(
                    f"Could not place CNV size {size}"
                )

    events.sort(
        key=lambda x: (
            int(x[0]),
            x[1]
        )
    )

    output = os.path.join(
        truth_dir,
        f"rep{rep}.truth.tsv"
    )

    with open(output, "w") as out:

        out.write(
            "#chrom\tstart\tend\tcnv_type\t"
            "copy_number\tsize_class\t"
            "cnv_id\tcnv_size\n"
        )

        for event in events:
            out.write(
                "\t".join(
                    map(str, event)
                ) + "\n"
            )

# ------------------------------------------------------------
# Sample/depth design
# ------------------------------------------------------------

samples = {

    "IH01": {
        "native": 97.656559,
        "parent": 20.0,
        "depths": [
            0.1,
            0.25,
            0.5,
            1,
            2,
            5,
            10,
            20
        ]
    },

    "IH02": {
        "native": 19.416836,
        "parent": 10.0,
        "depths": [
            0.1,
            0.25,
            0.5,
            1,
            2,
            5,
            10
        ]
    },

    "IH03": {
        "native": 1.085300,
        "parent": 1.0,
        "depths": [
            0.1,
            0.25,
            0.5,
            1
        ]
    }
}

# ------------------------------------------------------------
# Parent BAM manifest
# ------------------------------------------------------------

parent_manifest = os.path.join(
    manifest_dir,
    "parent_tasks.tsv"
)

with open(
    parent_manifest,
    "w",
    newline=""
) as f:

    w = csv.writer(
        f,
        delimiter="\t"
    )

    w.writerow([
        "sample",
        "native_depth",
        "parent_depth",
        "source_bam",
        "output_bam",
        "seed"
    ])

    for i, (sample, d) in enumerate(
        samples.items()
    ):

        w.writerow([
            sample,
            d["native"],
            d["parent"],
            os.path.join(
                SHENDURE_DIR,
                f"{sample}.bam"
            ),
            os.path.join(
                root,
                "parents",
                f"{sample}_{d['parent']:g}x.bam"
            ),
            101 + i
        ])

# ------------------------------------------------------------
# Spike manifest: 3 backgrounds x 5 replicates
# ------------------------------------------------------------

spike_manifest = os.path.join(
    manifest_dir,
    "spike_tasks.tsv"
)

with open(
    spike_manifest,
    "w",
    newline=""
) as f:

    w = csv.writer(
        f,
        delimiter="\t"
    )

    w.writerow([
        "sample",
        "replicate",
        "parent_depth",
        "parent_bam",
        "truth",
        "output_bam",
        "seed"
    ])

    for sample, d in samples.items():

        parent = os.path.join(
            root,
            "parents",
            f"{sample}_{d['parent']:g}x.bam"
        )

        for rep in range(1, 6):

            truth = os.path.join(
                root,
                "truth",
                f"rep{rep}.truth.tsv"
            )

            output = os.path.join(
                root,
                "spiked",
                sample,
                f"rep{rep}",
                f"{sample}_rep{rep}_"
                f"{d['parent']:g}x.spiked.bam"
            )

            w.writerow([
                sample,
                rep,
                d["parent"],
                parent,
                truth,
                output,
                1000 + rep
            ])

# ------------------------------------------------------------
# Positive simulated BAM depth manifest
# ------------------------------------------------------------

depth_manifest = os.path.join(
    manifest_dir,
    "depth_tasks.tsv"
)

seed = 2000

with open(
    depth_manifest,
    "w",
    newline=""
) as f:

    w = csv.writer(
        f,
        delimiter="\t"
    )

    w.writerow([
        "sample",
        "replicate",
        "target_depth",
        "source_depth",
        "source_bam",
        "output_bam",
        "seed",
        "truth"
    ])

    for sample, d in samples.items():

        for rep in range(1, 6):

            source = os.path.join(
                root,
                "spiked",
                sample,
                f"rep{rep}",
                f"{sample}_rep{rep}_"
                f"{d['parent']:g}x.spiked.bam"
            )

            truth = os.path.join(
                root,
                "truth",
                f"rep{rep}.truth.tsv"
            )

            for depth in d["depths"]:

                output = os.path.join(
                    root,
                    "depth_series",
                    sample,
                    f"rep{rep}",
                    f"{sample}_rep{rep}_"
                    f"{depth:g}x.spiked.bam"
                )

                w.writerow([
                    sample,
                    rep,
                    depth,
                    d["parent"],
                    source,
                    output,
                    seed,
                    truth
                ])

                seed += 1

# ------------------------------------------------------------
# Negative-control manifest
# ------------------------------------------------------------

negative_manifest = os.path.join(
    manifest_dir,
    "negative_tasks.tsv"
)

seed = 5000

with open(
    negative_manifest,
    "w",
    newline=""
) as f:

    w = csv.writer(
        f,
        delimiter="\t"
    )

    w.writerow([
        "sample",
        "replicate",
        "target_depth",
        "source_depth",
        "source_bam",
        "output_bam",
        "seed",
        "truth"
    ])

    for sample, d in samples.items():

        source = os.path.join(
            root,
            "parents",
            f"{sample}_{d['parent']:g}x.bam"
        )

        for depth in d["depths"]:

            output = os.path.join(
                root,
                "negative_controls",
                sample,
                f"{sample}_{depth:g}x.noCNV.bam"
            )

            w.writerow([
                sample,
                0,
                depth,
                d["parent"],
                source,
                output,
                seed,
                "NA"
            ])

            seed += 1

# ------------------------------------------------------------
# Master analysis manifest
# ------------------------------------------------------------

analysis_manifest = os.path.join(
    manifest_dir,
    "analysis_manifest.tsv"
)

with open(
    analysis_manifest,
    "w",
    newline=""
) as f:

    w = csv.writer(
        f,
        delimiter="\t"
    )

    w.writerow([
        "sample_id",
        "background",
        "class",
        "replicate",
        "depth",
        "bam",
        "truth"
    ])

    for sample, d in samples.items():

        # Positive simulations
        for rep in range(1, 6):

            truth = os.path.join(
                root,
                "truth",
                f"rep{rep}.truth.tsv"
            )

            for depth in d["depths"]:

                bam = os.path.join(
                    root,
                    "depth_series",
                    sample,
                    f"rep{rep}",
                    f"{sample}_rep{rep}_"
                    f"{depth:g}x.spiked.bam"
                )

                w.writerow([
                    f"{sample}_rep{rep}_{depth:g}x",
                    sample,
                    "spiked",
                    rep,
                    depth,
                    bam,
                    truth
                ])

        # Negative controls
        for depth in d["depths"]:

            bam = os.path.join(
                root,
                "negative_controls",
                sample,
                f"{sample}_{depth:g}x.noCNV.bam"
            )

            w.writerow([
                f"{sample}_negative_{depth:g}x",
                sample,
                "negative",
                0,
                depth,
                bam,
                "NA"
            ])

print("Truth files and manifests created successfully.")

PY

echo
echo "======================================================"
echo "PREPARATION COMPLETE"
echo "======================================================"

echo
echo "Truth replicates:"
for REP in 1 2 3 4 5
do
    N=$(grep -vc '^#' "$SIM1/truth/rep${REP}.truth.tsv")
    echo "rep${REP}: $N CNVs"
done

echo
echo "Manifest sizes:"
echo "Parents:       $(($(wc -l < "$MAN/parent_tasks.tsv")-1))"
echo "Spiked:        $(($(wc -l < "$MAN/spike_tasks.tsv")-1))"
echo "Positive BAMs: $(($(wc -l < "$MAN/depth_tasks.tsv")-1))"
echo "Negative BAMs: $(($(wc -l < "$MAN/negative_tasks.tsv")-1))"
echo "Total BAMs:    $(($(wc -l < "$MAN/analysis_manifest.tsv")-1))"

echo
echo "Expected:"
echo "Parents        = 3"
echo "Spiked parents = 15"
echo "Positive BAMs  = 95"
echo "Negative BAMs  = 19"
echo "Analysis BAMs  = 114"
echo
