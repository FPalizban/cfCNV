#!/usr/bin/env python3

import csv
import os
from pathlib import Path

PROJECT = Path(
    "/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project"
)

SIM = Path(os.environ["SIM"])
SIM1_MAN = Path(os.environ["SIM1_MAN"])

DUDE = PROJECT / "07_tool_outputs/simulation1/dudeML"

META = DUDE / "metadata"
LABEL = META / "labels"

META.mkdir(parents=True, exist_ok=True)
LABEL.mkdir(parents=True, exist_ok=True)

MASTER = SIM1_MAN / "master.tsv"

BINS = [
    10000,
    30000,
    50000,
    100000,
    500000,
    1000000,
]


# ------------------------------------------------------------
# CHROMOSOME LIST
# ------------------------------------------------------------

with (META / "chromosomes.txt").open("w") as out:
    for c in range(1, 23):
        out.write(f"{c}\n")


# ------------------------------------------------------------
# CONVERT SIMULATION TRUTH → DUDEML FORMAT
#
# gain:
#     dup 0.5 2
#     expected relative copy = 1.5
#
# loss:
#     del 0.5 0
#     expected relative copy = 0.5
# ------------------------------------------------------------

truth_dir = (
    SIM
    / "healthy_cfDNA_spikein"
    / "truth"
)

for rep in range(1, 6):

    infile = truth_dir / f"rep{rep}.truth.tsv"

    dupfile = LABEL / f"rep{rep}.dup.bed"
    delfile = LABEL / f"rep{rep}.del.bed"

    nd = 0
    nl = 0

    with infile.open() as fh, \
         dupfile.open("w") as dup, \
         delfile.open("w") as dele:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        for row in reader:

            chrom = row["#chrom"].replace("chr", "")
            start = int(row["start"])
            end = int(row["end"])
            typ = row["cnv_type"].lower()

            if typ == "gain":

                dup.write(
                    f"{chrom}\t{start}\t{end}"
                    f"\tdup\t0.5\t2\n"
                )
                nd += 1

            elif typ == "loss":

                dele.write(
                    f"{chrom}\t{start}\t{end}"
                    f"\tdel\t0.5\t0\n"
                )
                nl += 1

            else:
                raise RuntimeError(
                    f"Unknown CNV type: {typ}"
                )

    if nd != 15 or nl != 15:
        raise RuntimeError(
            f"rep{rep}: gain={nd}, loss={nl}; "
            "expected 15/15"
        )

    print(
        f"rep{rep}: "
        f"{nd} gains + {nl} losses"
    )


# ------------------------------------------------------------
# READ MASTER MANIFEST
# ------------------------------------------------------------

with MASTER.open() as fh:

    samples = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )

if len(samples) != 114:
    raise RuntimeError(
        f"Expected 114 samples; found {len(samples)}"
    )


# ------------------------------------------------------------
# PREDICTION TASKS
#
# 114 samples × 6 bins = 684
# ------------------------------------------------------------

pred_file = META / "prediction_tasks.tsv"

with pred_file.open("w", newline="") as fh:

    fields = [
        "task_id",
        "sample_id",
        "bam",
        "background",
        "target_depth",
        "replicate",
        "class",
        "bin_bp",
    ]

    w = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writeheader()

    task = 0

    for s in samples:

        for b in BINS:

            w.writerow({
                "task_id": task,
                "sample_id": s["sample_id"],
                "bam": s["bam"],
                "background": s["background"],
                "target_depth": s["target_depth"],
                "replicate": s["replicate"],
                "class": s["class"],
                "bin_bp": b,
            })

            task += 1


# ------------------------------------------------------------
# MODEL TASKS
#
# One model for each:
#
# depth × bin × held-out replicate
#
# 8 depths × 6 bins × 5 folds = 240
# ------------------------------------------------------------

positive_depths = sorted(
    {
        s["target_depth"]
        for s in samples
        if s["class"].lower() == "positive"
    },
    key=float
)

model_file = META / "model_tasks.tsv"

with model_file.open("w", newline="") as fh:

    fields = [
        "task_id",
        "target_depth",
        "bin_bp",
        "heldout_replicate",
    ]

    w = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writeheader()

    task = 0

    for depth in positive_depths:

        for b in BINS:

            for fold in range(1, 6):

                w.writerow({
                    "task_id": task,
                    "target_depth": depth,
                    "bin_bp": b,
                    "heldout_replicate": fold,
                })

                task += 1


print()
print("====================================")
print("dudeML metadata ready")
print("====================================")
print("Samples          :", len(samples))
print("Prediction tasks :", 114 * 6)
print("Model tasks      :", len(positive_depths) * 6 * 5)
print("Depths           :", positive_depths)
print("Metadata         :", META)
