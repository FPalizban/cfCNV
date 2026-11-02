#!/usr/bin/env python3

import csv
import os
from pathlib import Path

PROJECT = Path(
    "/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project"
)

sim1_man_env = os.environ.get("SIM1_MAN", "").strip()

if not sim1_man_env:
    raise SystemExit(
        "ERROR: SIM1_MAN is not defined. "
        "Source config.sh and 00_env.sh before running this script."
    )

MAN = Path(sim1_man_env)

master = MAN / "master.tsv"
out = MAN / "bicseq2_tasks.tsv"

if not master.exists():
    raise SystemExit(
        f"ERROR: master.tsv not found at:\n{master}"
    )

bins = [
    10_000,
    30_000,
    50_000,
    100_000,
    500_000,
    1_000_000,
]

with master.open() as fh:
    samples = list(csv.DictReader(fh, delimiter="\t"))

print(f"Found {len(samples)} samples")
print(f"Master manifest: {master}")

if len(samples) != 114:
    raise SystemExit(
        f"ERROR: expected 114 Simulation-1 samples, found {len(samples)}"
    )

required = [
    "sample_id",
    "bam",
    "background",
    "target_depth",
    "replicate",
    "class",
]

missing = [
    x for x in required
    if x not in samples[0]
]

if missing:
    raise SystemExit(
        f"ERROR: master.tsv missing columns: {missing}"
    )

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

rows = []

task_id = 0

for s in samples:
    for bin_bp in bins:

        rows.append({
            "task_id": task_id,
            "sample_id": s["sample_id"],
            "bam": s["bam"],
            "background": s["background"],
            "target_depth": s["target_depth"],
            "replicate": s["replicate"],
            "class": s["class"],
            "bin_bp": bin_bp,
        })

        task_id += 1

if len(rows) != 684:
    raise SystemExit(
        f"ERROR: expected 684 BIC-seq2 tasks, generated {len(rows)}"
    )

with out.open("w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(rows)

print()
print("BIC-seq2 task manifest created successfully")
print(f"Samples : {len(samples)}")
print(f"Bins    : {len(bins)}")
print(f"Tasks   : {len(rows)}")
print(f"Output  : {out}")
