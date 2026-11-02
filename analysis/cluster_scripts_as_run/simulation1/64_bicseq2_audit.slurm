#!/bin/bash
#SBATCH --job-name=bic_audit
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/bic_audit_%j.out
#SBATCH --error=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/logs/simulation1/bic_audit_%j.err

set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project

source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

SIM1_OUT="${SIM1_OUT:-$PROJECT/07_tool_outputs/simulation1}"

BICOUT="$SIM1_OUT/BICseq2"
TASKS="$SIM1_MAN/bicseq2_tasks.tsv"

export BICOUT TASKS

python3 <<'PY'
import csv
import os
from pathlib import Path
from collections import Counter

root = Path(os.environ["BICOUT"])
tasks = Path(os.environ["TASKS"])

rows = []

with tasks.open() as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    for r in reader:

        sample = r["sample_id"]
        bin_bp = r["bin_bp"]

        od = root / sample / f"bin_{bin_bp}"
        cnv = od / "seg" / f"{sample}.cnv"
        done = od / ".done"

        if done.exists() and cnv.exists():
            status = "OK"
        else:
            status = "MISSING"

        rows.append({
            **r,
            "status": status,
            "cnv_file": str(cnv),
        })

out = root / "bicseq2_run_status.tsv"

with out.open("w", newline="") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=list(rows[0].keys()),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)

counts = Counter(r["status"] for r in rows)

unique_done = len(
    list((root / "unique").glob("*/.done"))
)

print("======================================")
print("BIC-seq2 Simulation-1 audit")
print("======================================")
print(f"Expected configurations : {len(rows)}")
print(f"OK                      : {counts['OK']}")
print(f"MISSING                 : {counts['MISSING']}")
print(f"Unique samples complete : {unique_done} / 114")
print(f"Status file             : {out}")

if (
    len(rows) == 684
    and counts["OK"] == 684
    and unique_done == 114
):
    print()
    print("ALL BICSEQ2 OUTPUTS COMPLETE")
else:
    print()
    print("BICSEQ2 OUTPUTS INCOMPLETE")
PY
