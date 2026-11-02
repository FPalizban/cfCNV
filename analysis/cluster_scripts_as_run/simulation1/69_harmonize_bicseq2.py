#!/usr/bin/env python3

import csv
import os
import re
from pathlib import Path

PROJECT = Path(
    "/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project"
)

sim1_man = Path(os.environ["SIM1_MAN"])

sim1_out = Path(
    os.environ.get(
        "SIM1_OUT",
        str(PROJECT / "07_tool_outputs/simulation1")
    )
)

tasks_file = sim1_man / "bicseq2_tasks.tsv"

bic_root = sim1_out / "BICseq2"

outdir = PROJECT / "08_harmonized_calls/simulation1/BICseq2"
outdir.mkdir(parents=True, exist_ok=True)

segments_out = outdir / "BICseq2_all_segments.tsv"
calls_out = outdir / "BICseq2_calls_log2_0.2.tsv"
status_out = outdir / "BICseq2_harmonization_status.tsv"


def norm(s):
    return re.sub(r"[^a-z0-9]", "", s.lower())


def find_col(header, candidates):
    h = {norm(x): i for i, x in enumerate(header)}

    for c in candidates:
        if norm(c) in h:
            return h[norm(c)]

    return None


with tasks_file.open() as fh:
    tasks = list(csv.DictReader(fh, delimiter="\t"))


segment_rows = []
call_rows = []
status_rows = []


for task in tasks:

    sample = task["sample_id"]
    bin_bp = task["bin_bp"]

    cnv = (
        bic_root
        / sample
        / f"bin_{bin_bp}"
        / "seg"
        / f"{sample}.cnv"
    )

    if not cnv.exists():

        status_rows.append({
            "task_id": task["task_id"],
            "sample_id": sample,
            "bin_bp": bin_bp,
            "status": "MISSING",
            "message": str(cnv),
        })

        continue

    with cnv.open() as fh:

        lines = [
            line.rstrip("\n")
            for line in fh
            if line.strip()
        ]

    if not lines:

        status_rows.append({
            "task_id": task["task_id"],
            "sample_id": sample,
            "bin_bp": bin_bp,
            "status": "EMPTY",
            "message": str(cnv),
        })

        continue

    header = re.split(r"\s+", lines[0].strip())

    chr_i = find_col(
        header,
        ["chrom", "chromosome", "chromName"]
    )

    start_i = find_col(
        header,
        ["start", "startPos"]
    )

    end_i = find_col(
        header,
        ["end", "endPos"]
    )

    log2_i = find_col(
        header,
        [
            "log2.copyRatio",
            "log2copyRatio",
            "log2ratio",
            "log2"
        ]
    )

    p_i = find_col(
        header,
        [
            "pvalue",
            "p.value",
            "pval"
        ]
    )

    if None in (chr_i, start_i, end_i, log2_i):

        raise RuntimeError(
            "\nCould not identify BIC-seq2 columns.\n"
            f"File: {cnv}\n"
            f"Header: {header}\n"
        )

    nseg = 0
    ncalls = 0

    for line in lines[1:]:

        fields = re.split(r"\s+", line.strip())

        if len(fields) <= max(chr_i, start_i, end_i, log2_i):
            continue

        try:
            chrom = fields[chr_i]
            start1 = int(float(fields[start_i]))
            end1 = int(float(fields[end_i]))
            log2ratio = float(fields[log2_i])
        except ValueError:
            continue

        # BIC-seq2 coordinates are represented as genomic
        # segment coordinates. Store both original coordinates
        # and BED-compatible 0-based start.
        start0 = max(0, start1 - 1)

        pvalue = ""

        if p_i is not None and p_i < len(fields):
            pvalue = fields[p_i]

        if log2ratio >= 0.2:
            state = "GAIN"
        elif log2ratio <= -0.2:
            state = "LOSS"
        else:
            state = "NEUTRAL"

        row = {
            "tool": "BICseq2",
            "task_id": task["task_id"],
            "sample_id": sample,
            "background": task["background"],
            "target_depth": task["target_depth"],
            "replicate": task["replicate"],
            "class": task["class"],
            "bin_bp": bin_bp,
            "chrom": chrom,
            "start_1based": start1,
            "start_0based": start0,
            "end": end1,
            "log2_ratio": log2ratio,
            "pvalue": pvalue,
            "state_log2_0.2": state,
            "source_file": str(cnv),
        }

        segment_rows.append(row)
        nseg += 1

        if state != "NEUTRAL":
            call_rows.append(row)
            ncalls += 1

    status_rows.append({
        "task_id": task["task_id"],
        "sample_id": sample,
        "bin_bp": bin_bp,
        "status": "OK",
        "message": f"segments={nseg};calls={ncalls}",
    })


fields = [
    "tool",
    "task_id",
    "sample_id",
    "background",
    "target_depth",
    "replicate",
    "class",
    "bin_bp",
    "chrom",
    "start_1based",
    "start_0based",
    "end",
    "log2_ratio",
    "pvalue",
    "state_log2_0.2",
    "source_file",
]


with segments_out.open("w", newline="") as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writeheader()
    w.writerows(segment_rows)


with calls_out.open("w", newline="") as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writeheader()
    w.writerows(call_rows)


with status_out.open("w", newline="") as fh:

    fields_status = [
        "task_id",
        "sample_id",
        "bin_bp",
        "status",
        "message",
    ]

    w = csv.DictWriter(
        fh,
        fieldnames=fields_status,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writeheader()
    w.writerows(status_rows)


print("====================================")
print("BIC-seq2 harmonization complete")
print("====================================")
print(f"Tasks           : {len(tasks)}")
print(f"Segments        : {len(segment_rows)}")
print(f"GAIN/LOSS calls : {len(call_rows)}")
print()
print(segments_out)
print(calls_out)
print(status_out)
