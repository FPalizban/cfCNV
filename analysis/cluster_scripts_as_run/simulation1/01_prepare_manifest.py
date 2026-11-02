#!/usr/bin/env python3

import csv
import os
import re
from pathlib import Path

sim = Path(os.environ["SIM1_ROOT"])
src = sim / "manifests" / "analysis_manifest.tsv"
out = Path(os.environ["SIM1_MAN"])
out.mkdir(parents=True, exist_ok=True)

rows = []

def val(row, *names):
    low = {k.lower(): k for k in row}
    for n in names:
        if n.lower() in low:
            x = row[low[n.lower()]]
            if x and x.strip():
                return x.strip()
    return ""

def infer(text, pattern):
    m = re.search(pattern, text, re.I)
    return m.group(1) if m else ""

with open(src) as f:
    rdr = csv.DictReader(f, delimiter="\t")

    for r in rdr:

        bam = val(
            r,
            "bam",
            "bam_path",
            "analysis_bam",
            "path",
            "file"
        )

        if not bam:
            for x in r.values():
                if x and x.strip().endswith(".bam"):
                    bam = x.strip()
                    break

        if not bam:
            continue

        if not Path(bam).exists():
            continue

        sid = (
            val(r, "sample_id", "sample", "id", "name")
            or Path(bam).stem
        )

        text = " ".join(
            [sid, bam] +
            [str(x) for x in r.values() if x is not None]
        )

        bg = (
            val(r, "background")
            or infer(text, r"\b(IH0[123])\b")
        )

        dep = val(
            r,
            "target_depth",
            "depth",
            "coverage"
        )

        if dep:
            m = re.search(r"([0-9]+(?:\.[0-9]+)?)", dep)
            if m:
                dep = m.group(1)

        if not dep:
            dep = infer(
                text,
                r"(0[._p]1|0[._p]25|0[._p]5|1|2|5|10|20)[xX]"
            )

            dep = (
                dep
                .replace("_", ".")
                .replace("p", ".")
            )

        rep = (
            val(r, "replicate", "rep", "truth_rep")
            or infer(text, r"(?:rep|r)[_-]?([1-5])\b")
        )

        cls = val(
            r,
            "class",
            "status",
            "sample_class"
        )

        if not cls:
            if "negative" in text.lower():
                cls = "negative"
            else:
                cls = "positive"

        rows.append({
            "sample_id": sid,
            "bam": str(Path(bam).resolve()),
            "background": bg,
            "target_depth": dep,
            "replicate": rep,
            "class": cls.lower()
        })


if len(rows) != 114:
    raise SystemExit(
        f"ERROR: parsed {len(rows)} BAMs; expected 114. "
        f"Inspect {src}"
    )


fields = [
    "sample_id",
    "bam",
    "background",
    "target_depth",
    "replicate",
    "class"
]


# -----------------------------
# Master manifest
# -----------------------------

with open(
    out / "master.tsv",
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    w.writeheader()
    w.writerows(rows)


# -----------------------------
# Caller task tables
# -----------------------------

bins = [
    10000,
    30000,
    50000,
    100000,
    500000,
    1000000
]


for name, blist in [
    ("qdnaseq", bins),
    ("cnvkit", bins),
    ("ichor", [1000000])
]:

    tasks = []

    for r in rows:
        for b in blist:
            x = dict(r)
            x["bin_size"] = b
            tasks.append(x)

    with open(
        out / f"{name}_tasks.tsv",
        "w",
        newline=""
    ) as f:

        w = csv.DictWriter(
            f,
            fieldnames=fields + ["bin_size"],
            delimiter="\t",
            lineterminator="\n"
        )

        w.writeheader()
        w.writerows(tasks)


# -----------------------------
# Spectre
# -----------------------------

with open(
    out / "spectre_tasks.tsv",
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    w.writeheader()
    w.writerows(rows)


print("====================================")
print("Simulation 1 manifests created")
print("====================================")
print("Samples:   ", len(rows))
print("QDNAseq:   ", len(rows) * 6)
print("CNVkit:    ", len(rows) * 6)
print("ichorCNA:  ", len(rows))
print("Spectre:   ", len(rows))
print("====================================")
