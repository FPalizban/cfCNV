#!/usr/bin/env python3
"""Simulation 2: task lists for every tool, from the mixtures that were actually built
(mixtures_v2/*.done). Columns of every task file: sample_id  bam  depth  param."""
import csv, os
from pathlib import Path

M = Path(os.environ["S2_META"]); MIX = Path(os.environ["S2_MIX"])
man = {r["sample_id"]: r for r in csv.DictReader(open(M / "HG008_mixtures_v2.tsv"), delimiter="\t")}
rows = []
for d in sorted(MIX.glob("*.done")):
    lines = d.read_text().strip().splitlines()
    rec = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
    s = rec["sample_id"]; bam = MIX / f"{s}.bam"
    if s in man and bam.exists() and Path(str(bam) + ".bai").exists():
        rows.append({**man[s], "bam": str(bam), "realised_fraction": rec["realised_fraction"]})
with open(M / "HG008_mixtures_realised.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(rows)

B = [10000, 30000, 50000, 100000, 500000, 1000000]
T = {
    "qdnaseq":     [(r, b) for r in rows for b in B],
    "cnvkit":      [(r, b) for r in rows for b in B],
    "ichor":       [(r, c) for r in rows for c in ("1000kb_default", "1000kb_lowTF", "500kb_default")],
    "spectre":     [(r, b) for r in rows for b in B],
    "ace":         [(r, b) for r in rows for b in (100000, 500000, 1000000)],
    "dudeml_feat": [(r, "all") for r in rows],
    "dudeml_pred": [(r, b) for r in rows for b in B],
    "meann_feat":  [(r, b) for r in rows for b in B],
    "bic_unique":  [(r, "all") for r in rows],
    "bic_norm":    [(r, b) for r in rows for b in B],
    "lyceum":      [(r, "exon") for r in rows],
}
for k, tt in T.items():
    with open(M / f"tasks_{k}.tsv", "w") as f:
        f.write("sample_id\tbam\tdepth\tparam\n")
        for r, p in tt:
            f.write(f"{r['sample_id']}\t{r['bam']}\t{r['depth_x']}\t{p}\n")
    print(f"  {k:12s} {len(tt):5d} tasks")
print(len(rows), "mixtures ->", M / "HG008_mixtures_realised.tsv")
