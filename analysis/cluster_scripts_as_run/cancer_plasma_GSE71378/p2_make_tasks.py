#!/usr/bin/env python3
"""GSE71378: analysis BAM manifest + tool task lists (columns: sample_id bam depth param).

Only the configurations that were best in Simulations 1 and 2 are run on the real data:
QDNAseq 30/50/100 kb (+500 kb, 1 Mb for ACE), CNVkit 30/100 kb, BIC-seq2 10/30 kb,
Spectre 100/500 kb, ichorCNA (1 Mb low-TF, 1 Mb tumour defaults, 500 kb), ACE, Fragle.
Also writes p3_readlen.txt (most common read length) for BIC-seq2 and p3_qc.tsv."""
import csv, os, re
from collections import Counter
from pathlib import Path

M = Path(os.environ["S2_META"]); B = Path(os.environ["P3_BAM"])
rows, qc, rl = [], [], Counter()
for d in sorted(B.glob("*/derived.tsv")):
    S = d.parent.name
    for line in d.read_text().splitlines():
        sid, base, bam, dep, level, layout = line.split("\t")
        rows.append({"sample_id": sid, "base_sample": base, "bam": bam, "depth_x": dep, "level": level, "layout": layout})
    st = (d.parent / f"{S}.stats.txt").read_text()
    sn = dict(re.findall(r"^SN\t([^:]+):\t(\S+)", st, re.M))
    rlh = re.findall(r"^RL\t(\d+)\t(\d+)", st, re.M)
    for L, n in rlh: rl[int(L)] += int(n)
    mk = (d.parent / f"{S}.markdup.txt"); dup = ""
    if mk.exists():
        t = dict(re.findall(r"^([A-Z ]+):\s*(\d+)", mk.read_text(), re.M))
        if t.get("READ") and t.get("DUPLICATE TOTAL"):
            dup = round(int(t["DUPLICATE TOTAL"]) / int(t["READ"]), 4)
    qc.append({"sample": S, "reads_after_dedup": sn.get("raw total sequences", ""), "mapped": sn.get("reads mapped", ""),
               "insert_size_mean": sn.get("insert size average", ""), "insert_size_sd": sn.get("insert size standard deviation", ""),
               "dup_rate": dup, "depth_x": [r["depth_x"] for r in rows if r["sample_id"] == f"{S}_full"][0]})
with open(M / "p3_bams.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
with open(M / "p3_qc.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(qc[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(qc)
readlen = rl.most_common(1)[0][0] if rl else 150
(M / "p3_readlen.txt").write_text(f"{readlen}\n")

T = {
    "qdnaseq":    [(r, b) for r in rows for b in (30000, 50000, 100000, 500000, 1000000)],
    "cnvkit":     [(r, b) for r in rows for b in (30000, 100000)],
    "ichor":      [(r, c) for r in rows for c in ("1000kb_lowTF", "1000kb_default", "500kb_default")],
    "spectre":    [(r, b) for r in rows for b in (100000, 500000)],
    "ace":        [(r, b) for r in rows for b in (100000, 500000, 1000000)],
    "bic_unique": [(r, "all") for r in rows],
    "bic_norm":   [(r, b) for r in rows for b in (10000, 30000)],
    "fragle":     [(r, "R") for r in rows if r["layout"] == "PAIRED"],
}
for k, tt in T.items():
    with open(M / f"tasks_{k}.tsv", "w") as f:
        f.write("sample_id\tbam\tdepth\tparam\n")
        for r, p in tt:
            f.write(f"{r['sample_id']}\t{r['bam']}\t{r['depth_x']}\t{p}\n")
    print(f"  {k:11s} {len(tt):5d} tasks")
print(f"{len(rows)} analysis BAMs from {len(qc)} samples; read length {readlen} -> BIC-seq2 resources {readlen}mer")
