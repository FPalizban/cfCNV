#!/usr/bin/env python3
"""Simulation 1 truth set: 30 clonal single-copy events per replicate (50 kb - 50 Mb), autosomes only,
>= 1 Mb apart, >= 5 Mb from chromosome ends, alternating loss/gain. Identical to the design used in
the cfCNV paper: with the GRCh37 autosome sizes and seed 1000 + replicate it reproduces
assets/truth/rep{1..5}.truth.tsv exactly.

usage: make_truth.py --chrom-sizes sizes.tsv --replicate 1 [--seed-base 1000] --output rep1.truth.tsv
"""
import argparse, random

DESIGN = {50_000: 4, 100_000: 4, 250_000: 4, 500_000: 4, 1_000_000: 4,
          5_000_000: 4, 10_000_000: 2, 25_000_000: 2, 50_000_000: 2}

p = argparse.ArgumentParser()
p.add_argument("--chrom-sizes", required=True, help="chrom<TAB>length (a .fai works)")
p.add_argument("--replicate", type=int, required=True)
p.add_argument("--seed-base", type=int, default=1000)
p.add_argument("--output", required=True)
a = p.parse_args()

chroms = []
for ln in open(a.chrom_sizes):
    c, L = ln.split()[:2]
    k = c[3:] if c.startswith("chr") else c
    if k.isdigit() and 1 <= int(k) <= 22:
        chroms.append((c, int(L)))
chroms.sort(key=lambda x: int(x[0].replace("chr", "")))

random.seed(a.seed_base + a.replicate)
events, n = [], 1

def overlaps(c, s, e, pad=1_000_000):
    return any(ec == c and s < ee + pad and e > es - pad for ec, es, ee, *_ in events)

for size, k in DESIGN.items():
    for _ in range(k):
        for _attempt in range(10000):
            c, L = random.choice(chroms)
            margin = max(5_000_000, int(size * 0.2))
            if L <= size + 2 * margin:
                continue
            s = random.randint(margin, L - size - margin); e = s + size
            if overlaps(c, s, e):
                continue
            t = "loss" if n % 2 else "gain"
            sc = "small" if size <= 250_000 else "medium" if size <= 5_000_000 else "large"
            events.append((c, s, e, t, 1 if t == "loss" else 3, sc, f"REP{a.replicate}_CNV{n}", size))
            n += 1
            break
        else:
            raise RuntimeError(f"could not place an event of {size} bp")

events.sort(key=lambda x: (int(x[0].replace("chr", "")), x[1]))
with open(a.output, "w") as out:
    out.write("#chrom\tstart\tend\tcnv_type\tcopy_number\tsize_class\tcnv_id\tcnv_size\n")
    for ev in events:
        out.write("\t".join(map(str, ev)) + "\n")
