#!/usr/bin/env python
"""A  10-kb bin annotation of each reference: GC fraction, N fraction, ENCODE blacklist overlap.

Output: $MLO/reference/{GRCh37,GRCh38}.bins10kb.tsv.gz   columns chrom start end gc n_frac blacklist_frac
Chromosome names are written without 'chr' (1..22) for both builds. Autosomes only.
The FASTA files are only read (their .fai index must already exist; this script never creates one).
"""
import os, re, sys
from pathlib import Path
import numpy as np, pandas as pd, pysam

BIN = 10000
E = os.environ
OUT = Path(E["MLO"]) / "reference"; OUT.mkdir(parents=True, exist_ok=True)
AUTO = [str(i) for i in range(1, 23)]
nc = lambda c: re.sub(r"^chr", "", str(c))


def blacklist(path):
    bl = {}
    if path and Path(path).exists():
        b = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], dtype={0: str})
        for c, s, e in zip(b[0].map(nc), b[1], b[2]):
            bl.setdefault(c, []).append((int(s), int(e)))
    return bl


for build, fa_path, bl_path in (("GRCh37", E["REF37"], None), ("GRCh38", E["REF38"], E.get("BLACKLIST38"))):
    out = OUT / f"{build}.bins10kb.tsv.gz"
    if out.exists():
        print("SKIP (exists)", out); continue
    if not Path(fa_path + ".fai").exists():
        sys.exit(f"{fa_path}.fai not found - the FASTA index is needed (this script will not create it)")
    fa = pysam.FastaFile(fa_path)
    bl = blacklist(bl_path)
    parts = []
    for contig in fa.references:
        c = nc(contig)
        if c not in AUTO:
            continue
        seq = np.frombuffer(fa.fetch(contig).upper().encode(), dtype=np.uint8)
        L = len(seq); nb = (L + BIN - 1) // BIN
        m = np.full(nb * BIN, 78, dtype=np.uint8); m[:L] = seq; m = m.reshape(nb, BIN); del seq   # pad with 'N'
        gcn = ((m == 67) | (m == 71)).sum(axis=1)                                           # C, G
        acgt = gcn + ((m == 65) | (m == 84)).sum(axis=1)                                    # + A, T
        del m
        start = np.arange(nb) * BIN; end = np.minimum(start + BIN, L); ln = end - start
        blf = np.zeros(nb)
        for s, e in bl.get(c, []):
            for b in range(s // BIN, min((e - 1) // BIN, nb - 1) + 1):
                blf[b] += max(0, min(e, end[b]) - max(s, start[b]))
        with np.errstate(invalid="ignore", divide="ignore"):
            gc = np.where(acgt > 0, gcn / np.maximum(acgt, 1), np.nan)
        parts.append(pd.DataFrame({"chrom": c, "start": start, "end": end, "gc": np.round(gc, 4),
                                   "n_frac": np.round(1 - acgt / ln, 4), "blacklist_frac": np.round(np.minimum(blf / ln, 1), 4)}))
        print(build, c, nb, "bins", flush=True)
    g = pd.concat(parts, ignore_index=True)
    g["chrom"] = pd.Categorical(g["chrom"], AUTO, ordered=True)
    g = g.sort_values(["chrom", "start"]).astype({"chrom": str})
    g.to_csv(out, sep="\t", index=False, compression="gzip")
    print("wrote", out, len(g), "bins")
