#!/usr/bin/env python3
"""Reproducibility check: compare the workflow's harmonised calls with the calls reported in the paper.

usage: compare_with_paper.py <new all_calls.tsv> <paper calls.tsv(.gz)> --map NEW_SAMPLE=PAPER_SAMPLE [...]
example (GSE71378 0.1x BAMs given to the workflow as samples IC37 and IH03):
  python analysis/compare_with_paper.py results/summary/all_calls.tsv \
      analysis/real_plasma/paper_calls_GSE71378_primary.tsv.gz --map IC37=IC37_d0p1 IH03=IH03_d0p1
Prints, per sample and caller, the number of calls in each run and the base-pair Jaccard index of the altered
genome (gains and losses separately, then summed). 1.00 = identical calls.
"""
import argparse
import pandas as pd


def merged(rows):
    out = {}
    for c, s, e in sorted(rows):
        L = out.setdefault(c, [])
        if L and s <= L[-1][1]: L[-1][1] = max(L[-1][1], e)
        else: L.append([s, e])
    return out


def length(m): return sum(e - s for L in m.values() for s, e in L)


def inter(A, B):
    t = 0
    for c in A:
        a, b, i, j = A[c], B.get(c, []), 0, 0
        while i < len(a) and j < len(b):
            lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
            if hi > lo: t += hi - lo
            if a[i][1] < b[j][1]: i += 1
            else: j += 1
    return t


ap = argparse.ArgumentParser(); ap.add_argument("new"); ap.add_argument("paper"); ap.add_argument("--map", nargs="+", required=True)
a = ap.parse_args()
N = pd.read_csv(a.new, sep="\t", dtype={"chrom": str}); P = pd.read_csv(a.paper, sep="\t", dtype={"chrom": str})
for d in (N, P): d["chrom"] = d.chrom.str.replace("^chr", "", regex=True)
pairs = set(zip(P.tool, P.config))                      # compare the same caller settings only
N = N[[(t, c) in pairs for t, c in zip(N.tool, N.config)]]
rows = []
for pair in a.map:
    ns, ps = pair.split("=")
    for tool in sorted(set(N.tool) | set(P.tool)):
        n = N[(N.sample_id == ns) & (N.tool == tool)]; p = P[(P.sample_id == ps) & (P.tool == tool)]
        if n.empty and p.empty: continue
        I = U = 0
        for t in ("gain", "loss"):
            A = merged(zip(n[n.type == t].chrom, n[n.type == t].start, n[n.type == t].end))
            B = merged(zip(p[p.type == t].chrom, p[p.type == t].start, p[p.type == t].end))
            i = inter(A, B); I += i; U += length(A) + length(B) - i
        rows.append({"sample": ns, "paper_sample": ps, "tool": tool, "calls_workflow": len(n), "calls_paper": len(p),
                     "bp_jaccard": round(I / U, 3) if U else 1.0})
print(pd.DataFrame(rows).to_string(index=False))
