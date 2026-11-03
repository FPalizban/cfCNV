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
ap.add_argument("--tf", help="workflow summary/tumour_fraction.tsv"); ap.add_argument("--paper-tf", help="paper tumour-fraction table")
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

# tumour fraction: same method and setting in the workflow and in the paper
METHOD = {"ichorCNA 1Mb low-TF": "ichorCNA 1000kb_lowTF", "ACE ploidy 2 100kb": "ACE ploidy2 100kb", "Fragle raw WGS": "Fragle"}
ok = all(r["bp_jaccard"] >= 0.95 for r in rows)
if a.tf and a.paper_tf:
    T = pd.read_csv(a.tf, sep="\t"); PT = pd.read_csv(a.paper_tf, sep="\t"); tr = []
    for pair in a.map:
        ns, ps = pair.split("=")
        for m, pm in METHOD.items():
            x = T[(T.sample_id == ns) & (T.method == m)]; y = PT[(PT.sample_id == ps) & (PT.method == pm)]
            if len(x) and len(y):
                d = abs(float(x.TF.iloc[0]) - float(y.TF.iloc[0]))
                tr.append({"sample": ns, "method": m, "TF_workflow": round(float(x.TF.iloc[0]), 4),
                           "TF_paper": round(float(y.TF.iloc[0]), 4), "abs_difference": round(d, 4)})
                ok = ok and d <= 0.01
    print("\nTumour fraction"); print(pd.DataFrame(tr).to_string(index=False))
print("\nRESULT:", "SAME AS THE PAPER (calls bp_jaccard >= 0.95, tumour fraction within 0.01)" if ok
      else "DIFFERENCES FOUND - send this table")
