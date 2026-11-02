#!/usr/bin/env python3
"""Score CNV calls against a truth set, with the matching rule of the cfCNV paper.

Event level: one-to-one matching, same chromosome, same direction (gain/loss), reciprocal overlap >= 0.5
(greedy by overlap) -> TP, FP, FN, precision, recall, F1; sensitivity also at 10% and 30% overlap.
Base-pair level: precision, recall, specificity, MCC over the autosomes.
Negative controls (no truth file): number of calls and false calls per Gb.
Per truth event: detected yes/no (for size-stratified sensitivity).

usage: score_calls.py --calls all_calls.tsv --truth-map truth_map.tsv --fai ref.fa.fai [--overlap 0.5] [--outdir .]
  all_calls.tsv : sample_id ... tool config chrom start end type   (output of collect_calls.py), or a plain
                  BED-like file 'chrom start end type' together with --sample NAME
  truth_map.tsv : sample_id <TAB> truth_file   (truth_file empty or NA = negative control)
outputs: event_scores.tsv, truth_detection.tsv
"""
import argparse, math
from pathlib import Path
import numpy as np, pandas as pd

AUTO = [str(i) for i in range(1, 23)]
nc = lambda c: str(c)[3:] if str(c).startswith("chr") else str(c)


def read_truth(f):
    d = pd.read_csv(f, sep="\t", dtype=str, comment=None)
    d.columns = [c.lstrip("#") for c in d.columns]
    tcol = "cnv_type" if "cnv_type" in d.columns else "type"
    return [(nc(c), int(s), int(e), t.lower()) for c, s, e, t in zip(d.chrom, d.start, d.end, d[tcol]) if nc(c) in AUTO]


def match(truth, calls, thr):
    pairs = []
    for i, t in enumerate(truth):
        for j, p in enumerate(calls):
            if p[0] != t[0] or p[3] != t[3]:
                continue
            ov = min(t[2], p[2]) - max(t[1], p[1])
            if ov <= 0:
                continue
            ro = min(ov / (t[2] - t[1]), ov / (p[2] - p[1]))
            if ro >= thr:
                pairs.append((ro, i, j))
    ut, up, m = set(), set(), {}
    for ro, i, j in sorted(pairs, reverse=True):
        if i in ut or j in up:
            continue
        ut.add(i); up.add(j); m[i] = j
    return m


def merged(iv):
    out = {}
    for c, s, e, *_ in sorted(iv):
        L = out.setdefault(c, [])
        if L and s <= L[-1][1]:
            L[-1][1] = max(L[-1][1], e)
        else:
            L.append([s, e])
    return out


def inter(A, B):
    tot = 0
    for c in A:
        a, b, i, j = A[c], B.get(c, []), 0, 0
        while i < len(a) and j < len(b):
            lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
            if hi > lo: tot += hi - lo
            if a[i][1] < b[j][1]: i += 1
            else: j += 1
    return tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True); ap.add_argument("--truth-map", required=True)
    ap.add_argument("--fai", required=True); ap.add_argument("--overlap", type=float, default=0.5)
    ap.add_argument("--sample"); ap.add_argument("--outdir", default=".")
    a = ap.parse_args(); O = Path(a.outdir); O.mkdir(parents=True, exist_ok=True)
    G = sum(int(l.split()[1]) for l in open(a.fai) if nc(l.split()[0]) in AUTO)

    C = pd.read_csv(a.calls, sep="\t", dtype={"chrom": str})
    if "sample_id" not in C.columns:  # plain BED-like file
        C = pd.read_csv(a.calls, sep="\t", header=None, comment="#", dtype={0: str}).iloc[:, :4]
        C.columns = ["chrom", "start", "end", "type"]
        C = C[pd.to_numeric(C.start, errors="coerce").notna()].astype({"start": int, "end": int})
        C["sample_id"] = a.sample or "sample"; C["tool"] = "user"; C["config"] = "-"
    C["chrom"] = C.chrom.map(nc); C["type"] = C.type.str.lower().map(
        lambda t: "loss" if t.startswith(("loss", "del")) else "gain" if t.startswith(("gain", "dup", "amp")) else t)
    TM = pd.read_csv(a.truth_map, sep="\t", header=None, names=["sample_id", "truth"], dtype=str).fillna("")
    truth = {s: (read_truth(f) if f not in ("", "NA") else None) for s, f in zip(TM.sample_id, TM.truth)}

    runs = C[["sample_id", "tool", "config"]].drop_duplicates()
    if Path(a.calls).with_name("call_burden.tsv").exists():   # include caller runs that made no calls
        b = pd.read_csv(Path(a.calls).with_name("call_burden.tsv"), sep="\t", dtype=str)
        runs = pd.concat([runs, b[["sample_id", "tool", "config"]]]).drop_duplicates()
    rows, det = [], []
    for r in runs.itertuples():
        if r.sample_id not in truth:
            continue
        g = C[(C.sample_id == r.sample_id) & (C.tool == r.tool) & (C.config == r.config)]
        calls = [(c, s, e, t) for c, s, e, t in zip(g.chrom, g.start, g.end, g.type) if c in AUTO and e > s]
        T = truth[r.sample_id]
        base = {"sample_id": r.sample_id, "tool": r.tool, "config": r.config, "n_calls": len(calls)}
        if T is None:
            rows.append({**base, "class": "negative", "n_truth": 0, "false_calls_per_Gb": len(calls) / (G / 1e9),
                         "false_altered_fraction": sum(e - s for _, s, e, _ in calls) / G})
            continue
        m = match(T, calls, a.overlap); tp = len(m); fp = len(calls) - tp; fn = len(T) - tp
        p = tp / len(calls) if calls else np.nan; rc = tp / len(T) if T else np.nan
        f1 = 2 * p * rc / (p + rc) if calls and p + rc > 0 else 0.0
        sens = {th: len(match(T, calls, th)) / len(T) for th in (0.1, 0.3)}
        tm, cm = {}, {}
        for typ in ("gain", "loss"):
            tm[typ] = merged([x for x in T if x[3] == typ]); cm[typ] = merged([x for x in calls if x[3] == typ])
        tb = sum(e - s for typ in tm for L in tm[typ].values() for s, e in L)
        cb = sum(e - s for typ in cm for L in cm[typ].values() for s, e in L)
        ob = sum(inter(tm[t], cm[t]) for t in tm)
        bfp, bfn = cb - ob, tb - ob; btn = G - tb - bfp
        den = math.sqrt(max((ob + bfp) * (ob + bfn) * (btn + bfp) * (btn + bfn), 1))
        rows.append({**base, "class": "positive", "n_truth": len(T), "TP": tp, "FP": fp, "FN": fn,
                     "precision": p, "recall": rc, "F1": f1, "recall_10pct": sens[0.1], "recall_30pct": sens[0.3],
                     "false_calls_per_Gb": fp / (G / 1e9),
                     "bp_precision": ob / cb if cb else np.nan, "bp_recall": ob / tb if tb else np.nan,
                     "bp_specificity": btn / (btn + bfp) if btn + bfp else np.nan, "bp_MCC": (ob * btn - bfp * bfn) / den})
        for i, t in enumerate(T):
            det.append({"sample_id": r.sample_id, "tool": r.tool, "config": r.config, "chrom": t[0], "start": t[1],
                        "end": t[2], "type": t[3], "size_bp": t[2] - t[1], "detected": int(i in m)})
    pd.DataFrame(rows).to_csv(O / "event_scores.tsv", sep="\t", index=False)
    pd.DataFrame(det).to_csv(O / "truth_detection.tsv", sep="\t", index=False)
    print(f"scored {len(rows)} caller runs")


if __name__ == "__main__":
    main()
