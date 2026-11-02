#!/usr/bin/env python3
"""
cfCNV Simulation 1 - background-artifact analysis (IH01 and the others).

  python 04_background_mask.py --bundle <cfCNV_Sim1_results_for_local> --out mask/

1. Characterise negative-control calls per background: count, size, gain/loss,
   chromosomes, and recurrence across depths (does a call at depth d1 reappear
   at d2?). Recurrent calls = signal intrinsic to the background, not noise.
2. Build a "self-normal" mask per tool x bin x background: regions that the
   SAME tool called in that background's unmodified negative controls at two
   or more depths (reproducible background signal, not one-off noise) - the
   analogue of a per-tool panel of normals from the same donor.
   Two variants: 'all_depths' (all controls) and 'lodo' (leave-one-depth-out:
   the mask for a sample at depth d is built only from controls at other
   depths, so the sample's own-depth control never contributes).
3. Re-score positives with and without the mask, using the bundle's rules
   (one-to-one, same chromosome + direction, reciprocal overlap >= 0.5).
   Calls with > 50% of their length in the mask are removed; truth events are
   all kept. Masking can also remove TRUE calls that fall in masked regions -
   the change in recall measures that cost.
   The unmasked re-score is compared with the bundle's F1 as a validation.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ALLB = {10000, 30000, 50000, 100000, 500000, 1000000}
VALID = {"QDNAseq": ALLB, "CNVkit": ALLB, "ichorCNA": {1000000}, "Spectre": {1000}}


def merge(iv):
    out = []
    for c, s, e in sorted(iv):
        if out and out[-1][0] == c and s <= out[-1][2]:
            out[-1][2] = max(out[-1][2], e)
        else:
            out.append([c, s, e])
    return out


def frac_in(mask, c, s, e):
    L = max(e - s, 1); cov = 0
    for mc, ms, me in mask.get(c, []):
        if me <= s:
            continue
        if ms >= e:
            break
        cov += min(e, me) - max(s, ms)
    return cov / L


def match(truth, calls, thr=0.5):
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
    used_t, used_p, n = set(), set(), 0
    for ro, i, j in sorted(pairs, reverse=True):
        if i in used_t or j in used_p:
            continue
        used_t.add(i); used_p.add(j); n += 1
    return n


def prf(tp, nc, nt):
    p = tp / nc if nc else np.nan
    r = tp / nt if nt else np.nan
    f = 2 * p * r / (p + r) if (nc and p + r > 0) else 0.0
    return p, r, f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    B, O = Path(a.bundle), Path(a.out); O.mkdir(parents=True, exist_ok=True)
    h = pd.read_csv(B / "09_metrics/simulation1/harmonized_calls.csv", dtype={"chrom": str})
    h = h[[b in VALID.get(t, set()) for t, b in zip(h.tool, h.analysis_bin_bp)]]
    h["type"] = h["type"].str.lower()
    m = pd.read_csv(B / "09_metrics/simulation1/simulation1_overlap_master.csv")
    m = m[[b in VALID.get(t, set()) for t, b in zip(m.tool, m.analysis_bin_bp)]]

    # ---- 1. negative-control characterisation
    neg = h[h.sample_class == "negative"]
    char = (neg.groupby(["tool", "analysis_bin_bp", "background", "depth_x"])
            .agg(n_calls=("size_bp", "size"), median_size_kb=("size_bp", lambda x: x.median() / 1e3),
                 loss_fraction=("type", lambda x: (x == "loss").mean()), n_chromosomes=("chrom", "nunique"))
            .reset_index())
    char.to_csv(O / "negative_call_characteristics.tsv", sep="\t", index=False)
    rec = []
    for (t, b, bg), g in neg.groupby(["tool", "analysis_bin_bp", "background"]):
        ds = sorted(g.depth_x.unique())
        for d1, d2 in zip(ds[:-1], ds[1:]):
            x, y = g[g.depth_x == d1], g[g.depth_x == d2]
            if not len(x):
                continue
            mk = {}
            for c, s, e in merge(list(zip(y.chrom, y.start, y.end))):
                mk.setdefault(c, []).append((c, s, e))
            hit = np.mean([frac_in(mk, r.chrom, r.start, r.end) > 0 for r in x.itertuples()])
            rec.append({"tool": t, "bin_bp": b, "background": bg, "depth_from": d1, "depth_to": d2,
                        "n_from": len(x), "n_to": len(y), "fraction_recurring": hit})
    rec = pd.DataFrame(rec); rec.to_csv(O / "negative_call_recurrence.tsv", sep="\t", index=False)

    # ---- 2. self-normal masks (per tool x bin x background, recurrent at >= 2 depths)
    def build(g):
        keep = []
        for d, x in g.groupby("depth_x"):
            other = g[g.depth_x != d]
            om = {}
            for c, s_, e in merge(list(zip(other.chrom, other.start, other.end))):
                om.setdefault(c, []).append((c, s_, e))
            keep += [(r.chrom, r.start, r.end) for r in x.itertuples() if frac_in(om, r.chrom, r.start, r.end) > 0]
        mk = {}
        for c, s_, e in merge(keep):
            mk.setdefault(c, []).append((c, s_, e))
        return mk

    lodo = {}
    for (t, b, bg), g in neg.groupby(["tool", "analysis_bin_bp", "background"]):
        for d in g.depth_x.unique():
            lodo[(t, b, bg, d)] = build(g[g.depth_x != d])
    masks, msum = {}, []
    for (t, b, bg), g in neg.groupby(["tool", "analysis_bin_bp", "background"]):
        keep = []
        for d, x in g.groupby("depth_x"):
            other = g[g.depth_x != d]
            om = {}
            for c, s_, e in merge(list(zip(other.chrom, other.start, other.end))):
                om.setdefault(c, []).append((c, s_, e))
            keep += [(r.chrom, r.start, r.end) for r in x.itertuples() if frac_in(om, r.chrom, r.start, r.end) > 0]
        mg = merge(keep)
        mk = {}
        for c, s_, e in mg:
            mk.setdefault(c, []).append((c, s_, e))
        masks[(t, b, bg)] = mk
        msum.append({"tool": t, "bin_bp": b, "background": bg, "mask_regions": len(mg),
                     "mask_Mb": sum(e - s_ for _, s_, e in mg) / 1e6})
        pd.DataFrame(mg, columns=["chrom", "start", "end"]).to_csv(
            O / f"self_normal_mask_{t}_{b}_{bg}.bed", sep="\t", index=False, header=False)
    msum = pd.DataFrame(msum); msum.to_csv(O / "mask_summary.tsv", sep="\t", index=False)

    # ---- 3. re-score
    truth = {}
    for r in range(1, 6):
        t = pd.read_csv(B / f"truth/rep{r}.truth.tsv", sep="\t", dtype={"#chrom": str})
        truth[r] = [(str(c), int(s), int(e), ty.lower()) for c, s, e, ty in
                    zip(t["#chrom"], t.start, t.end, t.cnv_type)]
    pos = h[h.sample_class == "positive"]
    keys = m[m.sample_class == "positive"][["tool", "sample", "analysis_bin_bp", "background", "depth_x",
                                             "replicate", "F1_50pct"]]
    rows = []
    for k in keys.itertuples(index=False):
        g = pos[(pos.tool == k.tool) & (pos["sample"] == k.sample) & (pos.analysis_bin_bp == k.analysis_bin_bp)]
        calls = [(c, s, e, ty) for c, s, e, ty in zip(g.chrom, g.start, g.end, g["type"])]
        tr = truth[int(k.replicate)]
        kept = [c for c in calls if frac_in(masks.get((k.tool, k.analysis_bin_bp, k.background), {}), *c[:3]) <= 0.5]
        lm = lodo.get((k.tool, k.analysis_bin_bp, k.background, k.depth_x), {})
        kept_lodo = [c for c in calls if frac_in(lm, *c[:3]) <= 0.5]
        out = {"tool": k.tool, "bin_bp": k.analysis_bin_bp, "sample": k.sample, "background": k.background,
               "depth": k.depth_x, "replicate": k.replicate, "F1_bundle": k.F1_50pct}
        for lab, cl in (("unmasked", calls), ("masked", kept), ("masked_lodo", kept_lodo)):
            p, r, f = prf(match(tr, cl), len(cl), len(tr))
            out.update({f"n_calls_{lab}": len(cl), f"precision_{lab}": p, f"recall_{lab}": r, f"F1_{lab}": f})
        rows.append(out)
    rs = pd.DataFrame(rows)
    rs.to_csv(O / "masked_rescore_per_sample.tsv", sep="\t", index=False)
    agree = np.isclose(rs.F1_unmasked, rs.F1_bundle.fillna(0), atol=1e-6).mean()
    print(f"validation: unmasked re-score equals bundle F1 in {agree:.1%} of {len(rs)} configurations")
    print(msum.to_string(index=False))
    s = rs.groupby(["tool", "bin_bp", "background"])[["F1_unmasked", "F1_masked", "F1_masked_lodo",
                                                       "precision_unmasked", "precision_masked_lodo",
                                                       "recall_unmasked", "recall_masked_lodo"]].median().round(3)
    s.to_csv(O / "masked_rescore_summary.tsv", sep="\t")
    print(s.to_string())


if __name__ == "__main__":
    main()
