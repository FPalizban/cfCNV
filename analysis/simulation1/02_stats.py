#!/usr/bin/env python3
"""
cfCNV Simulation 1 - summary statistics (Methods 2.10).

  python 02_stats.py --tables tables/ --out stats/

1. Primary configuration per tool = bin size with the highest median event F1
   on the matched-depth set (0.1-1x, all three backgrounds, positives).
2. Table 1: per tool at its primary bin - median [95% bootstrap CI] precision,
   recall, F1, bin-level MCC; negative-control false calls/Gb; breakpoint error;
   runtime per run.
3. Paired two-sided Wilcoxon signed-rank tests of F1 between tools (paired on
   sample, matched depth), Benjamini-Hochberg FDR. Repeated within IH02+IH03.
4. Friedman test of the bin-size effect within each tool (paired on sample).
5. Spearman correlation of F1 with depth within each tool x background, BH FDR.
6. Size limit: smallest simulated CNV size with recall >= 0.5 per tool x bin.
"""
import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare, spearmanr, wilcoxon

SEG_TOOLS = ["QDNAseq", "ichorCNA", "BICseq2", "CNV-MEANN", "dudeML", "CNVkit", "Spectre", "LYCEUM",
             "ichorCNA_default", "CNVkit_flat", "Spectre_1kb"]
MAIN_TOOLS = SEG_TOOLS[:8]


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def boot(x, f=np.median, n=5000, seed=11):
    x = np.asarray(pd.Series(x).dropna(), float)
    if len(x) < 2:
        return np.nan, np.nan
    r = np.random.default_rng(seed)
    v = f(r.choice(x, (n, len(x))), axis=1)
    return np.percentile(v, 2.5), np.percentile(v, 97.5)


def fmt(x, lo, hi, d=2):
    return "NA" if pd.isna(x) else f"{x:.{d}f} [{lo:.{d}f}-{hi:.{d}f}]"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    T, O = Path(a.tables), Path(a.out); O.mkdir(parents=True, exist_ok=True)
    pc = pd.read_csv(T / "per_config.tsv", sep="\t")
    td = pd.read_csv(T / "truth_detection.tsv", sep="\t")
    rt = pd.read_csv(T / "runtime.tsv", sep="\t")
    pos = pc[(pc["class"] == "positive") & pc.tool.isin(SEG_TOOLS)]
    neg = pc[(pc["class"] == "negative") & pc.tool.isin(SEG_TOOLS)]

    # 1. primary configuration
    best = (pos[pos.matched_depth].groupby(["tool", "bin_bp"]).F1.median().reset_index()
            .sort_values("F1").groupby("tool").tail(1).set_index("tool").bin_bp.astype(int))
    best.rename("bin_bp").to_csv(O / "primary_configs.tsv", sep="\t")
    prim = lambda d: d[[best.get(t) == b for t, b in zip(d.tool, d.bin_bp)]]
    pp, nn = prim(pos), prim(neg)
    pm = pp[pp.matched_depth]

    # 2. Table 1
    rows = []
    for t in SEG_TOOLS:
        g, gn = pm[pm.tool == t], nn[(nn.tool == t) & nn.matched_depth]
        r = {"tool": t, "primary_bin_kb": best[t] / 1000, "n_positive": len(g), "n_negative": len(gn)}
        for m in ("precision", "recall", "F1", "bin_MCC"):
            v = g[m].median(); r[m] = fmt(v, *boot(g[m])); r[f"{m}_median"] = v
        v = gn.false_calls_per_Gb.median(); r["neg_false_calls_per_Gb"] = fmt(v, *boot(gn.false_calls_per_Gb), 1)
        r["neg_false_altered_pct"] = round(100 * gn.false_altered_fraction.median(), 2)
        r["breakpoint_error_kb_median"] = round(g.breakpoint_error_bp.median() / 1000, 1)
        w = rt[(rt.tool == t) & rt.state.str.startswith("COMPLETED")].wall_s
        r["runtime_min_median"] = round(w.median() / 60, 1) if len(w) else np.nan
        r["failed_runs"] = int((rt[rt.tool == t].state == "FAILED").sum())
        rows.append(r)
    t1 = pd.DataFrame(rows)
    t1.to_csv(O / "table1_primary_summary.tsv", sep="\t", index=False)
    print(t1[["tool", "primary_bin_kb", "precision", "recall", "F1", "bin_MCC",
              "neg_false_calls_per_Gb", "neg_false_altered_pct", "breakpoint_error_kb_median",
              "runtime_min_median", "failed_runs"]].to_string(index=False))

    # 3. pairwise Wilcoxon
    def pairwise(d, label):
        w = d.pivot_table(index="sample_id", columns="tool", values="F1")
        out = []
        for x, y in combinations([t for t in MAIN_TOOLS if t in w], 2):
            dd = (w[x] - w[y]).dropna()
            if len(dd) < 6 or (dd == 0).all():
                continue
            s = wilcoxon(dd)
            out.append({"subset": label, "tool_A": x, "tool_B": y, "n": len(dd),
                        "median_F1_diff": dd.median(), "ci_low": boot(dd)[0], "ci_high": boot(dd)[1],
                        "W": s.statistic, "p": s.pvalue})
        r = pd.DataFrame(out); r["q_BH"] = bh(r.p); return r
    pw = pd.concat([pairwise(pm, "matched depth, all backgrounds"),
                    pairwise(pm[pm.background != "IH01"], "matched depth, IH02+IH03")])
    pw.to_csv(O / "pairwise_wilcoxon_F1.tsv", sep="\t", index=False)
    print("\n", pw.round(4).to_string(index=False))

    # 4. Friedman bin effect
    fr = []
    for t, g in pos[pos.matched_depth].groupby("tool"):
        w = g.pivot_table(index="sample_id", columns="bin_bp", values="F1").dropna()
        if w.shape[1] >= 3:
            s = friedmanchisquare(*[w[c] for c in w.columns])
            fr.append({"tool": t, "bins_kb": ",".join(f"{c/1000:g}" for c in w.columns), "n": len(w),
                       "chi2": s.statistic, "p": s.pvalue,
                       **{f"median_F1_{c/1000:g}kb": w[c].median() for c in w.columns}})
    fr = pd.DataFrame(fr); fr.to_csv(O / "friedman_bin_effect.tsv", sep="\t", index=False)
    print("\n", fr.round(4).to_string(index=False))

    # 5. depth trend within background
    sp = []
    for (t, bg), g in pp.groupby(["tool", "background"]):
        if g.depth.nunique() < 3:
            continue
        for m in ("F1", "precision", "recall"):
            r = spearmanr(g.depth, g[m], nan_policy="omit")
            sp.append({"tool": t, "background": bg, "metric": m, "n": g[m].notna().sum(),
                       "depths": f"{g.depth.min():g}-{g.depth.max():g}x", "rho": r.statistic, "p": r.pvalue})
    sp = pd.DataFrame(sp); sp["q_BH"] = bh(sp.p.fillna(1))
    sp.to_csv(O / "depth_trend_spearman.tsv", sep="\t", index=False)
    print("\n", sp[sp.metric == "F1"].round(4).to_string(index=False))

    # 6. size limit
    rec = td.groupby(["tool", "bin_bp", "size_bp"]).detected.mean().reset_index()
    lim = (rec[rec.detected >= 0.5].groupby(["tool", "bin_bp"]).size_bp.min()
           .rename("smallest_size_recall_ge_0.5_bp").reset_index())
    rec.to_csv(O / "recall_by_size_all_depths.tsv", sep="\t", index=False)
    lim.to_csv(O / "size_detection_limit.tsv", sep="\t", index=False)
    print("\n", lim.to_string(index=False))


if __name__ == "__main__":
    main()
