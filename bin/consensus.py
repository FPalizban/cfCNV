#!/usr/bin/env python3
"""Window-level consensus and caller agreement, as used for the cancer-plasma evaluation in the cfCNV paper.

 * every caller's calls are converted to a signed state per window (default 1 Mb): +1 gain / -1 loss when
   calls of that direction cover >= 50% of the window
 * majority consensus: >= --min-callers callers agree on the direction
 * leave-one-out agreement: each caller vs. the windows where >= 2 of the other callers agree and agreeing
   callers outnumber disagreeing ones -> precision, recall, F1
 * depth robustness: for down-sampled levels, F1 of the profile against the same caller's native ('full') profile
 * optional TCGA concordance: Pearson r between window states and the signed TCGA PanCanAtlas frequency
   (gain - loss) of the cohort given in the samplesheet column tcga_cohort (e.g. LUSC, COAD+READ);
   needs --tcga tcga_window_freq_1Mb.tsv (GRCh37 windows; lift GRCh38 calls first or use GRCh37 data)

usage: consensus.py --calls all_calls.tsv --fai ref.fa.fai --primary "QDNAseq:100kb,..." [--samples samplesheet.csv]
                    [--tcga tcga_window_freq_1Mb.tsv] [--window 1000000] [--min-callers 3] [--outdir .]
"""
import argparse
from pathlib import Path
import warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore", category=RuntimeWarning)

AUTO = [str(i) for i in range(1, 23)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True); ap.add_argument("--fai", required=True)
    ap.add_argument("--primary", required=True, help="tool:config pairs used for the consensus")
    ap.add_argument("--samples"); ap.add_argument("--tcga")
    ap.add_argument("--window", type=int, default=1_000_000); ap.add_argument("--min-callers", type=int, default=3)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args(); O = Path(a.outdir); O.mkdir(parents=True, exist_ok=True); W = a.window

    size = {}
    for ln in open(a.fai):
        c, L = ln.split()[:2]; c = c[3:] if c.startswith("chr") else c
        if c in AUTO: size[c] = int(L)
    win = [(c, s) for c in AUTO if c in size for s in range(0, size[c], W)]
    idx = {k: i for i, k in enumerate(win)}; N = len(win)
    prim = [tuple(x.split(":", 1)) for x in a.primary.split(",") if ":" in x]
    tools = [t for t, _ in prim]

    def states(df):
        g = np.zeros(N); l = np.zeros(N)
        for c, s, e, t in zip(df.chrom.astype(str), df.start, df.end, df.type):
            if c not in size: continue
            for w in range(s // W, (e - 1) // W + 1):
                k = idx.get((c, w * W))
                if k is None: continue
                ov = min(e, (w + 1) * W) - max(s, w * W)
                if ov > 0: (g if t == "gain" else l)[k] += ov
        st = np.zeros(N); hg, hl = g >= W / 2, l >= W / 2
        st[hg] = 1; st[hl] = -1; both = hg & hl; st[both] = np.sign(g - l)[both]
        return st

    def prf(x, ref):
        tp = ((x != 0) & (x == ref)).sum(); fp = ((x != 0) & (x != ref)).sum(); fn = ((ref != 0) & (x != ref)).sum()
        p = tp / (tp + fp) if tp + fp else np.nan; r = tp / (tp + fn) if tp + fn else np.nan
        f = 2 * p * r / (p + r) if (p == p and r == r and p + r) else (0.0 if tp + fp + fn else np.nan)
        return p, r, f

    C = pd.read_csv(a.calls, sep="\t", dtype={"chrom": str, "level": str, "base": str})
    C["chrom"] = C.chrom.str.replace("^chr", "", regex=True)
    runs = pd.read_csv(Path(a.calls).with_name("call_burden.tsv"), sep="\t", dtype=str) \
        if Path(a.calls).with_name("call_burden.tsv").exists() else C[["base", "level", "tool", "config"]].drop_duplicates()
    S = {}
    for r in runs[["base", "level", "tool", "config"]].drop_duplicates().itertuples():
        if (r.tool, r.config) in prim:
            S[(r.base, r.level, r.tool)] = states(C[(C.base == r.base) & (C.level == r.level) & (C.tool == r.tool) & (C.config == r.config)])

    cohort = {}
    if a.samples:
        sm = pd.read_csv(a.samples, dtype=str).fillna("")
        if "tcga_cohort" in sm.columns: cohort = {s: c for s, c in zip(sm["sample"], sm["tcga_cohort"]) if c}
    tc = {}
    if a.tcga and cohort:
        tw = pd.read_csv(a.tcga, sep="\t", dtype={"chrom": str})
        for name in set(cohort.values()):
            mem = name.split("+"); t = tw[tw.cohort.isin(mem)].copy()
            if t.empty: print(f"warning: TCGA cohort {name} not found"); continue
            t["g"] = t.gain_freq * t.n; t["l"] = t.loss_freq * t.n
            agg = t.groupby(["chrom", "start"])[["g", "l", "n"]].sum(); v = np.full(N, np.nan)
            for (c, s), x in agg.iterrows():
                k = idx.get((c, s))
                if k is not None: v[k] = (x.g - x.l) / x.n
            tc[name] = v

    def tcga_r(st, coh):
        if coh not in tc: return np.nan
        v = tc[coh]; ok = ~np.isnan(v)
        if np.all(st[ok] == 0): return 0.0
        return float(np.corrcoef(st[ok], v[ok])[0, 1])

    rows, cons_rows, summ, rob, states_out = [], [], [], [], []
    for (b, lvl) in sorted({k[:2] for k in S}):
        avail = [t for t in tools if (b, lvl, t) in S]
        M = np.array([S[(b, lvl, t)] for t in avail])
        ng, nl = (M == 1).sum(0), (M == -1).sum(0)
        cons = np.where(ng >= a.min_callers, 1, np.where(nl >= a.min_callers, -1, 0))
        for i, t in enumerate(avail):
            o = np.delete(M, i, 0); og, ol = (o == 1).sum(0), (o == -1).sum(0)
            ref = np.where((og >= 2) & (og > ol), 1, np.where((ol >= 2) & (ol > og), -1, 0))
            p, r, f = prf(M[i], ref)
            rows.append({"base": b, "level": lvl, "tool": t, "altered_pct": 100 * (M[i] != 0).mean(),
                         "loo_reference_altered_pct": 100 * (ref != 0).mean(), "precision": p, "recall": r, "F1": f,
                         "tcga_cohort": cohort.get(b, ""), "tcga_r": tcga_r(M[i], cohort.get(b))})
        summ.append({"base": b, "level": lvl, "n_callers": len(avail), "consensus_altered_pct": 100 * (cons != 0).mean(),
                     "tcga_cohort": cohort.get(b, ""), "consensus_tcga_r": tcga_r(cons, cohort.get(b))})
        # consensus segments (merge adjacent windows with the same state)
        cur = None
        for k, (c, s) in enumerate(win):
            v = cons[k]
            if cur and (cur[0] != c or cur[3] != v or cur[2] != s):
                cons_rows.append(cur + [(b, lvl)]); cur = None
            if v != 0 and cur is None: cur = [c, s, min(s + W, size[c]), v]
            elif v != 0: cur[2] = min(s + W, size[c])
        if cur: cons_rows.append(cur + [(b, lvl)])
        for t in avail:
            states_out.append(pd.Series(S[(b, lvl, t)], name=f"{b}|{lvl}|{t}"))
        states_out.append(pd.Series(cons, name=f"{b}|{lvl}|consensus"))
        if lvl != "full":
            for t in avail:
                if (b, "full", t) in S:
                    p, r, f = prf(S[(b, lvl, t)], S[(b, "full", t)])
                    rob.append({"base": b, "level": lvl, "tool": t, "precision": p, "recall": r, "F1_vs_native": f})

    pd.DataFrame(rows).to_csv(O / "caller_agreement.tsv", sep="\t", index=False)
    pd.DataFrame(summ).to_csv(O / "consensus_summary.tsv", sep="\t", index=False)
    pd.DataFrame([{"base": x[4][0], "level": x[4][1], "chrom": x[0], "start": x[1], "end": x[2],
                   "type": "gain" if x[3] > 0 else "loss"} for x in cons_rows],
                 columns=["base", "level", "chrom", "start", "end", "type"]).to_csv(O / "consensus_calls.tsv", sep="\t", index=False)
    pd.DataFrame(rob, columns=["base", "level", "tool", "precision", "recall", "F1_vs_native"]).to_csv(
        O / "depth_robustness.tsv", sep="\t", index=False)
    if states_out:
        st = pd.concat([pd.DataFrame({"chrom": [c for c, _ in win], "start": [s for _, s in win]})] +
                       [x.to_frame() for x in states_out], axis=1)
        st.to_csv(O / "window_states.tsv.gz", sep="\t", index=False)
    print(f"{len(summ)} sample/depth profiles; {len(rows)} caller comparisons; {len(rob)} depth comparisons")


if __name__ == "__main__":
    main()
