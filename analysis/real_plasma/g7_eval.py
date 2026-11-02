#!/usr/bin/env python3
"""Cancer plasma (GSE71378) evaluation: 1-Mb window states per caller, majority consensus,
leave-one-out agreement, TCGA concordance against the non-cancer baseline, cancer-gene direction,
depth robustness. usage: g7_eval.py <cancer_calls.tsv> <healthy_calls_hg19.tsv> <refs_dir> <out_dir>"""
import sys, numpy as np, pandas as pd
from pathlib import Path
CALLS, HEALTHY, R, O = sys.argv[1:5]; R = Path(R); O = Path(O); O.mkdir(parents=True, exist_ok=True)
W = 1_000_000; AUTO = [str(i) for i in range(1, 23)]
PRIM = [("QDNAseq", "100kb"), ("ichorCNA", "1Mb low-TF"), ("BIC-seq2", "10kb"), ("CNVkit", "30kb"), ("Spectre", "100kb")]
COHORT = {"IC15": "LUSC", "IC17": "LIHC", "IC20": "LUSC", "IC35": "BRCA", "IC37": "COADREAD"}   # IC15 (SCLC): LUSC as nearest TCGA proxy
TYPE = {"IC15": "small-cell lung", "IC17": "hepatocellular", "IC20": "lung squamous", "IC35": "breast (DCIS)", "IC37": "colorectal"}
arms = pd.read_csv(R / "arms_hg19.tsv", sep="\t", dtype={"chrom": str}); size = arms.groupby("chrom").end.max().to_dict()
win = [(c, s) for c in AUTO for s in range(0, size[c], W)]; idx = {k: i for i, k in enumerate(win)}; N = len(win)
armof = np.array([arms[(arms.chrom == c) & (arms.start <= s) & (arms.end > s)].arm.iloc[0] if len(arms[(arms.chrom == c) & (arms.start <= s) & (arms.end > s)]) else f"{c}?" for c, s in win])

def states(df):
    """signed state per 1-Mb window: +1 gain / -1 loss where calls cover >= 50% of the window"""
    g = np.zeros(N); l = np.zeros(N)
    for c, s, e, t in zip(df.chrom.astype(str).str.replace("chr", ""), df.start, df.end, df.type):
        if c not in size: continue
        for w in range(s // W, (e - 1) // W + 1):
            k = idx.get((c, w * W)); 
            if k is None: continue
            ov = min(e, (w + 1) * W) - max(s, w * W)
            if ov > 0: (g if t == "gain" else l)[k] += ov
    st = np.zeros(N); st[g >= W / 2] = 1; st[l >= W / 2] = -1; st[(g >= W / 2) & (l >= W / 2)] = np.sign(g - l)[(g >= W / 2) & (l >= W / 2)]
    return st

# TCGA signed frequency per cohort
tw = pd.read_csv(R / "tcga_window_freq_1Mb.tsv", sep="\t", dtype={"chrom": str})
tc = {}
for name, mem in {"LIHC": ["LIHC"], "LUSC": ["LUSC"], "COADREAD": ["COAD", "READ"], "BRCA": ["BRCA"]}.items():
    t = tw[tw.cohort.isin(mem)].copy(); t["g"] = t.gain_freq * t.n; t["l"] = t.loss_freq * t.n
    a = t.groupby(["chrom", "start"])[["g", "l", "n"]].sum()
    v = np.full(N, np.nan)
    for (c, s), r in a.iterrows():
        k = idx.get((c, s)); 
        if k is not None: v[k] = (r.g - r.l) / r.n
    tc[name] = v
def tcga_r(st, coh):
    v = tc[coh]; ok = ~np.isnan(v)
    if np.all(st[ok] == 0): return 0.0
    return float(np.corrcoef(st[ok], v[ok])[0, 1])

C = pd.read_csv(CALLS, sep="\t", dtype={"chrom": str})
S = {}
for (sid, base, lvl, t, cfg), g in C.groupby(["sample_id", "base", "level", "tool", "config"]):
    if (t, cfg) in PRIM: S[(base, lvl, t)] = states(g)
bases = sorted({k[0] for k in S}); cancers = [b for b in bases if b in COHORT]
for b in bases:
    for lvl in ("full", "d1", "d0p1"):
        for t, _ in PRIM:
            if (b, lvl, t) not in S and any(k[:2] == (b, lvl) for k in S): S[(b, lvl, t)] = np.zeros(N)   # caller made no call

# consensus and leave-one-out
rows, cons = [], {}
for b in bases:
    M = np.array([S[(b, "full", t)] for t, _ in PRIM])
    ng, nl = (M == 1).sum(0), (M == -1).sum(0)
    cons[b] = np.where(ng >= 3, 1, np.where(nl >= 3, -1, 0))
    for i, (t, _) in enumerate(PRIM):
        o = np.delete(M, i, 0); og, ol = (o == 1).sum(0), (o == -1).sum(0)
        ref = np.where((og >= 2) & (og > ol), 1, np.where((ol >= 2) & (ol > og), -1, 0)); x = M[i]
        tp = ((x != 0) & (x == ref)).sum(); fp = ((x != 0) & (x != ref)).sum(); fn = ((ref != 0) & (x != ref)).sum()
        p = tp / (tp + fp) if tp + fp else np.nan; r = tp / (tp + fn) if tp + fn else np.nan
        rows.append({"sample": b, "tool": t, "altered_pct": 100 * (x != 0).mean(), "loo_ref_altered_pct": 100 * (ref != 0).mean(),
                     "precision": p, "recall": r, "F1": 2 * p * r / (p + r) if p and r and p + r else np.nan,
                     "tcga_r": tcga_r(x, COHORT[b]) if b in COHORT else np.nan})
L = pd.DataFrame(rows); L.to_csv(O / "cancer_loo_and_tcga.tsv", sep="\t", index=False)
pd.DataFrame({"sample": b, "consensus_altered_pct": 100 * (cons[b] != 0).mean(), "consensus_tcga_r": tcga_r(cons[b], COHORT[b]) if b in COHORT else np.nan} for b in bases).to_csv(O / "consensus_summary.tsv", sep="\t", index=False)

# healthy baseline: same correlation for the 25 non-cancer samples (native depth, calls lifted to GRCh37)
H = pd.read_csv(HEALTHY, sep="\t", dtype={"chrom": str}); H = H[H.level == "full"]
hb = []
for (base, t, cfg), g in H.groupby(["base", "tool", "config"]):
    st = states(g)
    for coh in tc: hb.append({"sample": base, "tool": t, "cohort": coh, "tcga_r": tcga_r(st, coh), "altered_pct": 100 * (st != 0).mean()})
HB = pd.DataFrame(hb); HB.to_csv(O / "healthy_tcga_baseline.tsv", sep="\t", index=False)
# per tool, z-score of each cancer sample against the non-cancer distribution for its cohort
L["healthy_mean"] = [HB[(HB.tool == t) & (HB.cohort == COHORT.get(s, "LUSC"))].tcga_r.mean() for s, t in zip(L["sample"], L.tool)]
L["healthy_sd"] = [HB[(HB.tool == t) & (HB.cohort == COHORT.get(s, "LUSC"))].tcga_r.std() for s, t in zip(L["sample"], L.tool)]
L["healthy_max"] = [HB[(HB.tool == t) & (HB.cohort == COHORT.get(s, "LUSC"))].tcga_r.max() for s, t in zip(L["sample"], L.tool)]
L.to_csv(O / "cancer_loo_and_tcga.tsv", sep="\t", index=False)

# cancer-gene direction (consensus and per caller, native depth)
genes = pd.read_csv(R / "genes_hg19.tsv", sep="\t", dtype={"chrom": str}); cg = pd.read_csv(R / "cancer_genes.tsv", sep="\t")
cg = cg[(cg.oncogene + cg.tsg) == 1]; genes = genes.merge(cg, on="gene").drop_duplicates("gene")
gk = [idx.get((c, (s // W) * W)) for c, s in zip(genes.chrom, (genes.start + genes.end) // 2)]; genes["k"] = gk; genes = genes.dropna(subset=["k"]); genes["k"] = genes.k.astype(int)
dr = []
for b in cancers:
    for t, st in [("consensus", cons[b])] + [(t, S[(b, "full", t)]) for t, _ in PRIM]:
        s = st[genes.k.values]; onc = genes.oncogene.values == 1; tsg = genes.tsg.values == 1
        conc = ((s == 1) & onc).sum() + ((s == -1) & tsg).sum(); disc = ((s == -1) & onc).sum() + ((s == 1) & tsg).sum()
        dr.append({"sample": b, "tool": t, "concordant": int(conc), "discordant": int(disc), "concordant_frac": conc / (conc + disc) if conc + disc else np.nan,
                   "oncogenes_gained": ",".join(genes[(s == 1) & onc].gene.head(12)), "tsg_lost": ",".join(genes[(s == -1) & tsg].gene.head(12))})
pd.DataFrame(dr).to_csv(O / "cancer_gene_direction.tsv", sep="\t", index=False)

# depth robustness: agreement of 1x and 0.1x with the caller's own native-depth profile
dp = []
for b in cancers:
    for t, _ in PRIM:
        ref = S[(b, "full", t)]
        for lvl in ("d1", "d0p1"):
            x = S.get((b, lvl, t))
            if x is None: continue
            tp = ((x != 0) & (x == ref)).sum(); fp = ((x != 0) & (x != ref)).sum(); fn = ((ref != 0) & (x != ref)).sum()
            p = tp / (tp + fp) if tp + fp else np.nan; r = tp / (tp + fn) if tp + fn else np.nan
            dp.append({"sample": b, "tool": t, "level": lvl, "precision": p, "recall": r, "F1": 2 * p * r / (p + r) if p and r else 0.0,
                       "consensus_recall": ((x != 0) & (x == cons[b])).sum() / max(1, (cons[b] != 0).sum())})
pd.DataFrame(dp).to_csv(O / "depth_robustness.tsv", sep="\t", index=False)

# arm-level consensus events (>= 50% of the arm's windows altered in one direction)
ar = []
for b in bases:
    for a in pd.unique(armof):
        m = armof == a; v = cons[b][m]
        if m.sum() < 3: continue
        for d, lab in ((1, "gain"), (-1, "loss")):
            if (v == d).mean() >= 0.5: ar.append({"sample": b, "type": TYPE.get(b, "healthy"), "arm": a, "event": lab, "frac_arm": (v == d).mean()})
pd.DataFrame(ar).to_csv(O / "consensus_arm_events.tsv", sep="\t", index=False)
np.save(O / "window_states.npy", {"win": win, "S": {f"{k[0]}|{k[1]}|{k[2]}": v for k, v in S.items()}, "cons": cons, "tcga": tc}, allow_pickle=True)
print("done")
