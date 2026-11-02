#!/usr/bin/env python3
"""
cfCNV real-data annotation and evaluation (used for both the non-cancer and the cancer cohort).

usage: annotate_calls.py <calls.tsv> <refs_dir> <out_dir> [--lift hg38ToHg19.over.chain.gz]
calls.tsv columns: sample_id base level tool config chrom start end type(gain/loss)

For every call (GRCh37; GRCh38 calls are lifted over first):
  * explained fraction: bp overlapping common germline CNVs (gnomAD-SV AF>=1%), segmental duplications,
    the ENCODE blacklist or pericentromeric/acrocentric regions (+-3 Mb of the centromere; 13p 14p 15p 21p 22p)
  * cancer-likeness: overlap with OncoKB cancer genes, and the TCGA same-direction alteration frequency of
    the called 1-Mb windows (LIHC, LUSC, COAD+READ, BRCA)
Each observed value is compared with a size-preserving random placement of the same calls on the same
chromosome (N_PERM permutations) -> enrichment = observed / expected.
"""
import argparse, bisect, gzip
from collections import defaultdict
from pathlib import Path
import numpy as np, pandas as pd

AUTO = [str(i) for i in range(1, 23)]; N_PERM = 50; W = 1_000_000
COHORTS = {"LIHC": ["LIHC"], "LUSC": ["LUSC"], "COADREAD": ["COAD", "READ"], "BRCA": ["BRCA"]}

class IV:  # merged intervals with prefix sums for fast overlap
    def __init__(self, rows):
        d = defaultdict(list)
        for c, s, e in rows: d[str(c).replace("chr", "")].append((int(s), int(e)))
        self.m = {}
        for c, L in d.items():
            L.sort(); M = []
            for s, e in L:
                if M and s <= M[-1][1]: M[-1][1] = max(M[-1][1], e)
                else: M.append([s, e])
            S = np.array([x[0] for x in M]); E = np.array([x[1] for x in M]); cum = np.concatenate([[0], np.cumsum(E - S)])
            self.m[c] = (S, E, cum)
    def bp(self, c, s, e):
        if c not in self.m: return 0
        S, E, cum = self.m[c]
        i = bisect.bisect_right(E, s); j = bisect.bisect_left(S, e)
        if j <= i: return 0
        tot = cum[j] - cum[i]
        tot -= max(0, s - S[i]); tot -= max(0, E[j - 1] - e)
        return max(0, int(tot))

def load_refs(R):
    arms = pd.read_csv(R / "arms_hg19.tsv", sep="\t", dtype={"chrom": str})
    size = arms.groupby("chrom").end.max().to_dict(); cen = arms[arms.arm.str.endswith("q")].set_index("chrom").start.to_dict()
    peri = [(c, max(0, cen[c] - 3_000_000), cen[c] + 3_000_000) for c in AUTO] + [(c, 0, cen[c]) for c in ("13", "14", "15", "21", "22")]
    g = pd.read_csv(R / "gnomad_sv_common_cnv.bed", sep="\t", header=None, dtype={0: str})
    sd = pd.read_csv(R / "segdups_hg19.bed", sep="\t", header=None, dtype={0: str})
    bl = pd.read_csv(R / "blacklist_hg19.bed", sep="\t", header=None, dtype={0: str})
    cats = {"germline_CNV": IV(g[[0, 1, 2]].values), "segdup": IV(sd[[0, 1, 2]].values), "blacklist": IV(bl[[0, 1, 2]].values), "pericentromeric": IV(peri)}
    cats["any"] = IV([r for k in ("germline_CNV", "segdup", "blacklist", "pericentromeric") for r in
                      ([tuple(x) for x in g[[0, 1, 2]].values] if k == "germline_CNV" else [tuple(x) for x in sd[[0, 1, 2]].values] if k == "segdup"
                       else [tuple(x) for x in bl[[0, 1, 2]].values] if k == "blacklist" else peri)])
    genes = pd.read_csv(R / "genes_hg19.tsv", sep="\t", dtype={"chrom": str})
    cg = pd.read_csv(R / "cancer_genes.tsv", sep="\t"); cgs = set(cg.gene)
    cgene = genes[genes.gene.isin(cgs)]
    tw = pd.read_csv(R / "tcga_window_freq_1Mb.tsv", sep="\t", dtype={"chrom": str})
    freq = {}
    for name, members in COHORTS.items():
        t = tw[tw.cohort.isin(members)].copy(); t["g"] = t.gain_freq * t.n; t["l"] = t.loss_freq * t.n
        a = t.groupby(["chrom", "start"])[["g", "l", "n"]].sum()
        freq[name] = {(c, s): (r.g / r.n, r.l / r.n) for (c, s), r in a.iterrows()}
    return dict(size=size, cats=cats, genes=genes, cgene=IV(cgene[["chrom", "start", "end"]].values), cg=cg, freq=freq)

def tcga_score(rows, freq):
    """bp-weighted mean TCGA same-direction frequency over the 1-Mb windows touched by the calls"""
    num = den = 0.0
    for c, s, e, t in rows:
        for w in range(s // W, (e - 1) // W + 1):
            f = freq.get((c, w * W))
            if f is None: continue
            ov = min(e, (w + 1) * W) - max(s, w * W)
            num += ov * (f[0] if t == "gain" else f[1]); den += ov
    return num / den if den else np.nan

def metrics(rows, ref):
    L = sum(e - s for _, s, e, _ in rows)
    if not L: return {}
    m = {f"frac_{k}": sum(ref["cats"][k].bp(c, s, e) for c, s, e, _ in rows) / L for k in ref["cats"]}
    m["calls_with_cancer_gene"] = np.mean([ref["cgene"].bp(c, s, e) > 0 for c, s, e, _ in rows])
    for k, f in ref["freq"].items(): m[f"tcga_{k}"] = tcga_score(rows, f)
    return m

def shuffle(rows, size, rng):
    out = []
    for c, s, e, t in rows:
        n = e - s; hi = max(1, size[c] - n); a = int(rng.integers(0, hi)); out.append((c, a, a + n, t))
    return out

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("calls"); ap.add_argument("refs"); ap.add_argument("out"); ap.add_argument("--lift")
    ap.add_argument("--configs", default="QDNAseq:100kb,Spectre:100kb,BIC-seq2:10kb,CNVkit:30kb,ichorCNA:1Mb low-TF")
    a = ap.parse_args(); R = Path(a.refs); O = Path(a.out); O.mkdir(parents=True, exist_ok=True)
    keep = [tuple(x.split(":")) for x in a.configs.split(",")]
    C = pd.read_csv(a.calls, sep="\t", dtype={"chrom": str})
    C = C[[ (t, c) in keep for t, c in zip(C.tool, C.config)]].copy()
    C["chrom"] = C.chrom.str.replace("chr", "", regex=False)
    C["call_id"] = np.arange(len(C)); C["call_len"] = C.end - C.start; C["lifted_frac"] = 1.0
    if a.lift:
        # liftover piecewise: every aligned block becomes a piece of the call; sequence absent from GRCh37
        # (assembly differences) is reported as the unlifted fraction
        from liftover import Lift
        L = Lift(a.lift); P = []
        for r in C.itertuples(index=False):
            pcs = [p for p in L.pieces("chr" + r.chrom, r.start, r.end) if p[0].replace("chr", "") == "chr" + r.chrom or p[0].replace("chr", "") in AUTO]
            if not pcs: continue
            main = max({p[0] for p in pcs}, key=lambda n: sum(z - y for m, y, z in pcs if m == n))   # majority chromosome
            pcs = [p for p in pcs if p[0] == main]; fr = sum(z - y for _, y, z in pcs) / (r.end - r.start)
            mg = []
            for m, y, z in sorted(pcs, key=lambda p: p[1]):          # merge pieces separated by < 10 kb
                if mg and y - mg[-1][1] < 10_000: mg[-1][1] = max(mg[-1][1], z)
                else: mg.append([y, z])
            for y, z in mg: P.append({**r._asdict(), "chrom": main.replace("chr", ""), "start": y, "end": z, "lifted_frac": fr})
        print(f"liftover: {len({p['call_id'] for p in P})}/{len(C)} calls lifted (piecewise), "
              f"{sum(p['lifted_frac'] * p['call_len'] for p in {q['call_id']: q for q in P}.values()) / C.call_len.sum():.1%} of called bp")
        C = pd.DataFrame(P)
    C = C[C.chrom.isin(AUTO) & (C.end > C.start)]
    C.to_csv(O / "calls_hg19.tsv", sep="\t", index=False)
    ref = load_refs(R); rng = np.random.default_rng(7)

    # per call annotation (pieces summed per call)
    ann = []
    for cid, g in C.groupby("call_id"):
        r = g.iloc[0]; n = (g.end - g.start).sum()
        d = {k: sum(ref["cats"][k].bp(c, s_, e_) for c, s_, e_ in zip(g.chrom, g.start, g.end)) / n for k in ref["cats"]}
        G = ref["genes"]; hit = set()
        for c, s_, e_ in zip(g.chrom, g.start, g.end): hit |= set(G[(G.chrom == c) & (G.start < e_) & (G.end > s_)].gene)
        cgh = sorted(hit & set(ref["cg"].gene))
        ann.append({k: r[k] for k in ("sample_id", "base", "level", "tool", "config", "type", "call_len", "lifted_frac")} |
                   {"chrom": r.chrom, "start": int(g.start.min()), "end": int(g.end.max())} |
                   {f"frac_{k}": round(v, 3) for k, v in d.items()} |
                   {"explained": (d["any"] * r.lifted_frac + (1 - r.lifted_frac)) >= 0.5, "n_genes": len(hit), "cancer_genes": ",".join(cgh[:15])})
    A = pd.DataFrame(ann); A.to_csv(O / "calls_annotated.tsv", sep="\t", index=False)

    # tool-level observed vs size-matched random placement
    res = []
    for (tool, cfg, lvl), g in C.groupby(["tool", "config", "level"]):
        rows = list(zip(g.chrom, g.start, g.end, g.type)); obs = metrics(rows, ref)
        if not obs: continue
        perm = [metrics(shuffle(rows, ref["size"], rng), ref) for _ in range(N_PERM)]
        for k, v in obs.items():
            pv = np.array([p[k] for p in perm if k in p], float)
            res.append({"tool": tool, "config": cfg, "level": lvl, "metric": k, "observed": v, "expected": np.nanmean(pv),
                        "enrichment": v / np.nanmean(pv) if np.nanmean(pv) > 0 else np.nan,
                        "p_perm": (1 + np.sum(pv >= v)) / (1 + len(pv)), "n_calls": g.call_id.nunique(),
                        "unlifted_bp_frac": 1 - (g.end - g.start).sum() / g.drop_duplicates("call_id").call_len.sum(), "called_Mb": sum(e - s for _, s, e, _ in rows) / 1e6})
    pd.DataFrame(res).to_csv(O / "tool_annotation_vs_random.tsv", sep="\t", index=False)

    # per-sample TCGA score (for the cancer cohort and its healthy baseline)
    ps = []
    for (sid, tool, cfg, lvl, base), g in C.groupby(["sample_id", "tool", "config", "level", "base"]):
        rows = list(zip(g.chrom, g.start, g.end, g.type))
        ps.append({"sample_id": sid, "base": base, "tool": tool, "config": cfg, "level": lvl, "called_Mb": sum(e - s for _, s, e, _ in rows) / 1e6,
                   **{f"tcga_{k}": tcga_score(rows, f) for k, f in ref["freq"].items()}})
    pd.DataFrame(ps).to_csv(O / "sample_tcga_scores.tsv", sep="\t", index=False)

    # recurrent hotspots (full depth): 1-Mb windows altered in >= half of the people, with likely cause
    F = C[C.level == "full"]; nb = F.base.nunique(); win = defaultdict(lambda: defaultdict(set))
    for r in F.itertuples():
        for w in range(r.start // W, (r.end - 1) // W + 1): win[(r.chrom, w * W)][(r.tool, r.type)].add(r.base)
    hs = []
    for (c, s), d in win.items():
        tools = {f"{t} {ty}": len(v) for (t, ty), v in d.items() if len(v) >= nb / 2}
        if not tools: continue
        fr = {k: ref["cats"][k].bp(c, s, s + W) / W for k in ref["cats"]}
        g = ref["genes"]; hit = g[(g.chrom == c) & (g.start < s + W) & (g.end > s)].gene
        cause = max(("germline_CNV", "segdup", "blacklist", "pericentromeric"), key=lambda k: fr[k]) if fr["any"] >= 0.1 else "unexplained"
        hs.append({"chrom": c, "start": s, "end": s + W, "tools(n_people)": "; ".join(f"{k} ({v})" for k, v in sorted(tools.items())),
                   "n_tools": len({k.split()[0] for k in tools}), **{f"frac_{k}": round(v, 2) for k, v in fr.items()},
                   "likely_cause": cause, "n_genes": len(hit), "cancer_genes": ",".join(sorted(set(hit) & set(ref["cg"].gene)))})
    H = pd.DataFrame(hs)
    if len(H): H = H.sort_values(["n_tools", "chrom", "start"], ascending=[False, True, True], key=lambda x: x.map(lambda v: int(v) if isinstance(v, str) and v.isdigit() else v))
    H.to_csv(O / "recurrent_hotspots_masklist.tsv", sep="\t", index=False)
    print("done:", len(A), "calls;", len(H), "hotspot windows")

if __name__ == "__main__":
    main()
