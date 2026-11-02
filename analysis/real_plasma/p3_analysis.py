#!/usr/bin/env python3
"""Phase 3 (PRJNA795275 WGS plasma, all labelled CTR = non-cancer controls): per-sample call burden,
fraction of genome altered (FGA), tumour-fraction estimates, cross-caller concordance, depth stability."""
import sys, re, glob, gzip
from pathlib import Path
import numpy as np, pandas as pd
R = Path(sys.argv[1]); O = Path(sys.argv[2]); O.mkdir(parents=True, exist_ok=True)
OUT = R / "07_tool_outputs/phase3"; META = R / "00_metadata/phase3"
AUTO = {f"chr{i}" for i in range(1, 23)} | {str(i) for i in range(1, 23)}
GENOME = 2.875e9  # GRCh38 autosomes
nc = lambda c: c if str(c).startswith("chr") else f"chr{c}"

def qdnaseq(f):
    out = []
    for ln in open(f):
        if ln.startswith("#") or ln.count("\t") < 7: continue
        p = ln.split("\t"); info = dict(x.split("=") for x in p[7].split(";") if "=" in x)
        if nc(p[0]) in AUTO: out.append((nc(p[0]), int(p[1]), int(info["END"]), "gain" if info["SVTYPE"] == "DUP" else "loss"))
    return out
def cnvkit(f):
    d = pd.read_csv(f, sep="\t")
    return [(nc(c), s, e, "gain" if n > 2 else "loss") for c, s, e, n in zip(d.chromosome, d.start, d.end, d.cn) if nc(c) in AUTO and n != 2]
def ichor(f):
    d = pd.read_csv(f, sep="\t")
    return [(nc(c), s, e, "gain" if n > 2 else "loss") for c, s, e, n in zip(d.chrom, d.start, d.end, d.Corrected_Copy_Number) if nc(c) in AUTO and pd.notna(n) and n != 2]
def spectre(f):
    out = []
    with gzip.open(f, "rt") as fh:
        for ln in fh:
            p = ln.rstrip().split("\t")
            if len(p) >= 4 and nc(p[0]) in AUTO: out.append((nc(p[0]), int(p[1]), int(p[2]), "gain" if p[3].upper().startswith("DUP") else "loss"))
    return out
CFG = {("QDNAseq", "50kb"): ("QDNAseq/{s}/50000/{s}.calls.vcf", qdnaseq), ("QDNAseq", "100kb"): ("QDNAseq/{s}/100000/{s}.calls.vcf", qdnaseq),
       ("CNVkit", "30kb"): ("CNVkit/{s}/30000/*.call.cns", cnvkit), ("CNVkit", "100kb"): ("CNVkit/{s}/100000/*.call.cns", cnvkit),
       ("ichorCNA", "1Mb low-TF"): ("ichorCNA/{s}/1000kb_lowTF/{s}.seg.txt", ichor), ("ichorCNA", "1Mb defaults"): ("ichorCNA/{s}/1000kb_default/{s}.seg.txt", ichor),
       ("ichorCNA", "500kb defaults"): ("ichorCNA/{s}/500kb_default/{s}.seg.txt", ichor),
       ("Spectre", "100kb"): ("Spectre/{s}/100000/calls/{s}_cnv.bed.gz", spectre), ("Spectre", "500kb"): ("Spectre/{s}/500000/calls/{s}_cnv.bed.gz", spectre),
       ("BIC-seq2", "10kb"): ("BICseq2/{s}/bin_10000/seg/{s}.cnv", None), ("BIC-seq2", "30kb"): ("BICseq2/{s}/bin_30000/seg/{s}.cnv", None)}
def bic(f, thr=0.2):
    d = pd.read_csv(f, sep="\t"); lc = [c for c in d.columns if "log2" in c.lower()][0]
    return [(nc(c), s, e, "gain" if v >= thr else "loss") for c, s, e, v in zip(d.iloc[:, 0], d.iloc[:, 1], d.iloc[:, 2], d[lc]) if nc(c) in AUTO and abs(v) >= thr]

bams = pd.read_csv(META / "p3_bams.tsv", sep="\t")
rows, calls = [], {}
for b in bams.itertuples():
    for (tool, cfg), (pat, fn) in CFG.items():
        f = glob.glob(str(OUT / pat.format(s=b.sample_id)))
        if not f: continue
        c = (fn or bic)(f[0]); c = [x for x in c if x[2] > x[1]]
        calls[(b.sample_id, tool, cfg)] = c
        L = np.array([e - s for _, s, e, _ in c]) if c else np.array([])
        rows.append({"sample_id": b.sample_id, "base": b.base_sample, "level": b.level, "depth_x": b.depth_x, "tool": tool, "config": cfg,
                     "n_calls": len(c), "n_gain": sum(t == "gain" for *_, t in c), "n_loss": sum(t == "loss" for *_, t in c),
                     "altered_Mb": L.sum() / 1e6, "FGA_pct": 100 * L.sum() / GENOME, "median_call_kb": np.median(L) / 1e3 if len(L) else np.nan,
                     "calls_ge_1Mb": int((L >= 1e6).sum()) if len(L) else 0})
B = pd.DataFrame(rows); B.to_csv(O / "p3_burden.tsv", sep="\t", index=False)
lv = dict(zip(bams.sample_id, bams.level)); bs = dict(zip(bams.sample_id, bams.base_sample))
pd.DataFrame([{"sample_id": k[0], "base": bs[k[0]], "level": lv[k[0]], "tool": k[1], "config": k[2], "chrom": ch, "start": st, "end": e, "type": t}
              for k, c in calls.items() for ch, st, e, t in c]).to_csv(O / "p3_calls_all.tsv", sep="\t", index=False)

# tumour fraction estimates
tf = []
for b in bams.itertuples():
    for cfg in ("1000kb_lowTF", "1000kb_default", "500kb_default"):
        f = OUT / f"ichorCNA/{b.sample_id}/{cfg}/{b.sample_id}.params.txt"
        if f.exists():
            t = f.read_text(); m = re.search(r"Tumor Fraction:\s*([0-9.eE-]+)", t); p = re.search(r"Ploidy:\s*([0-9.eE-]+)", t)
            tf.append({"sample_id": b.sample_id, "base": b.base_sample, "level": b.level, "method": f"ichorCNA {cfg}", "TF": float(m.group(1)), "ploidy": float(p.group(1))})
    for bin_ in (100000, 500000, 1000000):
        f = OUT / f"ACE/{b.sample_id}/{bin_}/{b.sample_id}.ace.tsv"
        if f.exists():
            a = pd.read_csv(f, sep="\t").iloc[0]
            tf.append({"sample_id": b.sample_id, "base": b.base_sample, "level": b.level, "method": f"ACE ploidy2 {bin_//1000}kb", "TF": a.ace_cellularity_ploidy2, "ploidy": 2})
    f = OUT / f"Fragle/{b.sample_id}/out/Fragle.csv"
    if f.exists():
        tf.append({"sample_id": b.sample_id, "base": b.base_sample, "level": b.level, "method": "Fragle", "TF": float(pd.read_csv(f).iloc[0, 1]), "ploidy": np.nan})
TF = pd.DataFrame(tf); TF.to_csv(O / "p3_tf.tsv", sep="\t", index=False)

# cross-caller concordance on full-depth samples (bp Jaccard of altered sequence, direction-aware)
def iv(c):
    d = {}
    for ch, s, e, t in c: d.setdefault((ch, t), []).append((s, e))
    return d
def ovl(a, b):
    tot = 0
    for k, A in a.items():
        for s1, e1 in A:
            for s2, e2 in b.get(k, []):
                tot += max(0, min(e1, e2) - max(s1, s2))
    return tot
prim = [("QDNAseq", "100kb"), ("CNVkit", "30kb"), ("ichorCNA", "1Mb low-TF"), ("Spectre", "100kb"), ("BIC-seq2", "10kb")]
conc = []
for s in bams[bams.level == "full"].sample_id:
    for i, x in enumerate(prim):
        for y in prim[i + 1:]:
            if (s, *x) in calls and (s, *y) in calls:
                a, b_ = calls[(s, *x)], calls[(s, *y)]
                la = sum(e - st for _, st, e, _ in a); lb = sum(e - st for _, st, e, _ in b_); o = ovl(iv(a), iv(b_))
                conc.append({"sample_id": s, "A": f"{x[0]} {x[1]}", "B": f"{y[0]} {y[1]}", "A_Mb": la / 1e6, "B_Mb": lb / 1e6, "overlap_Mb": o / 1e6,
                             "jaccard": o / (la + lb - o) if (la + lb - o) > 0 else np.nan})
pd.DataFrame(conc).to_csv(O / "p3_concordance.tsv", sep="\t", index=False)

# recurrent regions across samples (full depth), 1-Mb windows
rec = []
for (tool, cfg) in CFG:
    win = {}
    for s in bams[bams.level == "full"].sample_id:
        seen = set()
        for ch, st, e, t in calls.get((s, tool, cfg), []):
            for w in range(st // 1_000_000, (e - 1) // 1_000_000 + 1): seen.add((ch, w, t))
        for k in seen: win[k] = win.get(k, 0) + 1
    for (ch, w, t), n in win.items():
        rec.append({"tool": tool, "config": cfg, "chrom": ch, "start": w * 1_000_000, "end": (w + 1) * 1_000_000, "type": t, "n_samples": n})
pd.DataFrame(rec).to_csv(O / "p3_recurrent_windows.tsv", sep="\t", index=False)

# insert-size summary
ins = []
for f in glob.glob(str(R / "09_metrics/phase3/qc/*.insert_sizes.tsv")):
    d = pd.read_csv(f, sep="\t", header=None, names=["size", "n"])
    if len(d) == 0 or d.n.sum() == 0: continue
    d = d[d["size"] <= 1000]; cum = d.n.cumsum() / d.n.sum()
    ins.append({"base": Path(f).name.split(".")[0], "mode_bp": int(d.loc[d.n.idxmax(), "size"]), "median_bp": int(d["size"][cum >= 0.5].iloc[0]),
                "frac_lt150": d.n[d["size"] < 150].sum() / d.n.sum(), "frac_gt250": d.n[d["size"] > 250].sum() / d.n.sum()})
pd.DataFrame(ins).to_csv(O / "p3_insert_sizes.tsv", sep="\t", index=False)
print("done")
