#!/usr/bin/env python3
"""
D - Simulation 2 (HG008) scoring.

HG008-T has a mean ploidy of ~1.44 (most of its genome is one copy). Read-depth callers
normalise to the sample's own median, so the one-copy state is their "neutral" baseline.
Scoring gain/loss against two copies would therefore be wrong. Three readouts:

1. TUMOUR FRACTION  ichorCNA (3 configs) and ACE (ploidy 2 / free ploidy, 3 bins):
   estimate vs the measured (realised) tumour fraction of each mixture.
2. PROFILE CORRELATION (primary, no thresholds, ploidy-free)
   1-Mb windows inside the NIST benchmark (>= 80 % covered).
   expected shape d_w = log2(CN_w / P) (CN_w = benchmark copy number, P = mean ploidy);
   observed o_w = length-weighted mean of the tool's segment values in the window
   (log2 ratio where the tool reports one; +1 / -1 / 0 for class-only tools).
   r = Pearson(o, d). The 0 % mixtures (no tumour DNA) give the null distribution;
   limit of detection (LoD) at each depth = lowest tumour fraction whose median r exceeds the
   highest r seen in any 0 % mixture of that tool/configuration.
3. EVENT LEVEL (secondary)  truth events = benchmark segments whose copy number differs from
   the modal copy number (the callers' baseline): above -> gain, below -> loss; >= 50 kb.
   Calls are clipped to benchmark regions; one-to-one matching, same direction, RO >= 0.5
   (same rule as Simulation 1). A diploid-baseline variant is reported for transparency.

  python d_score_sim2.py            (paths from s2_env.sh)
Outputs in 09_metrics/simulation2/: tf_estimates.csv, profile_correlation.csv,
lod_by_depth.csv, event_metrics.csv, summary_*.csv, SCORING_SUMMARY.txt
"""
import csv, gzip, importlib.util, os, re
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

P = Path(os.environ["P"]); OUT = Path(os.environ["S2_OUT"]); MET = Path(os.environ["S2_MET"]); MET.mkdir(parents=True, exist_ok=True)
spec = importlib.util.spec_from_file_location("sc", P / "00_metadata/rescue/50_score_sim1_v2.py")
sc = importlib.util.module_from_spec(spec); spec.loader.exec_module(sc)
AUTO = [str(i) for i in range(1, 23)]
BINS = [10000, 30000, 50000, 100000, 500000, 1000000]
nc = lambda c: re.sub(r"^chr", "", str(c))


# ---------------------------------------------------------------- truth
def load_truth(path):
    t = pd.read_csv(path, sep="\t", header=None, comment="#", dtype={0: str})
    t = pd.DataFrame({"chrom": t[0].map(nc), "start": t[1].astype(int), "end": t[2].astype(int), "cn": t[3].astype(float)})
    t = t[t.chrom.isin(AUTO) & (t.end > t.start)].sort_values(["chrom", "start"])
    L = (t.end - t.start)
    P_ = float((t.cn * L).sum() / L.sum())
    mode = float(t.groupby("cn").apply(lambda g: (g.end - g.start).sum()).idxmax())
    return t, P_, mode


def truth_events(t, base, minlen=50000):
    ev = []
    for c, g in t.groupby("chrom"):
        cur = None
        for r in g.itertuples():
            st = "gain" if r.cn > base else ("loss" if r.cn < base else None)
            if cur and st == cur["type"] and r.start - cur["end"] < 1000:
                cur["end"] = r.end
            else:
                if cur: ev.append(cur)
                cur = {"chrom": c, "start": r.start, "end": r.end, "type": st} if st else None
        if cur: ev.append(cur)
    return [e for e in ev if e["end"] - e["start"] >= minlen]


def grid(t, w=1_000_000, minfrac=0.8):
    rows = []
    for c, g in t.groupby("chrom"):
        s, e, cn = g.start.values, g.end.values, g.cn.values
        for ws in range(0, int(e.max()), w):
            we = ws + w; ov = np.clip(np.minimum(e, we) - np.maximum(s, ws), 0, None)
            if ov.sum() >= minfrac * w:
                rows.append((c, ws, we, float((cn * ov).sum() / ov.sum())))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "cn"])


def window_mean(G, segs):
    """length-weighted mean of segment values per grid window (0 where nothing is reported)"""
    o = np.zeros(len(G))
    if not segs: return o
    S = pd.DataFrame(segs, columns=["c", "s", "e", "v"]); S = S[np.isfinite(S.v)]
    for c, gi in G.groupby("chrom").groups.items():
        x = S[S.c == c]
        if not len(x): continue
        ws, we = G.loc[gi, "start"].values[:, None], G.loc[gi, "end"].values[:, None]
        ov = np.clip(np.minimum(x.e.values[None, :], we) - np.maximum(x.s.values[None, :], ws), 0, None)
        o[np.asarray([G.index.get_loc(i) for i in gi])] = (ov * x.v.values[None, :]).sum(1) / (we[:, 0] - ws[:, 0])
    return o


def clip(calls, t):
    out = []; by = defaultdict(list)
    for r in t.itertuples(): by[r.chrom].append((r.start, r.end))
    for c in calls:
        for s, e in by.get(c["chrom"], []):
            a, b = max(s, c["start"]), min(e, c["end"])
            if b > a: out.append({**c, "start": a, "end": b})
    merged = []
    for c in sorted(out, key=lambda x: (x["chrom"], x["start"])):
        if merged and merged[-1]["chrom"] == c["chrom"] and merged[-1]["type"] == c["type"] and c["start"] - merged[-1]["end"] < 1000:
            merged[-1]["end"] = max(merged[-1]["end"], c["end"])
        else: merged.append(dict(c))
    return [c for c in merged if c["end"] - c["start"] >= 50000]


# ---------------------------------------------------------------- readers -> (segments with values, calls)
def rd(path, **kw):
    return pd.read_csv(path, sep="\t", **kw)


def segs_from(df, c, s, e, v):
    return [(nc(a), int(b), int(x), float(y)) for a, b, x, y in zip(df[c], df[s], df[e], df[v]) if nc(a) in AUTO]


def calls_from_states(segs, thr):
    return [{"chrom": c, "start": s, "end": e, "type": "gain" if v >= thr else "loss"} for c, s, e, v in segs if abs(v) >= thr]


def calls_vcf(path):
    try: return [dict(x, chrom=nc(x["chrom"])) for x in sc.parse_bedlike(Path(path), explicit=True)]
    except Exception: return []


def classfile(path):
    segs, calls = [], []
    op = gzip.open if str(path).endswith(".gz") else open
    for ln in op(path, "rt"):
        p = ln.rstrip("\n").split("\t")
        if len(p) < 3 or not p[1].isdigit(): continue
        txt = " ".join(p[3:]).lower()
        t = "gain" if ("dup" in txt or "gain" in txt) else ("loss" if ("del" in txt or "loss" in txt) else None)
        if t and nc(p[0]) in AUTO:
            segs.append((nc(p[0]), int(p[1]), int(p[2]), 1.0 if t == "gain" else -1.0))
            calls.append({"chrom": nc(p[0]), "start": int(p[1]), "end": int(p[2]), "type": t})
    return segs, calls


def lyceum_events(path):
    d = pd.read_csv(path, sep="\t"); d.columns = ["s", "c", "st", "en", "p"][:len(d.columns)]
    d["k"] = np.where(d.p.str.contains("DEL"), 2, np.where(d.p.str.contains("DUP"), 1, 0)); segs, calls = [], []
    for c, g in d.groupby("c"):
        g = g.sort_values("st"); k, s, e = g.k.values, g.st.values, g.en.values; i = 0
        while i < len(k):
            if k[i] == 0: i += 1; continue
            j, last, gap, n = i, i, 0, 1
            while j + 1 < len(k):
                j += 1
                if k[j] == k[i]: last, gap, n = j, 0, n + 1
                elif k[j] == 0 and gap < 2: gap += 1
                else: break
            if n >= 2 and nc(c) in AUTO:
                t = "gain" if k[i] == 1 else "loss"
                segs.append((nc(c), int(s[i]), int(e[last]), 1.0 if t == "gain" else -1.0))
                calls.append({"chrom": nc(c), "start": int(s[i]), "end": int(e[last]), "type": t})
            i = last + 1
    return segs, calls


def first(pattern_dir, *globs):
    for g in globs:
        f = sorted(Path(pattern_dir).glob(g))
        if f: return f[0]
    return None


FAIL = []   # parse problems are logged here instead of stopping the whole scoring


def _pick(cols, *cands, contains=None, exclude=None):
    low = {c.lower(): c for c in cols}
    for c in cands:
        if c.lower() in low: return low[c.lower()]
    if contains:
        for c in cols:
            if contains.lower() in c.lower() and not (exclude and exclude.lower() in c.lower()): return c
    return None


def _qdnaseq(s, b, d):
    g = rd(d / f"{s}.segments.igv", comment="#"); vc = g.columns[4]
    return segs_from(g, g.columns[0], g.columns[1], g.columns[2], vc), calls_vcf(d / f"{s}.calls.vcf")


def _cnvkit(s, b, d):
    cns = sorted(f for f in d.glob("*.cns") if not f.name.endswith((".call.cns", ".bintest.cns")))
    cf = first(d, "*.call.cns")
    if not cns and cf is None: raise FileNotFoundError(f"no .cns in {d}")
    g = rd(cns[0] if cns else cf); cc = rd(cf) if cf is not None else None
    calls = []
    if cc is not None and "cn" in cc.columns:
        calls = [{"chrom": nc(a), "start": int(x), "end": int(y), "type": "gain" if n > 2 else "loss"}
                 for a, x, y, n in zip(cc.chromosome, cc.start, cc.end, cc.cn) if n != 2 and nc(a) in AUTO]
    sg = segs_from(g, "chromosome", "start", "end", "log2")
    return sg, (calls if cc is not None else calls_from_states(sg, 0.2))


def _spectre(s, b, d):
    f = first(d / "calls", "*.vcf.gz", "*.vcf"); calls = calls_vcf(f) if f else []
    return [(c["chrom"], c["start"], c["end"], 1.0 if c["type"] == "gain" else -1.0) for c in calls], calls


def _bicseq2(s, b, d):
    g = rd(d / "seg" / f"{s}.cnv")
    lc = _pick(g.columns, "log2.copyRatio", contains="log2")
    if lc is None: raise KeyError(f"no log2 column in {list(g.columns)}")
    sg = segs_from(g, g.columns[0], g.columns[1], g.columns[2], lc)
    return sg, calls_from_states(sg, 0.2)


def _ichor(s, cfg, d):
    f = first(d, f"{s}.seg.txt", "*.seg.txt", f"{s}.seg", "*.seg")
    if f is None: raise FileNotFoundError(f"no .seg in {d}")
    g = rd(f)
    ch = _pick(g.columns, "chrom", "chr", "chromosome"); st = _pick(g.columns, "start", "loc.start"); en = _pick(g.columns, "end", "loc.end")
    vc = _pick(g.columns, "seg.median.logR", "median", "seg.mean", contains="logr", exclude="copy")
    if vc is None: raise KeyError(f"no logR column in {f.name}: {list(g.columns)}")
    sg = segs_from(g, ch, st, en, vc)
    cn = _pick(g.columns, "Corrected_Copy_Number", "copy.number", "Copy_Number")
    if cn is not None:
        calls = [{"chrom": nc(a), "start": int(x), "end": int(y), "type": "gain" if n > 2 else "loss"}
                 for a, x, y, n in zip(g[ch], g[st], g[en], g[cn]) if pd.notna(n) and n != 2 and nc(a) in AUTO]
    else:
        calls = calls_from_states(sg, 0.2)
    return sg, calls


def tool_outputs(s):
    """yield (tool, config, bin, segments, calls) for one mixture; a broken output is logged in FAIL"""
    jobs = []
    for b in BINS:
        d = OUT / "QDNAseq" / s / str(b)
        if (d / ".done").exists(): jobs.append(("QDNAseq", "default", b, lambda s=s, b=b, d=d: _qdnaseq(s, b, d)))
        d = OUT / "CNVkit" / s / str(b)
        if (d / ".done").exists(): jobs.append(("CNVkit", "flat", b, lambda s=s, b=b, d=d: _cnvkit(s, b, d)))
        d = OUT / "Spectre" / s / str(b)
        if (d / ".done").exists(): jobs.append(("Spectre", "cancer", b, lambda s=s, b=b, d=d: _spectre(s, b, d)))
        d = OUT / "BICseq2" / s / f"bin_{b}"
        if (d / ".done").exists(): jobs.append(("BICseq2", "lambda2", b, lambda s=s, b=b, d=d: _bicseq2(s, b, d)))
        d = OUT / "dudeML" / "predictions" / s / f"bin_{b}"
        if (d / ".done").exists(): jobs.append(("dudeML", "sim1_models", b, lambda s=s, d=d: classfile(d / f"{s}.segments.bed")))
        f = OUT / "CNV-MEANN" / "calls" / f"bin_{b}" / f"{s}.events.tsv"
        if f.exists(): jobs.append(("CNV-MEANN", "v2_sim1_trained", b, lambda f=f: classfile(f)))
    for cfg in ("1000kb_default", "1000kb_lowTF", "500kb_default"):
        d = OUT / "ichorCNA" / s / cfg
        if (d / ".done").exists():
            jobs.append(("ichorCNA", cfg, 1000000 if cfg.startswith("1000") else 500000, lambda s=s, cfg=cfg, d=d: _ichor(s, cfg, d)))
    f = OUT / "LYCEUM" / s / "calls" / f"{s}_exonlevel.csv"
    if f.exists(): jobs.append(("LYCEUM", "native_hg38", 0, lambda f=f: lyceum_events(f)))
    for tool, cfg, b, fn in jobs:
        try:
            sg, cl = fn()
        except Exception as e:
            FAIL.append({"sample": s, "tool": tool, "config": cfg, "bin_bp": b, "error": f"{type(e).__name__}: {e}"[:300]})
            continue
        yield tool, cfg, b, sg, cl


def tf_rows(s):
    for cfg in ("1000kb_default", "1000kb_lowTF", "500kb_default"):
        f = OUT / "ichorCNA" / s / cfg / f"{s}.params.txt"
        if f.exists():
            m = re.search(r"Tumor Fraction:\s*([0-9.eE-]+)", f.read_text())
            if m: yield "ichorCNA", cfg, float(m.group(1))
    for b in (100000, 500000, 1000000):
        f = OUT / "ACE" / s / str(b) / f"{s}.ace.tsv"
        if f.exists():
            try:
                a = pd.read_csv(f, sep="\t").iloc[0]
                yield "ACE", f"ploidy2_{b // 1000}kb", float(a.ace_cellularity_ploidy2)
                if pd.notna(a.get("ace_cellularity_free")): yield "ACE", f"freeploidy_{b // 1000}kb", float(a.ace_cellularity_free)
            except Exception as e:
                FAIL.append({"sample": s, "tool": "ACE", "config": str(b), "bin_bp": b, "error": f"{type(e).__name__}: {e}"[:300]})


# ---------------------------------------------------------------- main
def main():
    t, ploidy, mode = load_truth(os.environ["TRUTH38"])
    G = grid(t); d_shape = np.log2(np.clip(G.cn.values, 0.05, None) / ploidy)
    TE = {"modal": truth_events(t, mode), "diploid": truth_events(t, 2.0)}
    man = pd.read_csv(Path(os.environ["S2_META"]) / "HG008_mixtures_realised.tsv", sep="\t")
    print(f"truth: {len(t)} segments, {(t.end - t.start).sum() / 1e6:.0f} Mb, mean ploidy {ploidy:.2f}, modal CN {mode:g}, "
          f"{len(G)} 1-Mb windows, events vs modal {len(TE['modal'])}, vs diploid {len(TE['diploid'])}", flush=True)
    tf, prof, evm = [], [], []
    for m in man.itertuples():
        s, f, dep = m.sample_id, float(m.realised_fraction), float(m.depth_x)
        base = {"sample": s, "depth_x": dep, "nominal_fraction": m.tumor_fraction, "realised_fraction": f, "replicate": m.replicate}
        for tool, cfg, est in tf_rows(s):
            tf.append({**base, "tool": tool, "config": cfg, "estimated_fraction": est, "error": est - f})
        for tool, cfg, b, segs, calls in tool_outputs(s):
          try:
            o = window_mean(G, segs)
            r = float(np.corrcoef(o, d_shape)[0, 1]) if o.std() > 0 else 0.0
            slope = float(np.polyfit(d_shape, o, 1)[0]) if o.std() > 0 else 0.0
            prof.append({**base, "tool": tool, "config": cfg, "bin_bp": b, "pearson_r": r, "slope": slope,
                         "windows_with_signal": int((o != 0).sum())})
            cl = clip(calls, t)
            for name, T in TE.items():
                truth = T if f > 0 else []
                mm = sc.match(truth, cl, .5, True); tp = len(mm); fp = len(cl) - tp; fn = len(truth) - tp
                p, rc = sc.div(tp, tp + fp), sc.div(tp, tp + fn)
                evm.append({**base, "tool": tool, "config": cfg, "bin_bp": b, "baseline": name, "n_truth": len(truth),
                            "n_calls": len(cl), "TP": tp, "FP": fp, "FN": fn, "precision": p, "recall": rc, "F1": sc.f1(p, rc)})
          except Exception as e:
            FAIL.append({"sample": s, "tool": tool, "config": cfg, "bin_bp": b, "error": f"scoring {type(e).__name__}: {e}"[:300]})
        print("scored", s, flush=True)
    fl = pd.DataFrame(FAIL, columns=["sample", "tool", "config", "bin_bp", "error"]); fl.to_csv(MET / "parse_failures.csv", index=False)
    fl["error"] = fl["error"].str.replace(r"/\S+/", "", regex=True).str.slice(0, 120)
    print(f"parse/scoring failures: {len(fl)}" + ("" if fl.empty else "\n" + fl.groupby(["tool", "error"]).size().head(20).to_string()), flush=True)
    tf, prof, evm = pd.DataFrame(tf), pd.DataFrame(prof), pd.DataFrame(evm)
    tf.to_csv(MET / "tf_estimates.csv", index=False); prof.to_csv(MET / "profile_correlation.csv", index=False)
    evm.to_csv(MET / "event_metrics.csv", index=False)

    # limit of detection from the profile correlation
    lod = []
    for (tool, cfg, b), g in prof.groupby(["tool", "config", "bin_bp"]):
        null = g[g.nominal_fraction == 0].pearson_r; thr = null.max() if len(null) else np.nan
        for dep, h in g[g.nominal_fraction > 0].groupby("depth_x"):
            med = h.groupby("nominal_fraction").pearson_r.median()
            ok = med[med > thr]
            lod.append({"tool": tool, "config": cfg, "bin_bp": b, "depth_x": dep, "null_max_r": thr,
                        "LoD_fraction": float(ok.index.min()) if len(ok) else np.nan,
                        **{f"median_r_f{int(round(k * 100)):02d}": v for k, v in med.items()}})
    lod = pd.DataFrame(lod); lod.to_csv(MET / "lod_by_depth.csv", index=False)

    s1 = prof.groupby(["tool", "config", "bin_bp", "nominal_fraction"]).pearson_r.median().unstack().round(3)
    s1.to_csv(MET / "summary_profile_r.csv")
    s2 = evm[evm.baseline == "modal"].groupby(["tool", "config", "bin_bp", "nominal_fraction"]).F1.median().unstack().round(3)
    s2.to_csv(MET / "summary_event_F1_modal.csv")
    if len(tf):
        s3 = tf.groupby(["tool", "config", "nominal_fraction"]).estimated_fraction.median().unstack().round(3)
        s3.to_csv(MET / "summary_tf_estimates.csv")
    with open(MET / "SCORING_SUMMARY.txt", "w") as fh:
        fh.write(f"mixtures: {len(man)}\nmean ploidy {ploidy:.3f}, modal CN {mode:g}, 1-Mb windows {len(G)}\n"
                 f"truth events vs modal {len(TE['modal'])}, vs diploid {len(TE['diploid'])}\n\n")
        fh.write("tool x config x bin: number of mixtures scored\n")
        fh.write(prof.groupby(["tool", "config", "bin_bp"]).size().to_string() + "\n\n")
        fh.write("median profile correlation by tumour fraction (all depths)\n" + s1.to_string() + "\n")
        fh.write(f"\nparse/scoring failures: {len(fl)} (details in parse_failures.csv)\n")
        if len(fl): fh.write(fl.groupby(["tool", "error"]).size().to_string() + "\n")
    print(open(MET / "SCORING_SUMMARY.txt").read())


if __name__ == "__main__":
    main()
