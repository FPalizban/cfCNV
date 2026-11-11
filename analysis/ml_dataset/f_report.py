#!/usr/bin/env python
"""F  Manuscript tables: the two baseline models (v1 = Simulation 1 only, v2 = augmented) next to the 8 callers,
scored with EXACTLY the manuscript's definitions (taken from d_score_sim2.py / figures/common.py):

HG008  profile r: 1-Mb windows with >= 80 % NIST benchmark coverage; truth d = log2(max(CN,0.05) / mean ploidy);
       observed = mean of the model value over the window (0 where nothing is reported); r = Pearson(o, d), 0 if o constant.
       model value: P(gain) - P(loss) for the gradient-boosting models, the call (+1/0/-1) for the threshold models,
       the GC-corrected log2 ratio for "signal (no model)".
       strict LoD (figures/common.py lod_strict): threshold = max(null max, null mean + 3 SD) over the 0 % mixtures;
       LoD at a depth = lowest tumour fraction at which every replicate exceeds it; not detected = 0.5 for the mean.
       median_r_1to10 = median r over all mixtures with 1-10 % tumour fraction.
Cancer F1 against the 5-caller consensus on the 2,897 GRCh37 1-Mb windows (window altered when >= 50 % of its 100-kb
       bins are called in one direction); callers: their leave-one-out F1 from cancer_loo_and_tcga.tsv (native depth).
Non-cancer % of the genome called altered (100-kb bins) next to the callers' primary configurations.
Output: $MLO/manuscript/  (tables *.tsv, key_numbers.md, methods_notes.md)
"""
import os, re, json
from pathlib import Path
import numpy as np, pandas as pd

E = os.environ
MLO = Path(E["MLO"]); ML = Path(E["ML"]); DD = ML / "data"
OUT = MLO / "manuscript"; OUT.mkdir(parents=True, exist_ok=True)
DS = MLO / "hf_dataset" / "data" / "bins_100kb"
R1, R2 = MLO / "results", MLO / "results_v2"
AUTO = [str(i) for i in range(1, 23)]
nc = lambda c: re.sub(r"^chr", "", str(c))
NAME = {"BICseq2": "BIC-seq2"}
HEALTHY_PRIMARY = [("QDNAseq", "100kb"), ("Spectre", "100kb"), ("BIC-seq2", "10kb"), ("CNVkit", "30kb"), ("ichorCNA", "1Mb low-TF")]
MODELS = {"ML v1 gradient boosting": (R1, "hgb"), "ML v1 threshold": (R1, "thr"),
          "ML v2 gradient boosting": (R2, "hgb"), "ML v2 threshold": (R2, "thr")}
note = []


def preds(rdir, split, cols):
    f = rdir / "predictions" / f"{split}.csv.gz"
    return pd.read_csv(f, dtype={"chrom": str}, usecols=lambda c: c in cols) if f.exists() else None


def meta(split):
    m = pd.read_csv(DS / f"{split}.csv.gz", dtype={"chrom": str, "sample_id": str, "base_sample": str, "replicate": str},
                    usecols=["sample_id", "base_sample", "replicate", "depth_x", "tumour_fraction_nominal", "chrom", "start", "lr"])
    return m


def value_col(model, d):
    kind = MODELS[model][1]
    return (d.p_gain - d.p_loss) if kind == "hgb" else d.pred_thr


# ================================================================ HG008 profile correlation + LoD
t = pd.read_csv(E["TRUTH38"], sep="\t", header=None, comment="#", dtype={0: str})
t = pd.DataFrame({"chrom": t[0].map(nc), "start": t[1].astype(int), "end": t[2].astype(int), "cn": t[3].astype(float)})
t = t[t.chrom.isin(AUTO) & (t.end > t.start)].sort_values(["chrom", "start"])
L = t.end - t.start; PLOIDY = float((t.cn * L).sum() / L.sum())
G = []
for c, g in t.groupby("chrom"):
    s, e, cn = g.start.values, g.end.values, g.cn.values
    for ws in range(0, int(e.max()), 1_000_000):
        ov = np.clip(np.minimum(e, ws + 1_000_000) - np.maximum(s, ws), 0, None)
        if ov.sum() >= 0.8 * 1_000_000:
            G.append((c, ws, float((cn * ov).sum() / ov.sum())))
G = pd.DataFrame(G, columns=["chrom", "win", "cn"]); G["d"] = np.log2(np.clip(G.cn, 0.05, None) / PLOIDY)
G["key"] = G.chrom + ":" + G.win.astype(str)
print(f"HG008 grid: {len(G)} 1-Mb windows, mean ploidy {PLOIDY:.3f}", flush=True)


def profile_r(d, val):
    """d: one sample's 100-kb bins; val: values per bin -> Pearson r on the 1-Mb grid (manuscript definition)"""
    w = (d.chrom + ":" + ((d.start // 1_000_000) * 1_000_000).astype(str)).values
    s = pd.Series(np.asarray(val, float) * 0.1, index=w).groupby(level=0).sum()     # 100 kb / 1 Mb = 0.1
    o = s.reindex(G.key).fillna(0).values
    return float(np.corrcoef(o, G.d.values)[0, 1]) if o.std() > 0 else 0.0


hm = meta("test_hg008")
rows = []
for model, (rdir, kind) in MODELS.items():
    p = preds(rdir, "test_hg008", {"sample_id", "chrom", "start", "p_loss", "p_gain", "pred_thr"})
    if p is None:
        note.append(f"missing predictions for {model}"); continue
    d = hm.merge(p, on=["sample_id", "chrom", "start"], how="inner")
    for s, g in d.groupby("sample_id"):
        rows.append({"tool": model, "sample": s, "depth_x": g.depth_x.iloc[0], "nominal_fraction": g.tumour_fraction_nominal.iloc[0],
                     "replicate": g.replicate.iloc[0], "pearson_r": profile_r(g, value_col(model, g))})
for s, g in hm.groupby("sample_id"):
    rows.append({"tool": "Signal (no model)", "sample": s, "depth_x": g.depth_x.iloc[0], "nominal_fraction": g.tumour_fraction_nominal.iloc[0],
                 "replicate": g.replicate.iloc[0], "pearson_r": profile_r(g, g.lr.values)})
mlr = pd.DataFrame(rows)
cal = pd.read_csv(DD / "sim2_profile_correlation.csv"); prim = pd.read_csv(DD / "sim2_primary_configs.tsv", sep="\t")
cal = cal.merge(prim[["tool", "config", "bin_bp"]], on=["tool", "config", "bin_bp"])
cal["tool"] = cal.tool.map(lambda x: NAME.get(x, x))
cal = cal.rename(columns={"sample": "sample"})[["tool", "sample", "depth_x", "nominal_fraction", "replicate", "pearson_r"]]
allr = pd.concat([cal.assign(kind="caller"), mlr.assign(kind="this study")], ignore_index=True)
allr.to_csv(OUT / "hg008_profile_r_per_sample.tsv", sep="\t", index=False)


def lod_strict(g):
    null = g[g.nominal_fraction == 0].pearson_r
    thr = max(null.max(), null.mean() + 3 * null.std(ddof=1)) if len(null) > 1 else (null.max() if len(null) else np.nan)
    out = []
    for dep, h in g[g.nominal_fraction > 0].groupby("depth_x"):
        mn = h.groupby("nominal_fraction").pearson_r.min(); ok = mn[mn > thr]
        out.append({"depth_x": dep, "threshold": thr, "LoD": float(ok.index.min()) if len(ok) else np.nan})
    return pd.DataFrame(out)


lod = pd.concat([lod_strict(g).assign(tool=t_) for t_, g in allr.groupby("tool")], ignore_index=True)
lod["LoD_or_0.5"] = lod.LoD.fillna(0.5)
lod.to_csv(OUT / "hg008_lod_by_depth.tsv", sep="\t", index=False)
r_tf = allr[allr.nominal_fraction > 0].groupby(["tool", "nominal_fraction"]).pearson_r.median().unstack().round(3)
r_tf.columns = [f"median_r_TF{c:g}" for c in r_tf.columns]
summ = pd.DataFrame({"mean_LoD": lod.groupby("tool")["LoD_or_0.5"].mean().round(4),
                     "median_r_1to10": allr[(allr.nominal_fraction >= 0.01) & (allr.nominal_fraction <= 0.10)].groupby("tool").pearson_r.median().round(4)})
summ = summ.join(r_tf).reset_index().rename(columns={"index": "tool"}).sort_values("median_r_1to10", ascending=False)
summ.to_csv(OUT / "table_hg008_summary.tsv", sep="\t", index=False)
lod.pivot(index="tool", columns="depth_x", values="LoD").to_csv(OUT / "table_hg008_lod_depth_wide.tsv", sep="\t")
allr[allr.nominal_fraction > 0].groupby(["tool", "depth_x", "nominal_fraction"]).pearson_r.median().round(4).reset_index() \
    .to_csv(OUT / "hg008_median_r_tf_depth.tsv", sep="\t", index=False)
chk = summ.set_index("tool").reindex(["QDNAseq", "ichorCNA"])
for tl in chk.index:
    ref = prim.assign(tool=prim.tool.map(lambda x: NAME.get(x, x))).set_index("tool")
    if tl in ref.index:
        note.append(f"check vs manuscript table: {tl} mean LoD {chk.loc[tl, 'mean_LoD']} (paper {ref.loc[tl, 'mean_LoD']:.4f}), "
                    f"median r 1-10% {chk.loc[tl, 'median_r_1to10']} (paper {ref.loc[tl, 'median_r_1to10']:.4f})")

# ================================================================ cancer plasma (GSE71378) vs consensus, 1-Mb windows
W = pd.read_csv(DD / "grch37_1Mb_windows.tsv", sep="\t", dtype={"chrom": str})
W["key"] = W.chrom + ":" + W.start.astype(str)
cons = pd.read_csv(E["GSE_LABELS"], sep="\t", dtype={"chrom": str})
cons = cons[(cons.tool == "consensus") & (cons.level == "full")]
cm = meta("test_cancer")[["sample_id", "base_sample", "depth_x"]].drop_duplicates("sample_id")


def f1(y, p):
    at, ap = y != 0, p != 0
    tp = int(((p == y) & at).sum()); fp = int((ap & (p != y)).sum()); fn = int((at & (p != y)).sum())
    pr = tp / (tp + fp) if tp + fp else np.nan; rc = tp / (tp + fn) if tp + fn else np.nan
    return pr, rc, (2 * pr * rc / (pr + rc) if tp else (0.0 if tp + fn else np.nan))


crow = []
for model, (rdir, kind) in MODELS.items():
    p = preds(rdir, "test_cancer", {"sample_id", "chrom", "start", "pred_hgb", "pred_thr"})
    if p is None:
        continue
    col = "pred_hgb" if kind == "hgb" else "pred_thr"
    p["key"] = p.chrom + ":" + ((p.start // 1_000_000) * 1_000_000).astype(str)
    for s, g in p.groupby("sample_id"):
        mm = cm[cm.sample_id == s]
        if not len(mm):
            continue
        base = mm.base_sample.iloc[0]; level = s[len(base) + 1:] if s.startswith(base + "_") else s
        fr = g.groupby("key")[col].agg(gain=lambda v: (v == 1).mean(), loss=lambda v: (v == -1).mean())
        st = pd.Series(np.where(fr.gain >= 0.5, 1, np.where(fr.loss >= 0.5, -1, 0)), index=fr.index).reindex(W.key).fillna(0).values
        cb = cons[cons["sample"] == base]
        y = pd.Series(cb.state.values, index=(cb.chrom.astype(str) + ":" + cb.start.astype(str))).reindex(W.key).fillna(0).values
        pr_, rc_, f_ = f1(y, st)
        crow.append({"tool": model, "sample": base, "level": level, "depth_x": mm.depth_x.iloc[0], "altered_pct": 100 * float((st != 0).mean()),
                     "consensus_altered_pct": 100 * float((y != 0).mean()), "precision": pr_, "recall": rc_, "F1": f_})
C = pd.DataFrame(crow)
loo = pd.read_csv(DD / "cancer_loo_and_tcga.tsv", sep="\t")
C.to_csv(OUT / "cancer_ml_vs_consensus_per_sample.tsv", sep="\t", index=False)
if len(C):
    full = C[C.level == "full"][["tool", "sample", "F1", "altered_pct"]]
    both = pd.concat([loo[["tool", "sample", "F1", "altered_pct"]].assign(kind="caller (leave-one-out)"),
                      full.assign(kind="this study")], ignore_index=True)
    both.to_csv(OUT / "table_cancer_f1_native_depth.tsv", sep="\t", index=False)
    canc = both[~both["sample"].str.startswith("IH")]
    canc.groupby("tool").agg(median_F1=("F1", "median"), mean_F1=("F1", "mean"), samples=("sample", "nunique")).round(3) \
        .sort_values("median_F1", ascending=False).reset_index().to_csv(OUT / "table_cancer_summary.tsv", sep="\t", index=False)
    C[~C["sample"].str.startswith("IH")].groupby(["tool", "level"]).F1.median().round(3).unstack() \
        .to_csv(OUT / "table_cancer_ml_f1_by_depth.tsv", sep="\t")

# ================================================================ non-cancer plasma (PRJNA795275): % genome called altered
hb = pd.read_csv(DD / "healthy_burden.tsv", sep="\t")
keep = pd.Series(False, index=hb.index)
for tl, cf in HEALTHY_PRIMARY:
    keep |= (hb.tool == tl) & (hb.config == cf)
hb = hb[keep][["tool", "base", "level", "FGA_pct"]].rename(columns={"base": "sample"})
nm_ = meta("test_noncancer")[["sample_id", "base_sample"]].drop_duplicates("sample_id")
nrow = []
for model, (rdir, kind) in MODELS.items():
    p = preds(rdir, "test_noncancer", {"sample_id", "pred_hgb", "pred_thr"})
    if p is None:
        continue
    col = "pred_hgb" if kind == "hgb" else "pred_thr"
    for s, g in p.groupby("sample_id"):
        b = nm_[nm_.sample_id == s].base_sample.iloc[0]
        nrow.append({"tool": model, "sample": b, "level": s[len(b) + 1:] if s.startswith(b + "_") else s,
                     "FGA_pct": 100 * float((g[col] != 0).mean())})
N = pd.concat([hb.assign(kind="caller"), pd.DataFrame(nrow).assign(kind="this study")], ignore_index=True)
N.to_csv(OUT / "noncancer_pct_altered_per_sample.tsv", sep="\t", index=False)
N.groupby(["tool", "level"]).FGA_pct.median().round(3).unstack().to_csv(OUT / "table_noncancer_median_pct_altered.tsv", sep="\t")

# ================================================================ all splits, both models (bin level, 100 kb)
pm = []
for v, rdir in (("v1", R1), ("v2", R2)):
    f = rdir / "per_sample_metrics.tsv"
    if f.exists():
        x = pd.read_csv(f, sep="\t"); x["model"] = x.model.map(lambda m: f"ML {v} " + ("gradient boosting" if m.startswith("HGB") else "threshold")); pm.append(x)
PM = pd.concat(pm, ignore_index=True) if pm else pd.DataFrame()
if len(PM):
    PM.groupby(["split", "model"]).agg(samples=("sample_id", "nunique"), F1=("f1", "median"), precision=("precision", "median"),
                                       recall=("recall", "median"), pct_altered_pred=("pga_pred", lambda x: 100 * x.median()),
                                       pct_altered_true=("pga_true", lambda x: 100 * x.median())).round(3).reset_index() \
        .to_csv(OUT / "table_all_splits_bin_level.tsv", sep="\t", index=False)
    s1 = PM[PM.split == "test_sim1_ih03"]
    if len(s1):
        s1.groupby(["model", "depth_x"]).f1.median().round(3).unstack().to_csv(OUT / "table_sim1_ih03_f1_by_depth.tsv", sep="\t")

# ================================================================ key numbers (facts for writing the Results; not manuscript text)
def g_(tab, tool, col):
    try:
        return tab.set_index("tool").loc[tool, col]
    except Exception:
        return "NA"

info1 = json.loads((R1 / "run_info.json").read_text()) if (R1 / "run_info.json").exists() else {}
info2 = json.loads((R2 / "run_info.json").read_text()) if (R2 / "run_info.json").exists() else {}
k = ["# Key numbers (auto-generated by f_report.py; facts for writing, not manuscript text)", ""]
k.append(f"Dataset: {len(pd.read_csv(MLO / 'manifest.tsv', sep=chr(9)))} BAMs; 100-kb and 1-Mb configs; splits fixed before training.")
k.append(f"v1 trained on {info1.get('train_samples', 'NA')} Simulation 1 samples ({info1.get('train_bins', 'NA')} bins); "
         f"v2 adds {info2.get('n_train_aug', 'NA')} augmented samples from {len(info2.get('augmentation_bases', []))} negative controls.")
k += ["", "## HG008 (manuscript definitions)", "", "| tool | mean LoD | median r (1-10 % TF) |", "|---|---|---|"]
for r_ in summ.itertuples():
    k.append(f"| {r_.tool} | {r_.mean_LoD} | {r_.median_r_1to10} |")
if len(C):
    cs = pd.read_csv(OUT / "table_cancer_summary.tsv", sep="\t")
    k += ["", "## Cancer plasma, F1 vs consensus (native depth; callers leave-one-out)", "", "| tool | median F1 |", "|---|---|"]
    k += [f"| {r_.tool} | {r_.median_F1} |" for r_ in cs.itertuples()]
nt = N[N.level == "full"].groupby("tool").FGA_pct.median().round(3)
k += ["", "## Non-cancer plasma, median % genome called altered (native depth)", "", "| tool | median % |", "|---|---|"]
k += [f"| {i} | {v} |" for i, v in nt.items()]
if len(PM):
    tb = pd.read_csv(OUT / "table_all_splits_bin_level.tsv", sep="\t")
    k += ["", "## Bin-level F1 (100 kb, median over samples)", "", "| split | model | F1 |", "|---|---|---|"]
    k += [f"| {r_.split} | {r_.model} | {r_.F1} |" for r_ in tb.itertuples()]
k += ["", "## Checks"] + [f"- {n}" for n in note]
(OUT / "key_numbers.md").write_text("\n".join(k) + "\n")

(OUT / "methods_notes.md").write_text(f"""# Methods notes (facts for writing the Methods; rewrite in your own words)

Dataset. 310 BAMs from the four arms. Counting unit: one fragment (single-end read or read 1 of a pair); unmapped,
secondary, supplementary, QC-fail and duplicate reads skipped; MAPQ >= 20 for the depth signal. 10-kb bins on autosomes,
masked when N > 5 %, ENCODE hg38 blacklist overlap > 10 % (GRCh38 only), or mappability proxy (MAPQ>=20 share of all
fragments, pooled over all samples of a build, no labels) < 0.75. Usable 10-kb bins summed into 100-kb and 1-Mb bins,
bins with < 50 % usable length dropped. GC correction: count divided by the median count of its GC-quantile (20 groups),
log2, sample median set to 0. Features: log2 ratio, rolling medians within chromosome (3/9/25 bins in v1; 3/9/25/51/101
in v2), rolling mean (9), chromosome median, sample noise (MAD of successive differences / sqrt 2), noise-scaled medians,
log depth, GC, mappability proxy, usable fraction; v1 also fragment-length features (short <= 150 bp fraction, mean length).
Labels: Simulation 1 spike-in truth (event covering >= 50 % of a bin); HG008 NIST draft somatic benchmark, gain/loss
relative to the modal copy number ({PLOIDY:.2f} mean ploidy), bins < 50 % covered unlabelled, 0 % mixtures neutral;
PRJNA795275 all neutral; GSE71378 consensus of the five callers (weak labels).
Splits: train = Simulation 1 IH01+IH02 replicates 1-4 + negative controls; validation = replicate 5; tests = IH03, HG008,
PRJNA795275, GSE71378. No test split used for any choice.
Models. scikit-learn {info1.get('sklearn', '')} HistGradientBoostingClassifier (3 classes), seed {info1.get('seed', '')} (v1) /
{info2.get('seed', '')} (v2); single altered bins between unaltered neighbours set to neutral. Threshold model: |noise-scaled
rolling median| above a cut-off ({info1.get('threshold_lr_z9', 'NA')} for v1 on 9 bins, {info2.get('threshold_lr_z25', 'NA')} for v2
on 25 bins) chosen on validation. v1: class-balanced weights, trained on Simulation 1 only. v2: {info2.get('n_train_aug', 'NA')}
augmented training + {info2.get('n_val_aug', 'NA')} augmented validation samples from the Simulation 1 negative controls:
binomial thinning to 0.1x-native depth; random layouts (0-60 % genome altered, 0.3 Mb to whole chromosome, CN 0/1/3/4/5)
diluted by a tumour fraction of 3-60 % (15 % of samples 100 %): expected ratio ((1-TF)*2 + TF*CN)/2, losses by binomial
thinning, gains by added Poisson counts; square-root class weights; decision P(altered) > {info2.get('tau', 'NA')} chosen on
augmented validation.
Scoring. HG008 profile correlation, strict LoD and median r (1-10 %) exactly as for the callers (1-Mb windows >= 80 %
benchmark coverage; model value P(gain)-P(loss)); cancer F1 vs consensus on 1-Mb windows (window altered when >= 50 %
of its 100-kb bins agree); non-cancer % genome called altered.
""")
print((OUT / "key_numbers.md").read_text())
