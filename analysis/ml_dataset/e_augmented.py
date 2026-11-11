#!/usr/bin/env python
"""E  Augmented baseline (model v2): the same feature pipeline, trained on tumour-fraction- and depth-augmented data.

Why: model v1 (step D) was trained only on full-strength Simulation 1 CNVs and learned to ignore the weak shifts
of low-tumour-fraction samples, although the smoothed signal itself tracks the HG008 truth well (r up to 0.95).

Augmentation (only from TRAINING data: Simulation 1 negative controls of IH01 + IH02, 100-kb counts):
  1. depth:   counts are binomially thinned to a random lower depth (0.1x .. the control's own depth)
  2. CNVs:    a random genome layout is drawn (0-60 % of the genome altered, events 0.3 Mb .. whole chromosome,
              copy number 0/1 or 3/4/5), diluted by a random tumour fraction (3-60 %, 15 % of samples at 100 %)
              expected ratio r = ((1-TF)*2 + TF*CN) / 2 ; losses: binomial thinning, gains: + Poisson extra reads
  3. labels:  gain / loss of the TUMOUR genome (independent of how weak the signal is)
Features are recomputed with the same formulas as step C (GC correction, rolling medians) for every split,
plus wider windows (51, 101 bins) and noise-scaled versions, because tumour CNAs are often arm-sized.
Fragment-length features are not used (the simulations carry no fragment signal).
The decision rule (call altered when max(P(gain), P(loss)) > tau) and the threshold model are tuned on an
augmented VALIDATION set (new random layouts). No test split is used for any choice.
Output: $MLO/results_v2/  and  $MLO/ml_results_v2_for_claude.tar.gz
"""
import os, sys, json, time, tarfile, platform
from pathlib import Path
import numpy as np, pandas as pd
import sklearn
try:
    from sklearn.ensemble import HistGradientBoostingClassifier
except ImportError:
    from sklearn.experimental import enable_hist_gradient_boosting  # noqa: F401
    from sklearn.ensemble import HistGradientBoostingClassifier

E = os.environ
MLO = Path(E["MLO"]); D = MLO / "hf_dataset" / "data" / "bins_100kb"
R = MLO / "results_v2"; (R / "predictions").mkdir(parents=True, exist_ok=True)
SEED = 20261111
N_TRAIN_AUG = int(E.get("N_TRAIN_AUG", 300)); N_VAL_AUG = int(E.get("N_VAL_AUG", 60))
B = 100000
TESTS = ["validation", "test_sim1_ih03", "test_hg008", "test_noncancer", "test_cancer"]
WIN = [3, 9, 25, 51, 101]
FEATS = (["lr", "lr_mean9", "lr_chrom", "noise", "log_depth_bin", "gc", "map_proxy", "frac_usable"]
         + [f"lr_med{w}" for w in WIN] + [f"lr_z{w}" for w in WIN])
KEEP_META = ["dataset", "split", "sample_id", "base_sample", "class", "replicate", "depth_x",
             "tumour_fraction_nominal", "tumour_fraction_realised", "build"]


# ---------------------------------------------------------------- features (same formulas as c_assemble.py)
def features(a):
    """a: one sample, 100-kb bins sorted by chrom/start, with n_q20, frac_usable, gc"""
    a = a.reset_index(drop=True)
    c = a.n_q20 / a.frac_usable
    q = pd.qcut(a.gc.rank(method="first"), 20, labels=False, duplicates="drop") if len(a) >= 40 else pd.Series(0, index=a.index)
    exp = c.groupby(q).transform("median")
    lr = np.log2((c + 0.5) / (exp + 0.5)); lr = (lr - np.nanmedian(lr)).values
    ch = a.chrom.values
    s = pd.Series(lr)
    f = {"lr": lr}
    for w in WIN:
        f[f"lr_med{w}"] = s.groupby(ch).transform(lambda v: v.rolling(w, center=True, min_periods=1).median()).values
    f["lr_mean9"] = s.groupby(ch).transform(lambda v: v.rolling(9, center=True, min_periods=1).mean()).values
    f["lr_chrom"] = s.groupby(ch).transform("median").values
    dl = np.diff(lr)
    noise = float(np.nanmedian(np.abs(dl - np.nanmedian(dl))) * 1.4826 / np.sqrt(2)) if len(dl) > 10 else np.nan
    f["noise"] = np.full(len(a), noise)
    for w in WIN:                                   # median of w bins has sd ~ 1.2533 * noise / sqrt(w)
        f[f"lr_z{w}"] = f[f"lr_med{w}"] / (1.2533 * noise / np.sqrt(w)) if noise and noise > 0 else np.nan
    f["log_depth_bin"] = np.full(len(a), np.log10(np.nanmedian(c) + 1))
    out = pd.DataFrame(f)
    out["gc"] = a.gc.values; out["map_proxy"] = a.map_proxy.values; out["frac_usable"] = a.frac_usable.values
    return out


def load(split):
    f = D / f"{split}.csv.gz"
    return pd.read_csv(f, dtype={"chrom": str, "sample_id": str, "base_sample": str, "replicate": str, "event_id": str}) if f.exists() else None


def refeature(d):
    parts = []
    for s, g in d.groupby("sample_id", sort=False):
        ff = features(g)
        g = g.reset_index(drop=True)
        base = g[KEEP_META + ["chrom", "start", "label_state", "truth_log2", "event_id", "event_size"]]
        parts.append(pd.concat([base, ff], axis=1))
    return pd.concat(parts, ignore_index=True)


# ---------------------------------------------------------------- augmentation
def layout(bins, rng):
    """random CNV layout -> per-bin copy number (baseline 2)"""
    cn = np.full(len(bins), 2.0)
    target = rng.uniform(0, 0.6)
    if rng.random() < 0.1:
        return cn
    chrom = bins.chrom.values; idx = {c: np.where(chrom == c)[0] for c in pd.unique(chrom)}
    tries = 0
    while (cn != 2).mean() < target and tries < 200:
        tries += 1
        c = rng.choice(list(idx)); ii = idx[c]; n = len(ii)
        u = rng.random()
        if u < 0.15:
            L = n                                            # whole chromosome
        elif u < 0.35:
            L = max(3, n // 2)                               # ~arm
        else:
            L = int(np.clip(np.exp(rng.uniform(np.log(3), np.log(1500))), 3, n))   # 0.3 Mb .. 150 Mb
        st = rng.integers(0, max(1, n - L + 1))
        v = rng.choice([0, 1, 3, 4, 5], p=[0.04, 0.44, 0.34, 0.13, 0.05])
        cn[ii[st:st + L]] = v
    return cn


def augment(neg, rng, k):
    g = neg.reset_index(drop=True)
    dep0 = float(g.depth_x.iloc[0])
    lo = min(0.1, dep0)
    dep = float(np.exp(rng.uniform(np.log(lo), np.log(dep0)))) if dep0 > lo else dep0
    c = rng.binomial(g.n_q20.values.astype(np.int64), min(1.0, dep / dep0))
    cn = layout(g, rng)
    tf = 1.0 if rng.random() < 0.15 else float(np.exp(rng.uniform(np.log(0.03), np.log(0.6))))
    r = ((1 - tf) * 2 + tf * cn) / 2
    out = c.copy()
    lo_m = r < 1; hi_m = r > 1
    out[lo_m] = rng.binomial(c[lo_m], r[lo_m])
    out[hi_m] = c[hi_m] + rng.poisson(c[hi_m] * (r[hi_m] - 1))
    a = g[["chrom", "start", "gc", "map_proxy", "frac_usable"]].copy(); a["n_q20"] = out
    ff = features(a)
    meta = pd.DataFrame({"dataset": "sim1_augmented", "sample_id": f"AUG{k:04d}_{g.sample_id.iloc[0]}",
                         "base_sample": g.base_sample.iloc[0], "depth_x": round(dep, 4), "tumour_fraction_nominal": round(tf, 4),
                         "chrom": a.chrom.values, "start": a.start.values, "label_state": np.sign(cn - 2),
                         "label_cn": cn}, index=ff.index)
    return pd.concat([meta, ff], axis=1)


def scores(y, p):
    m = ~np.isnan(y); y = y[m]; p = p[m]
    if len(y) == 0:
        return dict(n_bins=0)
    at, ap = y != 0, p != 0
    tp = int(((p == y) & at).sum()); fp = int((ap & (p != y)).sum()); fn = int((at & (p != y)).sum())
    pr = tp / (tp + fp) if tp + fp else np.nan; rc = tp / (tp + fn) if tp + fn else np.nan
    f1 = 2 * pr * rc / (pr + rc) if (tp and pr + rc) else (0.0 if (tp + fn) else np.nan)
    return dict(n_bins=int(len(y)), true_altered=int(at.sum()), precision=pr, recall=rc, f1=f1,
                pga_pred=float(ap.mean()), pga_true=float(at.mean()), fp_bin_rate=float((ap & ~at).sum() / max((~at).sum(), 1)))


def postproc(d, pred):
    p = pd.Series(pred, index=d.index); key = (d.sample_id + ":" + d.chrom).values
    prev = p.groupby(key).shift(1); nxt = p.groupby(key).shift(-1)
    p[(p != 0) & (prev != p) & (nxt != p)] = 0
    return p.values


def pooled_f1(d, pred):
    return scores(d.label_state.values.astype(float), pred.astype(float)).get("f1", np.nan)


# ---------------------------------------------------------------- build augmented train / validation
t0 = time.time()
tr = load("train")
if tr is None:
    sys.exit("train split missing - run step C first")
negs = [g for _, g in tr[tr["class"] == "negative"].groupby("sample_id")]
if not negs:
    sys.exit("no Simulation 1 negative controls in the train split - cannot augment")
print(f"{len(negs)} negative controls as augmentation bases: " + ", ".join(g.sample_id.iloc[0] for g in negs), flush=True)
rng = np.random.default_rng(SEED)
aug_tr = pd.concat([augment(negs[k % len(negs)], rng, k) for k in range(N_TRAIN_AUG)], ignore_index=True)
rng_v = np.random.default_rng(SEED + 1)
aug_va = pd.concat([augment(negs[k % len(negs)], rng_v, 100000 + k) for k in range(N_VAL_AUG)], ignore_index=True)
orig = refeature(tr[tr.label_state.notna()])
train = pd.concat([orig[["sample_id", "chrom", "label_state"] + FEATS], aug_tr[["sample_id", "chrom", "label_state"] + FEATS]], ignore_index=True)
train = train[train.label_state.notna()].reset_index(drop=True)
(MLO / "augmented").mkdir(exist_ok=True)
aug_tr.round(5).to_csv(MLO / "augmented" / "train_augmented_100kb.csv.gz", index=False, compression="gzip")
aug_va.round(5).to_csv(MLO / "augmented" / "validation_augmented_100kb.csv.gz", index=False, compression="gzip")
print(f"augmented: {N_TRAIN_AUG} train + {N_VAL_AUG} validation samples; training bins {len(train)} ({time.time() - t0:.0f}s)", flush=True)

# ---------------------------------------------------------------- train
y = (train.label_state.values + 1).astype(int)
cnt = np.bincount(y, minlength=3).clip(min=1)
cw = np.sqrt(len(y) / (3 * cnt))                        # softened class balance (keeps false positives down)
hgb = HistGradientBoostingClassifier(max_iter=500, learning_rate=0.08, max_leaf_nodes=63, min_samples_leaf=500,
                                     l2_regularization=1.0, early_stopping=True, validation_fraction=0.1,
                                     n_iter_no_change=25, random_state=SEED)
hgb.fit(train[FEATS].values, y, sample_weight=cw[y])
print(f"HGB v2 trained, {hgb.n_iter_} iterations, class counts {cnt.tolist()} ({time.time() - t0:.0f}s)", flush=True)
del train

# ---------------------------------------------------------------- tune decision rules on augmented validation
pv = hgb.predict_proba(aug_va[FEATS].values)


def decide(pr, tau):
    alt = np.maximum(pr[:, 0], pr[:, 2]) > tau
    return np.where(alt, np.where(pr[:, 2] >= pr[:, 0], 1, -1), 0)


best_tau = max(((pooled_f1(aug_va, postproc(aug_va, decide(pv, t))), t) for t in np.arange(0.2, 0.96, 0.05)), key=lambda x: (np.nan_to_num(x[0], nan=-1), -x[1]))
TAU = round(float(best_tau[1]), 2)
best_thr = max(((pooled_f1(aug_va, postproc(aug_va, np.where(aug_va.lr_z25 > t, 1, np.where(aug_va.lr_z25 < -t, -1, 0)))), t)
                for t in np.arange(1.0, 15.01, 0.5)), key=lambda x: (np.nan_to_num(x[0], nan=-1), x[1]))
THR = float(best_thr[1])
print(f"tuned on augmented validation: HGB tau {TAU} (F1 {best_tau[0]:.3f}); threshold |lr_z25| > {THR} (F1 {best_thr[0]:.3f})", flush=True)

# ---------------------------------------------------------------- evaluate on every split
per = []
for split in TESTS + ["validation_augmented"]:
    if split == "validation_augmented":
        d = aug_va.assign(split="validation_augmented", **{k: np.nan for k in ("truth_log2",)})
        for k in KEEP_META:
            if k not in d:
                d[k] = np.nan
        d["split"] = "validation_augmented"
    else:
        raw = load(split)
        if raw is None:
            print("no split", split); continue
        d = refeature(raw)
    pr = hgb.predict_proba(d[FEATS].values)
    d = d.assign(p_loss=pr[:, 0], p_gain=pr[:, 2], pred_hgb=postproc(d, decide(pr, TAU)),
                 pred_thr=postproc(d, np.where(d.lr_z25 > THR, 1, np.where(d.lr_z25 < -THR, -1, 0))))
    d[["sample_id", "chrom", "start", "label_state", "p_loss", "p_gain", "pred_hgb", "pred_thr"]].round(4) \
        .to_csv(R / "predictions" / f"{split}.csv.gz", index=False, compression="gzip")
    for s, g in d.groupby("sample_id", sort=False):
        meta = {k: g[k].iloc[0] for k in KEEP_META}
        yv = g.label_state.values.astype(float)
        for model, col in (("HGB_v2", "pred_hgb"), ("THR_v2", "pred_thr")):
            row = {**meta, "model": model, **scores(yv, g[col].values.astype(float))}
            if split == "test_hg008" and meta["tumour_fraction_nominal"] != 0:
                ok = (g.truth_log2.notna() & np.isfinite(g.truth_log2)).values
                if ok.sum() > 50:
                    sc = (g.p_gain - g.p_loss) if model == "HGB_v2" else g.lr_med25
                    row["pearson_r_profile"] = float(np.corrcoef(sc.values[ok], g.truth_log2.values[ok])[0, 1])
            per.append(row)
    print("scored", split, d.sample_id.nunique(), "samples", flush=True)

P = pd.DataFrame(per); P.to_csv(R / "per_sample_metrics.tsv", sep="\t", index=False)
agg = dict(samples=("sample_id", "nunique"), f1_median=("f1", "median"), precision_median=("precision", "median"),
           recall_median=("recall", "median"), pga_pred_median=("pga_pred", "median"), pga_true_median=("pga_true", "median"),
           fp_bin_rate_median=("fp_bin_rate", "median"))
P.groupby(["split", "model"]).agg(**agg).reset_index().round(4).to_csv(R / "summary_by_split.tsv", sep="\t", index=False)
P.groupby(["split", "model", "depth_x"]).agg(**agg).reset_index().round(4).to_csv(R / "summary_by_split_depth.tsv", sep="\t", index=False)
h = P[P.split == "test_hg008"]
if len(h):
    h.groupby(["model", "tumour_fraction_nominal", "depth_x"]).agg(**dict(agg, r_profile_median=("pearson_r_profile", "median"))) \
        .reset_index().round(4).to_csv(R / "summary_hg008_tf_depth.tsv", sep="\t", index=False)
c = P[P.split == "test_cancer"]
if len(c):
    c.groupby(["model", "base_sample", "depth_x"]).agg(**agg).reset_index().round(4).to_csv(R / "summary_cancer_sample_depth.tsv", sep="\t", index=False)
try:
    import joblib; joblib.dump(hgb, R / "model_v2_hgb.joblib")
except Exception as ex:
    print("model not saved:", ex)
info = {"date": time.strftime("%Y-%m-%d %H:%M"), "python": platform.python_version(), "sklearn": sklearn.__version__,
        "pandas": pd.__version__, "numpy": np.__version__, "seed": SEED, "features": FEATS, "n_train_aug": N_TRAIN_AUG,
        "n_val_aug": N_VAL_AUG, "augmentation_bases": [g.sample_id.iloc[0] for g in negs], "tau": TAU,
        "tau_val_f1": best_tau[0], "threshold_lr_z25": THR, "threshold_val_f1": best_thr[0], "hgb_iterations": int(hgb.n_iter_),
        "class_weights": cw.round(4).tolist(), "minutes": round((time.time() - t0) / 60, 1)}
(R / "run_info.json").write_text(json.dumps(info, indent=1) + "\n")
tb = MLO / "ml_results_v2_for_claude.tar.gz"
with tarfile.open(tb, "w:gz") as t:
    for f in sorted(R.glob("*.tsv")) + [R / "run_info.json"]:
        t.add(f, arcname=f"ml_results_v2/{f.name}")
print("\n=== summary_by_split (v2) ===")
print(pd.read_csv(R / "summary_by_split.tsv", sep="\t").to_string(index=False))
print("\nbundle to send back:", tb)
