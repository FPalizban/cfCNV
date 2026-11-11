#!/usr/bin/env python
"""D  Baseline models on the 100-kb dataset: trained on `train`, tuned on `validation`, scored on every test split.

Two baselines, same post-processing (runs of a single altered bin are set back to neutral):
  THR  no learning: |lr_z9| above a threshold -> gain/loss by sign; the threshold is chosen on validation
  HGB  scikit-learn HistGradientBoostingClassifier (3 classes, class-balanced weights, fixed seed)
Bin-level scores per sample (altered bins only count as correct when the direction matches):
  precision, recall, F1, predicted fraction of genome altered (PGA), false-positive bin rate
Simulation 1 also gets event-level recall (event found when >= 50 % of its bins are called in its direction).
HG008 also gets the Pearson r between the model score P(gain) - P(loss) and log2(CN / ploidy) of the benchmark,
and the same r for the plain smoothed log2 ratio (lr_med9) as a no-model reference.
Output: $MLO/results/   and   $MLO/ml_results_for_claude.tar.gz  (small tables only, to send back)
"""
import os, sys, json, time, tarfile, platform
from pathlib import Path
import numpy as np, pandas as pd
import sklearn
try:
    from sklearn.ensemble import HistGradientBoostingClassifier
except ImportError:                                   # scikit-learn < 1.0
    from sklearn.experimental import enable_hist_gradient_boosting  # noqa: F401
    from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.inspection import permutation_importance

E = os.environ
MLO = Path(E["MLO"]); D = MLO / "hf_dataset" / "data" / "bins_100kb"
R = MLO / "results"; (R / "predictions").mkdir(parents=True, exist_ok=True)
FEATS = ["lr", "lr_med3", "lr_med9", "lr_med25", "lr_mean9", "lr_chrom", "noise", "lr_z9", "log_depth_bin",
         "gc", "map_proxy", "frac_usable", "short_frac_lr", "short_frac_med9", "flen_dev"]
TESTS = ["validation", "test_sim1_ih03", "test_hg008", "test_noncancer", "test_cancer"]
SEED = 20261110


def load(split):
    f = D / f"{split}.csv.gz"
    if not f.exists():
        return None
    return pd.read_csv(f, dtype={"chrom": str, "sample_id": str, "base_sample": str, "replicate": str, "event_id": str})


def postproc(d, pred):
    """single altered bins (both neighbours differ) -> neutral, within a sample and chromosome"""
    p = pd.Series(pred, index=d.index)
    key = d.sample_id + ":" + d.chrom
    prev = p.groupby(key.values).shift(1); nxt = p.groupby(key.values).shift(-1)
    lone = (p != 0) & (prev != p) & (nxt != p)
    p[lone] = 0
    return p.values


def scores(y, p):
    m = ~np.isnan(y); y = y[m]; p = p[m]
    if len(y) == 0:
        return dict(n_bins=0)
    at, ap = y != 0, p != 0
    tp = int(((p == y) & at).sum()); fp = int((ap & (p != y)).sum()); fn = int((at & (p != y)).sum())
    pr = tp / (tp + fp) if tp + fp else np.nan; rc = tp / (tp + fn) if tp + fn else np.nan
    f1 = 2 * pr * rc / (pr + rc) if (tp and pr + rc) else (0.0 if (tp + fn) else np.nan)
    return dict(n_bins=int(len(y)), true_altered=int(at.sum()), precision=pr, recall=rc, f1=f1,
                pga_pred=float(ap.mean()), pga_true=float(at.mean()),
                fp_bin_rate=float((ap & ~at).sum() / max((~at).sum(), 1)))


# ---------------------------------------------------------------- train
tr = load("train"); va = load("validation")
if tr is None or va is None:
    sys.exit("train/validation split missing - run step C first")
tr = tr[tr.label_state.notna()].reset_index(drop=True)
y = (tr.label_state.values + 1).astype(int)
cw = len(y) / (3 * np.bincount(y, minlength=3).clip(min=1))
print(f"train: {tr.sample_id.nunique()} samples, {len(tr)} bins, class counts {np.bincount(y, minlength=3).tolist()}", flush=True)
t0 = time.time()
hgb = HistGradientBoostingClassifier(max_iter=400, learning_rate=0.08, max_leaf_nodes=31, min_samples_leaf=200,
                                     l2_regularization=1.0, early_stopping=True, validation_fraction=0.1,
                                     n_iter_no_change=20, random_state=SEED)
hgb.fit(tr[FEATS].values, y, sample_weight=cw[y])
print(f"HGB trained in {time.time() - t0:.0f}s, {hgb.n_iter_} iterations", flush=True)

# threshold baseline: choose t on validation
va = va[va.label_state.notna()].reset_index(drop=True)
best = (-1, None)
for t in np.round(np.arange(0.5, 12.01, 0.25), 2):
    z = va.lr_z9.fillna(0).values
    p = postproc(va, np.where(z > t, 1, np.where(z < -t, -1, 0)))
    f = scores(va.label_state.values, p).get("f1", np.nan)
    if f == f and f > best[0]:
        best = (f, t)
THR = best[1] if best[1] is not None else 3.0
print(f"threshold baseline: |lr_z9| > {THR} (validation F1 {best[0]:.3f})", flush=True)


def predict(d):
    pr = hgb.predict_proba(d[FEATS].values)
    pred_h = postproc(d, pr.argmax(1) - 1)
    z = d.lr_z9.fillna(0).values
    pred_t = postproc(d, np.where(z > THR, 1, np.where(z < -THR, -1, 0)))
    return pr, pred_h, pred_t


# ---------------------------------------------------------------- evaluate
per, events, imp_done = [], [], False
SKEYS = ["dataset", "split", "sample_id", "base_sample", "class", "replicate", "depth_x",
         "tumour_fraction_nominal", "tumour_fraction_realised", "build"]
for split in TESTS + ["train"]:
    d = va if split == "validation" else load(split)
    if d is None:
        print("no split", split); continue
    if split == "train":
        d = tr
    pr, ph, pt = predict(d)
    d = d.assign(p_loss=pr[:, 0], p_neutral=pr[:, 1], p_gain=pr[:, 2], pred_hgb=ph, pred_thr=pt)
    d[["sample_id", "chrom", "start", "label_state", "truth_log2", "lr_med9", "p_loss", "p_neutral", "p_gain", "pred_hgb", "pred_thr"]] \
        .round(4).to_csv(R / "predictions" / f"{split}.csv.gz", index=False, compression="gzip")
    for s, g in d.groupby("sample_id", sort=False):
        meta = {k: g[k].iloc[0] for k in SKEYS}
        yv = g.label_state.values.astype(float)
        for model, col in (("HGB", "pred_hgb"), ("THR", "pred_thr")):
            row = {**meta, "model": model, **scores(yv, g[col].values.astype(float))}
            if split == "test_hg008" and g.truth_log2.notna().sum() > 50 and meta["tumour_fraction_nominal"] != 0:
                ok = g.truth_log2.notna() & np.isfinite(g.truth_log2)
                sc = (g.p_gain - g.p_loss) if model == "HGB" else g.lr_med9
                row["pearson_r_profile"] = float(np.corrcoef(sc[ok], g.truth_log2[ok])[0, 1])
            per.append(row)
        if split in ("validation", "test_sim1_ih03", "train") and g.event_id.notna().any():
            for e, ge in g[g.event_id.notna()].groupby("event_id"):
                want = ge.label_state.iloc[0]
                for model, col in (("HGB", "pred_hgb"), ("THR", "pred_thr")):
                    events.append({**meta, "model": model, "event_id": e, "event_size": ge.event_size.iloc[0],
                                   "direction": "gain" if want > 0 else "loss", "bins": len(ge),
                                   "detected": bool((ge[col] == want).mean() >= 0.5)})
    if split == "validation" and not imp_done:
        sub = d.sample(min(len(d), 200000), random_state=SEED)
        pi = permutation_importance(hgb, sub[FEATS].values, (sub.label_state.values + 1).astype(int), n_repeats=3,
                                    random_state=SEED, scoring="f1_macro", n_jobs=1)
        pd.DataFrame({"feature": FEATS, "importance_mean": pi.importances_mean, "importance_sd": pi.importances_std}) \
            .sort_values("importance_mean", ascending=False).to_csv(R / "feature_importance_validation.tsv", sep="\t", index=False)
        imp_done = True
    print("scored", split, d.sample_id.nunique(), "samples", flush=True)

P = pd.DataFrame(per); P.to_csv(R / "per_sample_metrics.tsv", sep="\t", index=False)
agg = dict(samples=("sample_id", "nunique"), f1_median=("f1", "median"), precision_median=("precision", "median"),
           recall_median=("recall", "median"), pga_pred_median=("pga_pred", "median"), pga_true_median=("pga_true", "median"),
           fp_bin_rate_median=("fp_bin_rate", "median"))
P.groupby(["split", "model"]).agg(**agg).reset_index().round(4).to_csv(R / "summary_by_split.tsv", sep="\t", index=False)
P.groupby(["split", "model", "depth_x"]).agg(**agg).reset_index().round(4).to_csv(R / "summary_by_split_depth.tsv", sep="\t", index=False)
h = P[P.split == "test_hg008"]
if len(h):
    a2 = dict(agg, r_profile_median=("pearson_r_profile", "median"))
    h.groupby(["model", "tumour_fraction_nominal", "depth_x"]).agg(**a2).reset_index().round(4) \
        .to_csv(R / "summary_hg008_tf_depth.tsv", sep="\t", index=False)
c = P[P.split == "test_cancer"]
if len(c):
    c.groupby(["model", "base_sample", "depth_x"]).agg(**agg).reset_index().round(4).to_csv(R / "summary_cancer_sample_depth.tsv", sep="\t", index=False)
if events:
    ev = pd.DataFrame(events); ev.to_csv(R / "sim1_event_detection.tsv", sep="\t", index=False)
    ev["size_kb"] = pd.cut(ev.event_size / 1000, [0, 300, 1000, 5000, 1e9], labels=["<=300kb", "0.3-1Mb", "1-5Mb", ">5Mb"])
    ev.groupby(["split", "model", "depth_x", "size_kb"], observed=True).detected.mean().reset_index() \
        .round(4).to_csv(R / "summary_sim1_event_recall.tsv", sep="\t", index=False)

try:
    import joblib; joblib.dump(hgb, R / "baseline_hgb.joblib")
except Exception as ex:
    print("model not saved:", ex)
info = {"date": time.strftime("%Y-%m-%d %H:%M"), "python": platform.python_version(), "sklearn": sklearn.__version__,
        "pandas": pd.__version__, "numpy": np.__version__, "seed": SEED, "features": FEATS, "threshold_lr_z9": THR,
        "threshold_validation_f1": best[0], "hgb_iterations": int(hgb.n_iter_), "class_weights": cw.round(4).tolist(),
        "train_samples": int(tr.sample_id.nunique()), "train_bins": int(len(tr))}
(R / "run_info.json").write_text(json.dumps(info, indent=1) + "\n")

# ---------------------------------------------------------------- small bundle to send back
tb = MLO / "ml_results_for_claude.tar.gz"
with tarfile.open(tb, "w:gz") as t:
    for f in sorted(R.glob("*.tsv")) + [R / "run_info.json", MLO / "hf_dataset" / "README.md", MLO / "hf_dataset" / "samples.tsv",
                                        MLO / "manifest.tsv"]:
        if f.exists():
            t.add(f, arcname=f"ml_results/{f.name}")
print("\n=== summary_by_split ===")
print(pd.read_csv(R / "summary_by_split.tsv", sep="\t").to_string(index=False))
print("\nbundle to send back:", tb)
