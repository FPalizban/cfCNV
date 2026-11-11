#!/usr/bin/env python
"""C  Builds the Hugging Face-style dataset from the 10-kb counts: features + labels + splits.

Two resolutions (configs): bins_100kb (used by the baseline model) and bins_1Mb.
Per sample:
  1. 10-kb bins are masked when N > 5 %, ENCODE blacklist > 10 % (GRCh38), or the mappability proxy of
     that build (sum MAPQ>=20 / sum all fragments over ALL samples of that build; no labels used) < 0.75.
  2. Usable 10-kb bins are summed into 100-kb / 1-Mb bins; bins with < 50 % usable length are dropped.
  3. GC correction: counts / median count of the bin's GC-quantile (20 groups), log2, sample median = 0.
  4. Neighbourhood features (rolling medians within the chromosome), noise, fragment-length features.
Labels (state: -1 loss, 0 neutral, +1 gain; NaN = unknown):
  sim1          spike-in truth (rep{N}.truth.tsv); event covers >= 50 % of the bin; negatives all 0
  hg008         NIST HG008-T draft somatic CNV benchmark; gain/loss relative to the modal copy number;
                bins with < 50 % benchmark coverage = NaN; pure normal (tumour fraction 0) = all 0
  prjna795275   all 0 (non-cancer plasma, assumed neutral)
  gse71378      consensus of the callers (1-Mb windows) = weak labels
Output: $MLO/hf_dataset/  (data/bins_100kb/<split>.csv.gz, data/bins_1Mb/<split>.csv.gz, samples.tsv, README.md)
An earlier hf_dataset folder is renamed hf_dataset.prev_<time>, never deleted.
"""
import os, re, sys, time, json
from pathlib import Path
import numpy as np, pandas as pd

E = os.environ
MLO = Path(E["MLO"])
AUTO = [str(i) for i in range(1, 23)]
nc = lambda c: re.sub(r"^chr", "", str(c))
RES = {"bins_100kb": 100000, "bins_1Mb": 1000000}
SPLITS = ["train", "validation", "test_sim1_ih03", "test_hg008", "test_noncancer", "test_cancer"]

man = pd.read_csv(MLO / "manifest.tsv", sep="\t", dtype=str)
for c in ("depth_x", "tumour_fraction_nominal", "tumour_fraction_realised"):
    man[c] = pd.to_numeric(man[c], errors="coerce")
cpath = lambda r: MLO / "counts_10kb" / r["dataset"] / f"{r['sample_id']}.bins10kb.tsv.gz"
have = [(MLO / "counts_10kb" / r["dataset"] / f"{r['sample_id']}.done").exists() and cpath(r).exists() for _, r in man.iterrows()]
missing = man[[not h for h in have]]
man = man[have].reset_index(drop=True)
print(f"{len(man)} samples with counts, {len(missing)} without (skipped)", flush=True)
if len(man) == 0:
    sys.exit("no count files yet - run step B first")

# ---------------------------------------------------------------- reference annotation + mappability proxy
ref = {}
for build in sorted(man.build.unique()):
    g = pd.read_csv(MLO / "reference" / f"{build}.bins10kb.tsv.gz", sep="\t", dtype={"chrom": str})
    s_all = np.zeros(len(g)); s_q = np.zeros(len(g))
    key = g.chrom + ":" + g.start.astype(str)
    for _, r in man[man.build == build].iterrows():
        d = pd.read_csv(cpath(r), sep="\t", dtype={"chrom": str}, usecols=["chrom", "start", "n_all", "n_q20"])
        d.index = d.chrom + ":" + d.start.astype(str)
        d = d.reindex(key)
        s_all += d.n_all.fillna(0).values; s_q += d.n_q20.fillna(0).values
    with np.errstate(invalid="ignore", divide="ignore"):
        g["map_proxy"] = np.where(s_all > 0, s_q / s_all, np.nan)
    g["usable"] = (g.n_frac <= 0.05) & (g.blacklist_frac <= 0.10) & (g.map_proxy >= 0.75)
    ref[build] = g
    print(build, f"{g.usable.mean():.1%} of 10-kb bins usable", flush=True)

# ---------------------------------------------------------------- truth sets
def intervals_sim1(path):
    t = pd.read_csv(path, sep="\t", dtype=str, comment=None)
    t.columns = [c.lstrip("#") for c in t.columns]
    return pd.DataFrame({"chrom": t["chrom"].map(nc), "start": t["start"].astype(int), "end": t["end"].astype(int),
                         "cn": t["copy_number"].astype(float), "event_id": t["cnv_id"],
                         "event_size": t["end"].astype(int) - t["start"].astype(int)})

hg = None
if (man.dataset == "hg008").any():
    t = pd.read_csv(E["TRUTH38"], sep="\t", header=None, comment="#", dtype={0: str})
    hg = pd.DataFrame({"chrom": t[0].map(nc), "start": t[1].astype(int), "end": t[2].astype(int), "cn": t[3].astype(float)})
    hg = hg[hg.chrom.isin(AUTO)]
    Lh = hg.end - hg.start
    HG_PLOIDY = float((hg.cn * Lh).sum() / Lh.sum())
    HG_MODE = float(hg.assign(L=Lh).groupby("cn").L.sum().idxmax())
    print(f"HG008 benchmark: mean ploidy {HG_PLOIDY:.3f}, modal CN {HG_MODE:g}", flush=True)

gse = {}
if (man.dataset == "gse71378").any() and Path(E["GSE_LABELS"]).exists():
    w = pd.read_csv(E["GSE_LABELS"], sep="\t", dtype={"chrom": str})
    w = w[(w.tool == "consensus") & (w.level == "full")]
    for b, g in w.groupby("sample"):
        gse[b] = {(nc(c), int(s)): int(v) for c, s, v in zip(g.chrom, g.start, g.state)}

_sim1_cache = {}


def overlap(bins, iv, B):
    """per bin: covered bp, cn-weighted bp, id of the interval covering most of the bin"""
    cov = np.zeros(len(bins)); cnw = np.zeros(len(bins)); best = np.zeros(len(bins)); bid = np.full(len(bins), None, object)
    pos = {(c, s): i for i, (c, s) in enumerate(zip(bins.chrom, bins.start))}
    for r in iv.itertuples():
        for b0 in range((r.start // B) * B, r.end, B):
            i = pos.get((r.chrom, b0))
            if i is None:
                continue
            o = min(r.end, b0 + B) - max(r.start, b0)
            if o <= 0:
                continue
            cov[i] += o; cnw[i] += o * r.cn
            if o > best[i]:
                best[i] = o; bid[i] = getattr(r, "event_id", None)
    return cov, cnw, bid


def labels(r, bins, B):
    n = len(bins)
    out = pd.DataFrame({"label_state": np.zeros(n), "label_cn": np.full(n, np.nan), "truth_log2": np.full(n, np.nan),
                        "event_id": np.full(n, None, object), "event_size": np.full(n, np.nan)})
    ds = r["dataset"]
    if ds == "sim1":
        out["label_cn"] = 2.0
        if r["class"] == "spiked":
            if r["truth"] not in _sim1_cache:
                _sim1_cache[r["truth"]] = intervals_sim1(r["truth"])
            iv = _sim1_cache[r["truth"]]
            cov, cnw, bid = overlap(bins, iv, B)
            hit = cov >= 0.5 * B
            cn = np.where(hit, np.round(cnw / np.maximum(cov, 1)), 2.0)
            out["label_cn"] = cn; out["label_state"] = np.sign(cn - 2)
            out.loc[hit, "event_id"] = bid[hit]
            sz = dict(zip(iv.event_id, iv.event_size))
            out["event_size"] = [sz.get(e, np.nan) if e is not None else np.nan for e in out.event_id]
        out["truth_log2"] = np.log2(out.label_cn / 2)
        out["label_source"] = "simulation_truth"
    elif ds == "hg008":
        cov, cnw, _ = overlap(bins, hg, B)
        ok = cov >= 0.5 * B
        cn = np.where(ok, cnw / np.maximum(cov, 1), np.nan)
        st = np.where(cn > HG_MODE + 0.5, 1.0, np.where(cn < HG_MODE - 0.5, -1.0, 0.0))
        st[~ok] = np.nan
        out["label_cn"] = np.round(cn, 3); out["truth_log2"] = np.log2(cn / HG_PLOIDY)
        tf = r["tumour_fraction_nominal"]
        out["label_state"] = 0.0 if (tf == 0) else st          # pure normal: neutral everywhere
        out["label_source"] = "NIST_HG008_draft_benchmark"
    elif ds == "prjna795275":
        out["label_source"] = "assumed_neutral_noncancer"
    else:
        lab = gse.get(r["base_sample"])
        if lab is None:
            out["label_state"] = np.nan; out["label_source"] = "no_consensus_available"
        else:
            out["label_state"] = [float(lab.get((c, (s // 1000000) * 1000000), 0)) for c, s in zip(bins.chrom, bins.start)]
            out["label_source"] = "caller_consensus_weak"
    return out


# ---------------------------------------------------------------- features
def roll(x, chrom, w, how="median"):
    s = pd.Series(x)
    f = (lambda v: v.rolling(w, center=True, min_periods=1).median()) if how == "median" else \
        (lambda v: v.rolling(w, center=True, min_periods=1).mean())
    return s.groupby(chrom.values).transform(f).values


def features(d10, g, B):
    d = d10.merge(g[["chrom", "start", "gc", "n_frac", "map_proxy", "usable"]], on=["chrom", "start"], how="left")
    d["usable"] = d["usable"].fillna(False).astype(bool)
    d["acgt"] = (d.end - d.start) * (1 - d.n_frac.fillna(1))
    d["gcw"] = d.gc.fillna(0) * d.acgt
    d["mapw"] = d.map_proxy.fillna(0) * (d.end - d.start)
    d["bin"] = (d.start // B) * B
    u = d[d.usable]
    a = u.groupby(["chrom", "bin"], sort=False).agg(
        n_q20=("n_q20", "sum"), n_all=("n_all", "sum"), n_frag=("n_frag", "sum"), n_short=("n_short", "sum"),
        n_long=("n_long", "sum"), sum_flen=("sum_flen", "sum"), acgt=("acgt", "sum"), gcw=("gcw", "sum"),
        mapw=("mapw", "sum"), usable_bp=("end", "sum"), start_bp=("start", "sum")).reset_index()
    a["usable_bp"] = a.usable_bp - a.start_bp
    a = a.drop(columns="start_bp").rename(columns={"bin": "start"})
    a["frac_usable"] = a.usable_bp / B
    a = a[a.frac_usable >= 0.5].copy()
    a["chrom"] = pd.Categorical(a.chrom, AUTO, ordered=True)
    a = a.sort_values(["chrom", "start"]).reset_index(drop=True)
    a["chrom"] = a.chrom.astype(str)
    a["end"] = a.start + B
    a["gc"] = a.gcw / a.acgt.clip(lower=1)
    a["map_proxy"] = a.mapw / a.usable_bp.clip(lower=1)
    c = a.n_q20 / a.frac_usable                                   # counts per full bin length
    q = pd.qcut(a.gc.rank(method="first"), 20, labels=False, duplicates="drop") if len(a) >= 40 else pd.Series(0, index=a.index)
    exp = c.groupby(q).transform("median")
    lr = np.log2((c + 0.5) / (exp + 0.5))
    lr = lr - np.nanmedian(lr)
    a["lr"] = lr
    a["lr_med3"] = roll(lr, a.chrom, 3); a["lr_med9"] = roll(lr, a.chrom, 9); a["lr_med25"] = roll(lr, a.chrom, 25)
    a["lr_mean9"] = roll(lr, a.chrom, 9, "mean")
    a["lr_chrom"] = pd.Series(lr).groupby(a.chrom.values).transform("median").values
    dl = np.diff(lr)
    noise = float(np.nanmedian(np.abs(dl - np.nanmedian(dl))) * 1.4826 / np.sqrt(2)) if len(dl) > 10 else np.nan
    a["noise"] = noise
    a["lr_z9"] = a.lr_med9 / (noise / 3 if noise and noise > 0 else np.nan)   # median of 9 ~ noise/3
    a["log_depth_bin"] = np.log10(np.nanmedian(c) + 1)
    if a.n_frag.sum() > 1000:
        sf = (a.n_short + 0.5) / (a.n_short + a.n_long + 1)
        a["short_frac_lr"] = np.log2(sf / np.nanmedian(sf))
        mf = a.sum_flen / a.n_frag.replace(0, np.nan)
        a["flen_dev"] = mf - np.nanmedian(mf)
        a["short_frac_med9"] = roll(a.short_frac_lr.values, a.chrom, 9)
    else:                                                          # single-end data: unknown, kept as NaN
        a["short_frac_lr"] = np.nan; a["flen_dev"] = np.nan; a["short_frac_med9"] = np.nan
    keep = ["chrom", "start", "end", "gc", "map_proxy", "frac_usable", "n_q20", "n_frag", "n_short", "n_long",
            "lr", "lr_med3", "lr_med9", "lr_med25", "lr_mean9", "lr_chrom", "noise", "lr_z9", "log_depth_bin",
            "short_frac_lr", "short_frac_med9", "flen_dev"]
    a = a[keep]
    for col in keep[3:]:
        if a[col].dtype.kind == "f":
            a[col] = a[col].round(5)
    return a


# ---------------------------------------------------------------- write
final = MLO / "hf_dataset"; tmpd = MLO / "hf_dataset.building"
if tmpd.exists():
    tmpd.rename(MLO / f"hf_dataset.failed_{time.strftime('%Y%m%d_%H%M%S')}")
for cfg in RES:
    (tmpd / "data" / cfg).mkdir(parents=True, exist_ok=True)
META = ["dataset", "split", "sample_id", "base_sample", "class", "replicate", "depth_x",
        "tumour_fraction_nominal", "tumour_fraction_realised", "build"]
started = set(); stats = []
for i, r in man.iterrows():
    d10 = pd.read_csv(cpath(r), sep="\t", dtype={"chrom": str})
    for cfg, B in RES.items():
        f = features(d10, ref[r["build"]], B)
        lab = labels(r, f, B)
        out = pd.concat([pd.DataFrame({k: [r[k]] * len(f) for k in META}), f.reset_index(drop=True), lab.reset_index(drop=True)], axis=1)
        fn = tmpd / "data" / cfg / f"{r['split']}.csv.gz"
        out.to_csv(fn, index=False, compression="gzip", mode="a", header=(fn not in started))
        started.add(fn)
        if cfg == "bins_100kb":
            ls = lab.label_state
            stats.append({**{k: r[k] for k in META}, "bins": len(f), "noise": f.noise.iloc[0] if len(f) else np.nan,
                          "bins_gain": int((ls == 1).sum()), "bins_loss": int((ls == -1).sum()),
                          "bins_unlabelled": int(ls.isna().sum()), "label_source": lab.label_source.iloc[0] if len(lab) else ""})
    if i % 10 == 0:
        print(f"[{i + 1}/{len(man)}] {r['sample_id']}", flush=True)

st = pd.DataFrame(stats)
st.to_csv(tmpd / "samples.tsv", sep="\t", index=False)
if len(missing):
    missing[["dataset", "sample_id", "split"]].to_csv(tmpd / "samples_without_counts.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- dataset card (Hugging Face README.md)
present = [s for s in SPLITS if (tmpd / "data" / "bins_100kb" / f"{s}.csv.gz").exists()]
yaml = ["---", "license: other", "pretty_name: cfCNV benchmark - bin-level cfDNA WGS copy-number dataset",
        "tags:", "- biology", "- genomics", "- copy-number-variation", "- cell-free-DNA", "configs:"]
for cfg in RES:
    yaml += [f"- config_name: {cfg}"] + (["  default: true"] if cfg == "bins_100kb" else []) + ["  data_files:"]
    for s in present:
        yaml += [f"  - split: {s}", f"    path: data/{cfg}/{s}.csv.gz"]
yaml.append("---")
def md_table(t):
    h = "| " + " | ".join(t.columns) + " |\n|" + "---|" * len(t.columns) + "\n"
    return h + "\n".join("| " + " | ".join(str(v) for v in row) + " |" for row in t.itertuples(index=False))


tab = st.groupby(["split", "dataset"]).agg(samples=("sample_id", "size"), depths=("depth_x", lambda x: ", ".join(f"{v:g}" for v in sorted(set(x.dropna())))),
                                            gain_bins=("bins_gain", "sum"), loss_bins=("bins_loss", "sum")).reset_index()
card = "\n".join(yaml) + f"""

# cfCNV benchmark dataset (bin level)

Bin-level read-depth and fragment-length features from low-coverage cell-free DNA whole-genome
sequencing, with copy-number labels, for training and testing CNV/CNA callers.
Built by `c_assemble.py` of the cfCNV workflow on {time.strftime('%Y-%m-%d')}.

**Status: internal draft. Before making this public, confirm that the data-use terms of every source
(GSE71378, PRJNA795275, GIAB HG008) allow redistribution of derived bin-level counts, and choose a licence
(`license: other` is a placeholder).**

## Splits (fixed before any model was trained)

| split | arm | labels |
|---|---|---|
| train | Simulation 1, healthy cfDNA backgrounds IH01 + IH02, replicates 1-4 + negative controls | simulated spike-ins (exact) |
| validation | Simulation 1, IH01 + IH02, replicate 5 (CNV layout not seen in train) | simulated spike-ins (exact) |
| test_sim1_ih03 | Simulation 1, unseen individual IH03 | simulated spike-ins (exact) |
| test_hg008 | HG008 tumour/normal mixtures (GRCh38) | NIST GIAB HG008-T draft somatic CNV benchmark |
| test_noncancer | PRJNA795275 non-cancer plasma (GRCh38) | assumed neutral everywhere |
| test_cancer | GSE71378 cancer plasma (GRCh37) | consensus of the benchmarked callers (weak labels) |

{md_table(tab)}

## Columns

* sample: `dataset split sample_id base_sample class replicate depth_x tumour_fraction_nominal tumour_fraction_realised build`
  (for test_cancer, `tumour_fraction_realised` holds the ichorCNA low-TF estimate of the full-depth sample)
* bin: `chrom` (1-22, no 'chr'), `start`, `end` (0-based, half-open, of the given build)
* annotation: `gc`, `map_proxy` (MAPQ>=20 share of fragments across all samples of the build), `frac_usable`
* counts: `n_q20` (fragments MAPQ>=20), `n_frag`, `n_short` (<=150 bp), `n_long` (>150 bp)
* features: `lr` (GC-corrected log2 ratio, sample median 0), `lr_med3/9/25` (rolling medians, bins), `lr_mean9`,
  `lr_chrom` (chromosome median), `noise` (sample-level), `lr_z9`, `log_depth_bin`, `short_frac_lr`,
  `short_frac_med9`, `flen_dev` (NaN for single-end data)
* labels: `label_state` (-1 loss, 0 neutral, +1 gain, empty = unknown), `label_cn`, `truth_log2`,
  `event_id`, `event_size` (Simulation 1), `label_source`

HG008: state is relative to the modal copy number of HG008-T ({'%g' % HG_MODE if hg is not None else 'NA'}); `truth_log2 = log2(CN / {('%.3f' % HG_PLOIDY) if hg is not None else 'NA'})`.
The label describes the tumour genome; the signal in a mixture scales with its tumour fraction.

## Loading

```python
from datasets import load_dataset
ds = load_dataset("path/or/hub-id", "bins_100kb")   # or "bins_1Mb"
```
"""
(tmpd / "README.md").write_text(card)
if final.exists():
    final.rename(MLO / f"hf_dataset.prev_{time.strftime('%Y%m%d_%H%M%S')}")
tmpd.rename(final)
print("\nDataset written:", final)
print(tab.to_string(index=False))
