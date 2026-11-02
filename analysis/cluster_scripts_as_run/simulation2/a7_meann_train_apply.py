#!/usr/bin/env python3
"""A7b - CNV-MEANN v2 for Simulation 2: train ONE network per bin size on ALL Simulation 1
positive samples (v2 features, all 5 replicates), apply it to every HG008 mixture.
Same network (4 -> 15 tanh -> 3), sampling (3000/1500/1500 bins per sample), prior
correction and bins->events rule (>= max(2, 50 kb / bin) consecutive bins) as Simulation 1.
No HG008 label is ever used for training.
  python a7_meann_train_apply.py --bin 500000"""
import argparse, os
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import MinMaxScaler

N_PER = {0: 3000, 1: 1500, 2: 1500}


def load(d, s, b):
    xs = []
    for c in range(1, 23):
        f = Path(d) / f"bin_{b}" / s / f"{s}.chr{c}.bin{b}.txt"
        if f.exists() and f.stat().st_size:
            x = np.loadtxt(f, ndmin=2); xs.append(np.column_stack([np.full(len(x), c), x]))
    x = np.vstack(xs); x[:, 6] = np.where(x[:, 6] >= 2, 2, x[:, 6]); return x


def events(x, cls, b, prefix="chr"):
    minb = max(2, int(np.ceil(50000 / b))); out = []
    for c in range(1, 23):
        m = x[:, 0] == c
        if not m.any(): continue
        bins = x[m, 1].astype(int); k = cls[m]; o = np.argsort(bins); bins, k = bins[o], k[o]; i = 0
        while i < len(k):
            j = i
            while j + 1 < len(k) and k[j + 1] == k[i] and bins[j + 1] == bins[j] + 1: j += 1
            if k[i] != 0 and j - i + 1 >= minb:
                out.append((f"{prefix}{c}", (bins[i] - 1) * b, bins[j] * b, "gain" if k[i] == 1 else "loss"))
            i = j + 1
    return pd.DataFrame(out, columns=["chrom", "start", "end", "type"])


ap = argparse.ArgumentParser(); ap.add_argument("--bin", type=int, required=True); a = ap.parse_args()
b = a.bin; rng = np.random.default_rng(20261023 + b)
m1 = pd.read_csv(Path(os.environ["P"]) / "00_metadata/simulation1_manifests/master.tsv", sep="\t")
X, y, prior = [], [], np.zeros(3)
for s in m1[m1["class"] == "positive"].sample_id:
    x = load(os.environ["SIM1_MEANN_V2"], s, b); lab = x[:, 6].astype(int)
    prior += np.bincount(lab, minlength=3)
    for k in (0, 1, 2):
        idx = np.where(lab == k)[0]
        if len(idx):
            p = rng.choice(idx, min(N_PER[k], len(idx)), replace=False); X.append(x[p, 2:6]); y.append(np.full(len(p), k))
X, y = np.vstack(X), np.concatenate(y); prior /= prior.sum(); w = prior / (np.bincount(y, minlength=3) / len(y))
sc = MinMaxScaler((-1, 1)).fit(X)
net = MLPClassifier(hidden_layer_sizes=(15,), activation="tanh", max_iter=400, early_stopping=True,
                    n_iter_no_change=15, random_state=b).fit(sc.transform(X), y)
print(f"bin {b}: trained on {len(y)} Sim1 bins, prior {np.round(prior, 3)}", flush=True)
S2 = Path(os.environ["S2_OUT"]) / "CNV-MEANN"; out = S2 / "calls" / f"bin_{b}"; out.mkdir(parents=True, exist_ok=True)
log = []
for d in sorted((S2 / "features_v2" / f"bin_{b}").iterdir()):
    if not (d / ".done").exists(): continue
    x = load(S2 / "features_v2", d.name, b)
    pr = net.predict_proba(sc.transform(x[:, 2:6])) * w
    ev = events(x, (pr / pr.sum(1, keepdims=True)).argmax(1), b)
    ev.to_csv(out / f"{d.name}.events.tsv", sep="\t", index=False); log.append((d.name, len(ev)))
pd.DataFrame(log, columns=["sample", "events"]).to_csv(out / "prediction_log.tsv", sep="\t", index=False)
(out / ".done").write_text("ok\n"); print("A7B_DONE bin", b, len(log), "samples")
