#!/usr/bin/env python3
"""CNV-MEANN network (Python port used in the cfCNV paper for Simulation 2): 4 features (rd, gc, rel, qua)
-> 15 tanh hidden units -> 3 classes (normal / gain / loss), min-max scaling to [-1, 1], balanced sampling
of 3000 normal / 1500 gain / 1500 loss bins per training sample, and prior correction at prediction time.
The original CNV-MEANN trains the same network in MATLAB (used for Simulation 1 in the paper).

usage: cnvmeann_train.py --bin 500000 --out model.joblib DIR1 DIR2 ...
  each DIR holds <sample>.chr<N>.bin<bin>.txt files (cnvmeann_features.py + cnvmeann_normalize.py, with --truth)
"""
import argparse
from pathlib import Path
import joblib
import numpy as np
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import MinMaxScaler

N_PER = {0: 3000, 1: 1500, 2: 1500}


def load_dir(d, b):
    xs = []
    for f in sorted(Path(d).glob(f"*.chr*.bin{b}.txt")):
        if f.stat().st_size:
            c = int(f.name.split(".chr")[1].split(".")[0])
            x = np.loadtxt(f, ndmin=2); xs.append(np.column_stack([np.full(len(x), c), x]))
    if not xs:
        raise SystemExit(f"no feature files for bin {b} in {d}")
    x = np.vstack(xs); x[:, 6] = np.where(x[:, 6] >= 2, 2, x[:, 6])
    return x


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("dirs", nargs="+")
    ap.add_argument("--bin", type=int, required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=20261023)
    a = ap.parse_args(); b = a.bin; rng = np.random.default_rng(a.seed + b)
    X, y, prior = [], [], np.zeros(3)
    for d in a.dirs:
        x = load_dir(d, b); lab = x[:, 6].astype(int); prior += np.bincount(lab, minlength=3)
        for k in (0, 1, 2):
            idx = np.where(lab == k)[0]
            if len(idx):
                p = rng.choice(idx, min(N_PER[k], len(idx)), replace=False); X.append(x[p, 2:6]); y.append(np.full(len(p), k))
    X, y = np.vstack(X), np.concatenate(y)
    if len(set(y)) < 3:
        raise SystemExit("training data must contain normal, gain and loss bins (give CNV-positive samples with truth)")
    prior /= prior.sum(); w = prior / (np.bincount(y, minlength=3) / len(y))
    sc = MinMaxScaler((-1, 1)).fit(X)
    net = MLPClassifier(hidden_layer_sizes=(15,), activation="tanh", max_iter=400, early_stopping=True,
                        n_iter_no_change=15, random_state=b).fit(sc.transform(X), y)
    joblib.dump({"scaler": sc, "net": net, "weights": w, "bin": b, "n_train_samples": len(a.dirs)}, a.out)
    print(f"CNV-MEANN bin {b}: trained on {len(y)} bins from {len(a.dirs)} samples, prior {np.round(prior, 3)}")


if __name__ == "__main__":
    main()
