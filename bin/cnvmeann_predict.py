#!/usr/bin/env python3
"""Apply a CNV-MEANN model (cnvmeann_train.py) to one sample and merge bins into events: runs of
>= max(2, 50 kb / bin) consecutive bins of the same predicted class (as in the cfCNV paper).
usage: cnvmeann_predict.py --model model.joblib --features DIR --out sample.events.tsv"""
import argparse
import joblib
import numpy as np
import pandas as pd
from cnvmeann_train import load_dir


def events(x, cls, b):
    minb = max(2, int(np.ceil(50000 / b))); out = []
    for c in range(1, 23):
        m = x[:, 0] == c
        if not m.any(): continue
        bins = x[m, 1].astype(int); k = cls[m]; o = np.argsort(bins); bins, k = bins[o], k[o]; i = 0
        while i < len(k):
            j = i
            while j + 1 < len(k) and k[j + 1] == k[i] and bins[j + 1] == bins[j] + 1: j += 1
            if k[i] != 0 and j - i + 1 >= minb:
                out.append((str(c), (bins[i] - 1) * b, bins[j] * b, "gain" if k[i] == 1 else "loss"))
            i = j + 1
    return pd.DataFrame(out, columns=["chrom", "start", "end", "type"])


ap = argparse.ArgumentParser()
for k in ("--model", "--features", "--out"): ap.add_argument(k, required=True)
a = ap.parse_args()
M = joblib.load(a.model); b = M["bin"]
x = load_dir(a.features, b)
pr = M["net"].predict_proba(M["scaler"].transform(x[:, 2:6])) * M["weights"]
ev = events(x, (pr / pr.sum(1, keepdims=True)).argmax(1), b)
ev.to_csv(a.out, sep="\t", index=False)
print(f"{len(ev)} CNV-MEANN events")
