#!/usr/bin/env python3
"""CNV-MEANN v2 normalisation for ONE sample x bin (identical rules to Simulation 1, step B1):
rd -> GC-stratum median correction (50 strata) -> / autosomal median ; rel re-computed on the
normalised rd (5 bins each side) ; gc, qua, label unchanged ; same 6-column layout."""
import argparse
from pathlib import Path
import numpy as np


def rel5(rd):
    n = len(rd); out = np.zeros(n); cs = np.concatenate([[0.0], np.cumsum(rd)])
    for i in range(n):
        lo, hi = max(0, i - 5), min(n, i + 6)
        l = (cs[i] - cs[lo]) / (i - lo) if i > lo else None
        r = (cs[hi] - cs[i + 1]) / (hi - i - 1) if hi > i + 1 else None
        if l is not None and r is not None: out[i] = abs(((rd[i] - l) + (rd[i] - r)) / 2)
        elif l is not None: out[i] = abs(rd[i] - l)
        elif r is not None: out[i] = abs(rd[i] - r)
    return out


ap = argparse.ArgumentParser()
for k in ("--raw", "--out", "--sample"): ap.add_argument(k, required=True)
ap.add_argument("--bin", type=int, required=True); a = ap.parse_args()
raw, out = Path(a.raw), Path(a.out); out.mkdir(parents=True, exist_ok=True)
per = {}
for c in range(1, 23):
    f = raw / f"{a.sample}.chr{c}.bin{a.bin}.txt"
    per[c] = np.loadtxt(f, ndmin=2) if f.exists() and f.stat().st_size else np.zeros((0, 6))
allx = np.vstack([v for v in per.values() if len(v)])
rd, gc = allx[:, 1], allx[:, 2]
ok = np.isfinite(rd) & (gc > 0.2) & (gc < 0.8)
edges = np.quantile(gc[ok], np.linspace(0, 1, 51))
glob = np.median(rd[ok & (rd > 0)]) if np.any(ok & (rd > 0)) else 1.0
idx = np.clip(np.searchsorted(edges, gc, side="right") - 1, 0, 49); med = np.ones(50)
for k in range(50):
    m = ok & (idx == k); v = np.median(rd[m]) if m.sum() >= 20 else 0
    med[k] = v if v > 0 else glob
for c, x in per.items():
    if len(x):
        x[:, 1] = x[:, 1] / med[np.clip(np.searchsorted(edges, x[:, 2], side="right") - 1, 0, 49)]
m0 = np.median(np.concatenate([x[:, 1] for x in per.values() if len(x)])); m0 = m0 if m0 > 0 else 1.0
for c, x in per.items():
    if len(x):
        x[:, 1] /= m0; x[:, 3] = rel5(x[:, 1]); x[:, 4] = np.nan_to_num(x[:, 4], nan=0.0)
    np.savetxt(out / f"{a.sample}.chr{c}.bin{a.bin}.txt", x, delimiter="\t",
               fmt=["%d", "%.5f", "%.5f", "%.5f", "%.5f", "%d"] if len(x) else "%g")
(out / ".done").write_text("ok\n")
