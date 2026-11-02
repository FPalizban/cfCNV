#!/usr/bin/env python3
"""
cfCNV Simulation 1 - Figure S4: background-specific artefact and self-normal masking.

  python 05_mask_figure.py --bundle <cfCNV_Sim1_results_for_local> --mask mask/ --out figures/

(a) Calls in unmodified negative controls (QDNAseq, 50 kb) by depth and background.
(b) Fraction of those calls that are losses.
(c) Event F1 per sample before and after masking regions the same tool called in
    that background's own negative controls (leave-one-depth-out), primary bins.
"""
import argparse
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

here = Path(__file__).parent
spec = importlib.util.spec_from_file_location("f3", here / "03_figures.py")
f3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f3)
MM, COLOR, MARK, INK, MUTED = f3.MM, f3.COLOR, f3.MARK, f3.INK, f3.MUTED
BG_MARK = {"IH01": "o", "IH02": "s", "IH03": "^"}
BG_SHADE = {"IH01": "#1a1a1a", "IH02": "#7a7a7a", "IH03": "#b5b5b5"}
PRIMARY = {"QDNAseq": 50000, "CNVkit": 100000, "ichorCNA": 1000000}  # = stats/primary_configs.tsv (v2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", required=True); ap.add_argument("--mask", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    f3.style()
    O = Path(a.out); O.mkdir(parents=True, exist_ok=True)
    ch = pd.read_csv(Path(a.mask) / "negative_call_characteristics.tsv", sep="\t")
    rs = pd.read_csv(Path(a.mask) / "masked_rescore_per_sample.tsv", sep="\t")
    q = ch[(ch.tool == "QDNAseq") & (ch.analysis_bin_bp == 50000)]

    fig = plt.figure(figsize=(180 * MM, 62 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.5], wspace=0.42)
    for k, (col, yl) in enumerate([("n_calls", "Calls in negative control"),
                                   ("loss_fraction", "Fraction of calls that are losses")]):
        ax = fig.add_subplot(gs[k])
        for bg in ("IH01", "IH02", "IH03"):
            g = q[q.background == bg].sort_values("depth_x")
            ax.plot(g.depth_x, g[col], color=BG_SHADE[bg], marker=BG_MARK[bg], ms=4, label=bg)
        ds = sorted(q.depth_x.unique())
        ax.set_xscale("log"); ax.set_xticks(ds); ax.set_xticklabels([f"{d:g}" for d in ds], rotation=45)
        ax.minorticks_off(); ax.set_xlabel("Sequencing depth (×)"); ax.set_ylabel(yl)
        if col == "loss_fraction":
            ax.set_ylim(-0.03, 1.03)
        else:
            ax.legend(loc="upper left")
        ax.set_title("QDNAseq 50 kb", loc="left"); f3.plab(ax, "ab"[k], -0.3)

    ax = fig.add_subplot(gs[2])
    rows = []
    x0 = 0
    ticks, labs = [], []
    for t, b in PRIMARY.items():
        for bg in ("IH01", "IH02", "IH03"):
            g = rs[(rs.tool == t) & (rs.bin_bp == b) & (rs.background == bg) & (rs.depth <= 1)]
            u, m = g.F1_unmasked.median(), g.F1_masked_lodo.median()
            ax.plot([x0, x0], [u, m], color=COLOR[t], lw=1.6, zorder=2)
            ax.scatter([x0], [u], s=26, facecolor="white", edgecolor=COLOR[t], marker=MARK[t], lw=1, zorder=3)
            ax.scatter([x0], [m], s=26, color=COLOR[t], marker=MARK[t], zorder=3)
            rows.append({"tool": t, "bin_bp": b, "background": bg, "n": len(g),
                         "F1_unmasked_median": u, "F1_masked_lodo_median": m})
            ticks.append(x0); labs.append(bg); x0 += 1
        x0 += 0.8
    ax.set_xticks(ticks); ax.set_xticklabels(labs, rotation=45, ha="right")
    for i, t in enumerate(PRIMARY):
        ax.text(i * 3.8 + 1, 1.04, f3.lab(t), ha="center", fontsize=6.5, color=COLOR[t],
                transform=ax.get_xaxis_transform())
    ax.set_ylim(-0.03, 1.0); ax.set_ylabel("Median event F1 (0.1–1×)"); ax.grid(axis="x", visible=False)
    h = [plt.Line2D([], [], color=MUTED, marker="o", mfc="white", ls="", ms=4),
         plt.Line2D([], [], color=MUTED, marker="o", ls="", ms=4)]
    ax.legend(h, ["unmasked", "self-normal mask"], loc="upper right", bbox_to_anchor=(1, 0.93))
    f3.plab(ax, "c", -0.18)
    f3.save(fig, O, "FigS4_background_artifact_masking", pd.DataFrame(rows))


if __name__ == "__main__":
    main()
