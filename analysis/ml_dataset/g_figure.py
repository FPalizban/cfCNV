#!/usr/bin/env python
"""G  Manuscript figure: baseline models vs the 8 callers (reads only the tables written by f_report.py).
a  HG008 profile correlation vs tumour fraction at 1x      b  strict limit of detection vs depth
c  non-cancer plasma: % genome called altered (native)     d  cancer plasma: F1 vs consensus (native depth)
Runs on the cluster (step F) or on a laptop:   python g_figure.py <folder with the manuscript tables>
Writes figure_ml_benchmark.pdf / .png / .svg into that folder.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

T = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
COL = {"QDNAseq": "#1b9e77", "ichorCNA": "#d95f02", "BIC-seq2": "#7570b3", "CNVkit": "#e7298a", "Spectre": "#66a61e",
       "CNV-MEANN": "#e6ab02", "dudeML": "#a6761d", "LYCEUM": "#666666",
       "ML v1 gradient boosting": "#9e9e9e", "ML v2 gradient boosting": "#000000", "Signal (no model)": "#000000"}
STY = {"CNV-MEANN": "--", "dudeML": "--", "LYCEUM": "--", "ML v1 gradient boosting": "-", "ML v2 gradient boosting": "-",
       "Signal (no model)": ":"}
CALLERS = ["QDNAseq", "ichorCNA", "BIC-seq2", "CNVkit", "Spectre", "CNV-MEANN", "dudeML", "LYCEUM"]
OURS = ["ML v1 gradient boosting", "ML v2 gradient boosting"]
MM = 1 / 25.4
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 7,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6, "legend.frameon": False,
                     "legend.fontsize": 6, "pdf.fonttype": 42, "svg.fonttype": "none", "savefig.bbox": "tight", "savefig.dpi": 600})


def panel(ax, s):
    ax.text(-0.16, 1.06, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


fig, axs = plt.subplots(2, 2, figsize=(180 * MM, 130 * MM)); fig.subplots_adjust(wspace=0.35, hspace=0.55)

# a
ax = axs[0, 0]; r = pd.read_csv(T / "hg008_median_r_tf_depth.tsv", sep="\t")
dep = 1.0 if (r.depth_x == 1.0).any() else sorted(r.depth_x.unique())[len(r.depth_x.unique()) // 2]
for tl in CALLERS + OURS + ["Signal (no model)"]:
    g = r[(r.tool == tl) & (r.depth_x == dep)].sort_values("nominal_fraction")
    if len(g):
        ax.plot(g.nominal_fraction * 100, g.pearson_r, STY.get(tl, "-"), color=COL.get(tl, "#333"), marker="o", ms=2.5,
                lw=1.6 if tl in OURS else 1, label=tl)
ax.set_xscale("log"); ax.set_xlabel("Tumour fraction (%)"); ax.set_ylabel("Profile correlation (r)")
ax.set_title(f"HG008 mixtures, {dep:g}x", loc="left"); ax.axhline(0, color="#ccc", lw=0.5); panel(ax, "a")

# b
ax = axs[0, 1]; l = pd.read_csv(T / "hg008_lod_by_depth.tsv", sep="\t")
for tl in CALLERS + OURS:
    g = l[l.tool == tl].sort_values("depth_x")
    if len(g):
        ax.plot(g.depth_x, g["LoD_or_0.5"] * 100, STY.get(tl, "-"), color=COL.get(tl, "#333"), marker="o", ms=2.5,
                lw=1.6 if tl in OURS else 1, label=tl)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("Depth (x)"); ax.set_ylabel("Limit of detection (% tumour)")
ax.set_yticks([1, 2, 5, 10, 20, 50]); ax.set_yticklabels(["1", "2", "5", "10", "20", "n.d."]); ax.minorticks_off()
ax.set_title("Strict limit of detection", loc="left"); panel(ax, "b")

# c
ax = axs[1, 0]; n = pd.read_csv(T / "noncancer_pct_altered_per_sample.tsv", sep="\t"); n = n[n.level == "full"]
order = [t for t in ["QDNAseq", "ichorCNA", "BIC-seq2", "CNVkit", "Spectre"] + OURS if t in set(n.tool)]
rng = np.random.default_rng(1)
for i, tl in enumerate(order):
    v = n[n.tool == tl].FGA_pct.values
    ax.scatter(i + rng.uniform(-.18, .18, len(v)), v, s=6, color=COL.get(tl, "#333"), lw=0, alpha=.8)
    ax.plot([i - .3, i + .3], [np.median(v)] * 2, color="k", lw=1)
ax.set_yscale("symlog", linthresh=0.1); ax.set_ylim(0, 100); ax.set_xticks(range(len(order)))
ax.set_xticklabels([t.replace(" gradient boosting", "") for t in order], rotation=35, ha="right")
ax.set_ylabel("Genome called altered (%)"); ax.set_title("Non-cancer plasma (native depth)", loc="left"); panel(ax, "c")

# d
ax = axs[1, 1]; c = pd.read_csv(T / "table_cancer_f1_native_depth.tsv", sep="\t"); c = c[~c["sample"].str.startswith("IH")]
order = [t for t in ["QDNAseq", "BIC-seq2", "CNVkit", "ichorCNA", "Spectre"] + OURS if t in set(c.tool)]
for i, tl in enumerate(order):
    v = c[c.tool == tl].F1.dropna().values
    ax.scatter(i + rng.uniform(-.18, .18, len(v)), v, s=8, color=COL.get(tl, "#333"), lw=0, alpha=.85)
    if len(v):
        ax.plot([i - .3, i + .3], [np.median(v)] * 2, color="k", lw=1)
ax.set_ylim(0, 1.02); ax.set_xticks(range(len(order))); ax.set_xticklabels([t.replace(" gradient boosting", "") for t in order], rotation=35, ha="right")
ax.set_ylabel("F1 vs consensus"); ax.set_title("Cancer plasma (native depth)", loc="left"); panel(ax, "d")

h, lab = [], []
for a_ in (axs[0, 0], axs[0, 1]):
    for hh, ll in zip(*a_.get_legend_handles_labels()):
        if ll not in lab:
            h.append(hh); lab.append(ll)
fig.legend(h, lab, loc="lower center", ncol=6, bbox_to_anchor=(0.5, -0.07), fontsize=6)
for ext in ("pdf", "png", "svg"):
    fig.savefig(T / f"figure_ml_benchmark.{ext}")
print("wrote", T / "figure_ml_benchmark.[pdf|png|svg]")
