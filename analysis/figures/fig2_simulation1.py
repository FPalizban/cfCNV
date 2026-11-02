"""Figure 2 - Simulation 1 (inserted events in healthy cfDNA).
a  median precision / recall / F1 per caller (primary bin, 0.1-1x, n = 60 positives)
b  false calls per Gb in unmodified negative controls (0.1-1x, n = 12)
c  median F1 by depth, one panel per background (IH01, IH02, IH03)
d  recall by event size (primary bin, 0.1-1x)
e  median F1 by bin size (heatmap, 0.1-1x)
Run:  python fig2_simulation1.py   -> output/Figure2_simulation1.[pdf|png|svg] + source-data tsv"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import style as S, common as C
S.apply()
pc, td, prim = C.sim1()
tools = [t for t in S.TOOL_ORDER if t in prim]
isprim = lambda d: d[[prim.get(t) == b for t, b in zip(d.tool, d.bin_bp)]]
P = isprim(pc); pos = P[(P["class"] == "positive")]; neg = P[(P["class"] == "negative") & (P.depth <= 1)]
m = pos[pos.depth <= 1]

fig = plt.figure(figsize=(S.FULL_WIDTH, 175 * S.MM))
gs = fig.add_gridspec(3, 6, height_ratios=[1, 0.9, 1], hspace=0.75, wspace=1.1)
src = []

# a
ax = fig.add_subplot(gs[0, :4]); x = np.arange(len(tools)); w = 0.26
for k, (met, alpha) in enumerate([("precision", 0.45), ("recall", 0.7), ("F1", 1.0)]):
    v = [m[m.tool == t][met].median() for t in tools]
    ax.bar(x + (k - 1) * w, v, w, color=[S.color(t) for t in tools], alpha=alpha, label=met if k else met)
    src += [{"panel": "a", "tool": t, "metric": met, "value": vv} for t, vv in zip(tools, v)]
ax.set_xticks(x); ax.set_xticklabels(tools, rotation=30, ha="right"); ax.set_ylim(0, 1.05); ax.set_ylabel("Median (0.1–1×)")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color="#555555", alpha=a, label=l) for l, a in [("Precision", .45), ("Recall", .7), ("F1", 1)]], ncol=3, loc="upper right")
S.panel(ax, "a", -0.08)

# b
ax = fig.add_subplot(gs[0, 4:]); rng = np.random.default_rng(1)
for i, t in enumerate(tools):
    v = neg[neg.tool == t].false_calls_per_Gb.values
    ax.scatter(np.full(len(v), i) + rng.uniform(-.2, .2, len(v)), v, s=7, color=S.color(t), lw=0, alpha=.8)
    ax.plot([i - .3, i + .3], [np.median(v)] * 2, color="k", lw=1)
    src += [{"panel": "b", "tool": t, "metric": "false_calls_per_Gb", "value": vv} for vv in v]
ax.set_yscale("symlog", linthresh=1); ax.set_xticks(range(len(tools))); ax.set_xticklabels(tools, rotation=60, ha="right")
ax.set_ylabel("False calls per Gb\n(negative controls)"); S.panel(ax, "b", -0.3)

# c
for k, bg in enumerate(["IH01", "IH02", "IH03"]):
    ax = fig.add_subplot(gs[1, 2 * k:2 * k + 2])
    for t in tools:
        g = pos[(pos.tool == t) & (pos.background == bg)].groupby("depth").F1.median()
        ax.plot(g.index, g.values, color=S.color(t), ls="--" if t in S.ML_TOOLS else "-", marker="o", ms=2.5, label=t)
        src += [{"panel": "c", "tool": t, "background": bg, "depth": d, "metric": "median_F1", "value": v} for d, v in g.items()]
    ax.set_xscale("log"); ax.set_xticks([0.1, 1, 10]); ax.set_xticklabels(["0.1", "1", "10"]); ax.minorticks_off()
    ax.set_ylim(-0.02, 0.85); ax.set_title(bg, loc="left"); ax.set_xlabel("Depth (×)")
    if k == 0: ax.set_ylabel("Median event F1"); S.panel(ax, "c", -0.3)
    if k == 2: ax.legend(fontsize=5, loc="center left", bbox_to_anchor=(1.0, 0.5))

# d
ax = fig.add_subplot(gs[2, :3]); T = td[(td.depth <= 1)]; T = T[[prim.get(t) == b for t, b in zip(T.tool, T.bin_bp)]]
sizes = sorted(T.size_bp.unique())
for t in tools:
    g = T[T.tool == t].groupby("size_bp").detected.mean().reindex(sizes)
    ax.plot(range(len(sizes)), g.values, color=S.color(t), ls="--" if t in S.ML_TOOLS else "-", marker="o", ms=2.5)
    src += [{"panel": "d", "tool": t, "size_bp": s, "metric": "recall", "value": v} for s, v in g.items()]
ax.set_xticks(range(len(sizes))); ax.set_xticklabels([f"{s/1e6:g}M" if s >= 1e6 else f"{s/1e3:g}k" for s in sizes])
ax.set_xlabel("Event size (bp)"); ax.set_ylabel("Recall (0.1–1×)"); ax.set_ylim(-0.03, 1.03); S.panel(ax, "d", -0.15)

# e
ax = fig.add_subplot(gs[2, 3:]); bins = [10000, 30000, 50000, 100000, 500000, 1000000]
ht = [t for t in tools if t != "LYCEUM"]
M = np.array([[pc[(pc.tool == t) & (pc.bin_bp == b) & (pc["class"] == "positive") & (pc.depth <= 1)].F1.median() for b in bins] for t in ht])
im = ax.imshow(M, cmap=S.HEATMAP_CMAP, vmin=0, vmax=0.75, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        ax.text(j, i, "–" if np.isnan(M[i, j]) else f"{M[i, j]:.2f}", ha="center", va="center", fontsize=5.5, color="white" if M[i, j] < 0.45 else "black")
        if prim.get(ht[i]) == bins[j]: ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False, lw=1.2, ec="k"))
src += [{"panel": "e", "tool": t, "bin_bp": b, "metric": "median_F1", "value": M[i, j]} for i, t in enumerate(ht) for j, b in enumerate(bins)]
ax.set_xticks(range(len(bins))); ax.set_xticklabels(["10k", "30k", "50k", "100k", "500k", "1M"]); ax.set_yticks(range(len(ht))); ax.set_yticklabels(ht)
ax.set_xlabel("Bin size"); cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label("Median F1"); S.panel(ax, "e", -0.3)
S.save(fig, "Figure2_simulation1"); pd.DataFrame(src).to_csv("output/Figure2_source_data.tsv", sep="\t", index=False)
