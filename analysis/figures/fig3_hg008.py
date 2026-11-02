"""Figure 3 - NIST HG008 tumour-normal mixtures (Simulation 2).
a  realised vs nominal tumour fraction (all 104 mixtures)
b  median profile correlation r vs tumour fraction, primary configuration, pooled over depths (IQR band)
c  strict limit of detection by depth (classical callers)
d  median r by depth x tumour fraction for QDNAseq and ichorCNA (heatmaps)
e  ichorCNA tumour-fraction estimates (tumour defaults vs low-TF settings) and ACE vs realised fraction
Run:  python fig3_hg008.py"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import style as S, common as C
S.apply()
pr, tf, mix, prim = C.sim2()
P = C.primary(pr, prim)
tools = [t for t in S.TOOL_ORDER if t in prim.index]
fr = sorted(P.nominal_fraction.unique()); deps = sorted(P.depth_x.unique())
fig = plt.figure(figsize=(S.FULL_WIDTH, 150 * S.MM)); gs = fig.add_gridspec(2, 4, hspace=0.6, wspace=0.85); src = []

# a
ax = fig.add_subplot(gs[0, 0])
ax.plot([0, 0.32], [0, 0.32], color="#999999", lw=0.8)
ax.scatter(mix.tumor_fraction / 100 if mix.tumor_fraction.max() > 1 else mix.tumor_fraction, mix.realised_fraction, s=8, color="#333333", lw=0)
ax.set_xlabel("Nominal tumour fraction"); ax.set_ylabel("Realised tumour fraction"); ax.set_title("104 mixtures", loc="left"); S.panel(ax, "a", -0.35)
src += [{"panel": "a", "sample": s, "nominal": n, "realised": r} for s, n, r in zip(mix.sample_id, mix.tumor_fraction, mix.realised_fraction)]

# b
ax = fig.add_subplot(gs[0, 1:3])
for t in tools:
    g = P[P.tool == t].groupby("nominal_fraction").pearson_r
    med, q1, q3 = g.median(), g.quantile(.25), g.quantile(.75)
    ax.plot(med.index * 100, med.values, color=S.color(t), ls="--" if t in S.ML_TOOLS else "-", marker="o", ms=2.5, label=t)
    ax.fill_between(med.index * 100, q1.values, q3.values, color=S.color(t), alpha=0.12, lw=0)
    src += [{"panel": "b", "tool": t, "fraction": f, "median_r": v} for f, v in med.items()]
ax.set_xscale("symlog", linthresh=1); ax.set_xticks([0, 1, 2, 5, 10, 20, 30]); ax.set_xticklabels(["0", "1", "2", "5", "10", "20", "30"])
ax.set_xlabel("Tumour DNA (%)"); ax.set_ylabel("Profile correlation r"); ax.set_ylim(-0.15, 1); ax.legend(fontsize=5, ncol=2, loc="upper left"); S.panel(ax, "b", -0.12)

# d (LoD)
L = C.primary(C.lod_strict(pr), prim)
ax = fig.add_subplot(gs[0, 3])
for k, t in enumerate(["QDNAseq", "ichorCNA", "BIC-seq2", "CNVkit"]):
    g = L[L.tool == t].set_index("depth_x").LoD.reindex(deps)
    ax.plot(deps, g.values * 100 * (1 + (k - 1.5) * 0.04), color=S.color(t), marker="o", ms=3, label=t)
    src += [{"panel": "c", "tool": t, "depth": d, "LoD": v} for d, v in g.items()]
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xticks([0.1, 1, 10]); ax.set_xticklabels(["0.1", "1", "10"]); ax.set_yticks([1, 2, 5, 10]); ax.set_yticklabels(["1", "2", "5", "10"]); ax.minorticks_off()
ax.set_xlabel("Depth (×)"); ax.set_ylabel("Limit of detection (% tumour)"); ax.legend(fontsize=5); S.panel(ax, "c", -0.4)

# c heatmaps
bot = gs[1, :].subgridspec(1, 5, width_ratios=[1, 1, 0.07, 0.55, 2.3], wspace=0.35)
for k, t in enumerate(["QDNAseq", "ichorCNA"]):
    ax = fig.add_subplot(bot[0, k]); f2 = [f for f in fr if f > 0]
    M = np.array([[P[(P.tool == t) & (P.nominal_fraction == f) & (P.depth_x == d)].pearson_r.median() for d in deps] for f in f2])
    im = ax.imshow(M, cmap=S.HEATMAP_CMAP, vmin=0, vmax=1, origin="lower", aspect="auto")
    ax.set_xticks(range(len(deps))); ax.set_xticklabels([f"{d:g}" for d in deps], fontsize=5, rotation=90); ax.set_yticks(range(len(f2))); ax.set_yticklabels([f"{f*100:g}" for f in f2])
    ax.set_xlabel("Depth (×)"); ax.set_ylabel("Tumour DNA (%)"); ax.set_title(t, loc="left")
    src += [{"panel": "d", "tool": t, "fraction": f, "depth": d, "median_r": M[i, j]} for i, f in enumerate(f2) for j, d in enumerate(deps)]
    if k == 0: S.panel(ax, "d", -0.6)
    else: ax.set_ylabel("")
cb = fig.colorbar(im, cax=fig.add_subplot(bot[0, 2])); cb.set_label("Median r")

# e TF estimates
ax = fig.add_subplot(bot[0, 4])
T = tf.copy(); T["x"] = T.realised_fraction
for lab, (tool, cfg), col, mk in [("ichorCNA, tumour defaults", ("ichorCNA", "1000kb_default"), S.color("ichorCNA"), "s"),
                                  ("ichorCNA, low-TF settings", ("ichorCNA", "1000kb_lowTF"), "#000000", "o"),
                                  ("ACE (100 kb, ploidy 2)", ("ACE", "ploidy2_100kb"), S.color("ACE"), "^")]:
    g = T[(T.tool == tool) & (T.config == cfg)]
    if g.empty: continue
    s = g.groupby("nominal_fraction").agg(x=("realised_fraction", "median"), y=("estimated_fraction", "median"), q1=("estimated_fraction", lambda v: v.quantile(.25)), q3=("estimated_fraction", lambda v: v.quantile(.75)))
    ax.errorbar(s.x, s.y, yerr=[s.y - s.q1, s.q3 - s.y], color=col, marker=mk, ms=3.5, lw=1, capsize=1.5, label=lab)
    src += [{"panel": "e", "method": lab, "realised": a, "estimate_median": b} for a, b in zip(s.x, s.y)]
xx = np.linspace(0, 0.31, 50); ax.plot(xx, xx, color="#999999", lw=0.8, label="identity (DNA fraction)")
ax.plot(xx, 2 * xx / (1.473 * (1 - xx) + 2 * xx), color="#999999", lw=0.8, ls=":", label="cellular fraction (ploidy 1.47)")
ax.set_ylim(-0.01, 0.45); ax.set_xlabel("Realised tumour DNA fraction"); ax.set_ylabel("Estimated tumour fraction"); ax.legend(fontsize=5, loc="upper left"); S.panel(ax, "e", -0.15)
S.save(fig, "Figure3_hg008_mixtures"); pd.DataFrame(src).to_csv("output/Figure3_source_data.tsv", sep="\t", index=False)
