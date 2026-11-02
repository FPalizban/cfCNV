"""Figure 5 - computational cost and summary across benchmark layers.
a  median wall time per run (bars) with 90th percentile (whisker), Simulation 1 (ACE: Simulation 2)
b  peak memory per run: median (bar) and maximum (whisker)
c  summary scorecard computed from the data files (colour = rank within column), including agreement with
   the other callers in cancer plasma (median leave-one-out F1 over 5 samples)
Run:  python fig5_cost_summary.py"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import style as S, common as C
S.apply()
R = pd.read_csv(C.D / "runtime_memory_summary.tsv", sep="\t")
R = R[(R.sim == "Sim1") | (R.tool == "ACE")].set_index("tool")
order = [t for t in S.TOOL_ORDER + ["ACE"] if t in R.index]
fig = plt.figure(figsize=(S.FULL_WIDTH, 150 * S.MM)); gs = fig.add_gridspec(2, 2, height_ratios=[0.8, 1.2], hspace=0.55, wspace=0.45); src = []
for k, (col, hi, lab, log) in enumerate([("med_min", "p90_min", "Wall time per run (min)", True), ("med_gb", "max_gb", "Peak memory per run (GB)", False)]):
    ax = fig.add_subplot(gs[0, k]); y = np.arange(len(order))
    ax.barh(y, R.loc[order, col], color=[S.color(t) for t in order]); ax.errorbar(R.loc[order, col], y, xerr=[np.zeros(len(order)), R.loc[order, hi] - R.loc[order, col]], fmt="none", ecolor="k", elinewidth=0.6, capsize=1.5)
    ax.set_yticks(y); ax.set_yticklabels(order); ax.invert_yaxis(); ax.set_xlabel(lab)
    if log: ax.set_xscale("log")
    S.panel(ax, "ab"[k], -0.3)
    src += [{"panel": "ab"[k], "tool": t, "median": R.loc[t, col], "upper": R.loc[t, hi]} for t in order]

# c scorecard from data
pc, td, prim = C.sim1(); pr, tf, mix, p2 = C.sim2(); B, *_ = C.healthy(); calls = C.healthy()[3]
pos = pc[(pc["class"] == "positive") & (pc.depth <= 1)]; pos = pos[[prim.get(t) == b for t, b in zip(pos.tool, pos.bin_bp)]]
P = C.primary(pr, p2); L = C.primary(C.lod_strict(pr), p2)
H = dict(C.HEALTHY_PRIMARY)
LC = pd.read_csv(C.D / "cancer_loo_and_tcga.tsv", sep="\t"); LC = LC[LC["sample"] != "IH03"]
tools = S.TOOL_ORDER; rows = []
for t in tools:
    f1 = pos[pos.tool == t].F1.median(); r5 = P[(P.tool == t) & (P.nominal_fraction == 0.05)].pearson_r.median()
    lod = L[(L.tool == t) & (L.depth_x.isin([0.25, 0.5, 1]))].LoD; lod = lod.max() if lod.notna().all() and len(lod) else np.nan
    if t in H:
        fga = B[(B.tool == t) & (B.config == H[t]) & (B.level == "full")].FGA_pct.median()
        c = calls[(calls.tool == t) & (calls.config == H[t]) & (calls.level == "full")]; big = c[(c.end - c.start) >= 10e6].base.nunique()
    else: fga, big = np.nan, np.nan
    cf = LC[LC.tool == t].F1.median() if t in set(LC.tool) else np.nan
    rt = R.loc[t, "med_min"] if t in R.index else np.nan
    rows.append([f1, r5, lod, fga, big, cf, rt])
V = np.array(rows, float)
better_high = [True, True, False, False, False, True, False]
G = np.full(V.shape, np.nan)
for j in range(V.shape[1]):
    col = V[:, j]; ok = ~np.isnan(col)
    if j == 2: col = np.where(np.isnan(col), 0.5, col); ok = np.ones(len(col), bool)   # not detected -> worst
    rk = pd.Series(col[ok]).rank(pct=True, ascending=better_high[j]).values; G[ok, j] = rk
cm = LinearSegmentedColormap.from_list("score", S.SCORE_CMAP)
fmt = [lambda v: f"{v:.2f}", lambda v: f"{abs(v):.2f}" if abs(v) < 0.005 else f"{v:.2f}", lambda v: "n.d." if np.isnan(v) else f"{v*100:g}%", lambda v: "not run" if np.isnan(v) else f"{v:.2g}%",
       lambda v: "not run" if np.isnan(v) else f"{int(v)}/25", lambda v: "not run" if np.isnan(v) else f"{v:.2f}", lambda v: f"{v:.1f}" if v < 10 else f"{v:.0f}"]
cols = ["Sim 1 event F1\n(0.1–1×)", "HG008 profile r\nat 5% tumour", "HG008 LoD\n(0.25–1×, worst)", "Non-cancer:\ngenome altered", "Non-cancer:\ncalls ≥10 Mb", "Cancer: F1 vs\nother callers", "Runtime per\nrun (min)"]
ax = fig.add_subplot(gs[1, :])
for i, t in enumerate(tools):
    for j in range(len(cols)):
        g = G[i, j]; ax.add_patch(plt.Rectangle((j, i), 0.96, 0.9, color="#e5e5e5" if np.isnan(g) else cm(g), lw=0))
        ax.text(j + .48, i + .45, fmt[j](V[i, j]), ha="center", va="center", fontsize=6.5)
    src += [{"panel": "c", "tool": t, "column": c.replace("\n", " "), "value": V[i, j]} for j, c in enumerate(cols)]
ax.set_xlim(0, len(cols)); ax.set_ylim(len(tools), 0); ax.set_xticks(np.arange(len(cols)) + .48); ax.set_xticklabels(cols, fontsize=5.6); ax.xaxis.tick_top()
ax.set_yticks(np.arange(len(tools)) + .45); ax.set_yticklabels(tools); ax.tick_params(length=0)
for s_ in ax.spines.values(): s_.set_visible(False)
S.panel(ax, "c", -0.1, 1.12)
S.save(fig, "Figure5_cost_and_summary"); pd.DataFrame(src).to_csv("output/Figure5_source_data.tsv", sep="\t", index=False)
