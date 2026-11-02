"""Figure 4 - real plasma.
a  % of the genome called altered in 25 non-cancer plasma samples, native depth (primary configurations)
b  median % altered at 0.1x, 1x and native depth
c  tumour fraction reported at 0% tumour (ichorCNA low-TF, ACE 100 kb ploidy 2, Fragle), native depth
d  cancer plasma (GSE71378): copy-number state per 1-Mb window for each caller and the majority consensus (native depth)
e  correlation of each caller's profile with TCGA gain/loss frequencies of the matching cancer type,
   cancer samples (coloured) vs the 25 non-cancer samples (grey box) scored against the same cohorts
f  agreement of the 1x and 0.1x profiles with the caller's own native-depth profile (F1 over altered windows)
Run:  python fig4_real_plasma.py"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import style as S, common as C
S.apply()
B, TF, INS, CALLS, REC, ANN = C.healthy()
CF = C.HEALTHY_PRIMARY
fig = plt.figure(figsize=(S.FULL_WIDTH, 225 * S.MM)); gs = fig.add_gridspec(3, 3, width_ratios=[1.35, 1, 1], height_ratios=[1, 2.0, 0.95], wspace=0.45, hspace=0.6); src = []; rng = np.random.default_rng(1)

ax = fig.add_subplot(gs[0, 0])
for i, (t, c) in enumerate(CF):
    v = B[(B.tool == t) & (B.config == c) & (B.level == "full")].FGA_pct.values
    ax.scatter(i + rng.uniform(-.18, .18, len(v)), v, s=7, color=S.color(t), lw=0, alpha=.8); ax.plot([i - .3, i + .3], [np.median(v)] * 2, color="k", lw=1)
    src += [{"panel": "a", "tool": t, "config": c, "pct_genome_altered": x} for x in v]
ax.set_yscale("symlog", linthresh=0.1); ax.set_ylim(0, 100); ax.set_xticks(range(len(CF))); ax.set_xticklabels([f"{t}\n{c}" for t, c in CF], fontsize=5.5, rotation=35, ha="right")
ax.set_ylabel("Genome called altered (%)"); ax.set_title("25 individuals without cancer", loc="left"); S.panel(ax, "a", -0.18)

ax = fig.add_subplot(gs[0, 1]); lv = [("d0p1", 0.1), ("d1", 1), ("full", 4)]
for t, c in CF + [("QDNAseq", "50kb")]:
    m = [B[(B.tool == t) & (B.config == c) & (B.level == l)].FGA_pct.median() for l, _ in lv]
    ax.plot([x for _, x in lv], m, color=S.color(t), marker="o", ms=3, ls=":" if c == "50kb" else "-", label=f"{t} {c}")
    src += [{"panel": "b", "tool": t, "config": c, "level": l, "median_pct_altered": v} for (l, _), v in zip(lv, m)]
ax.set_xscale("log"); ax.set_xticks([0.1, 1, 4]); ax.set_xticklabels(["0.1", "1", "native"]); ax.minorticks_off()
ax.set_yscale("symlog", linthresh=0.1); ax.set_ylim(0, 100); ax.set_xlabel("Depth (×)"); ax.set_ylabel("Median genome altered (%)"); ax.legend(fontsize=4.8, loc="upper left", bbox_to_anchor=(0, -0.3), ncol=2)
S.panel(ax, "b", -0.3)

ax = fig.add_subplot(gs[0, 2]); M = [("ichorCNA 1000kb_lowTF", "ichorCNA\nlow-TF", S.color("ichorCNA")), ("ACE ploidy2 100kb", "ACE\n100 kb", S.color("ACE")), ("Fragle", "Fragle", S.color("Fragle"))]
for i, (m, l, col) in enumerate(M):
    v = TF[(TF.method == m) & (TF.level == "full")].TF.values * 100
    ax.scatter(i + rng.uniform(-.18, .18, len(v)), v, s=7, color=col, lw=0, alpha=.8); ax.plot([i - .3, i + .3], [np.median(v)] * 2, color="k", lw=1)
    src += [{"panel": "c", "method": m, "tumour_fraction_pct": x} for x in v]
ax.set_xticks(range(len(M))); ax.set_xticklabels([l for _, l, _ in M]); ax.set_ylabel("Reported tumour fraction (%)\n(true value 0%)"); S.panel(ax, "c", -0.3)
# ---------------- cancer plasma
from matplotlib.colors import ListedColormap
W = pd.read_csv(C.D / "cancer_window_states.tsv", sep="\t", dtype={"chrom": str}); win = pd.read_csv(C.D / "grch37_1Mb_windows.tsv", sep="\t", dtype={"chrom": str})
key = {(c, s_): i for i, (c, s_) in enumerate(zip(win.chrom, win.start))}; N = len(win)
samples = ["IC15", "IC17", "IC20", "IC35", "IC37", "IH03"]; labels = {"IC15": "IC15 SCLC", "IC17": "IC17 HCC", "IC20": "IC20 LUSC", "IC35": "IC35 breast", "IC37": "IC37 CRC", "IH03": "IH03 healthy"}
tools = ["consensus", "QDNAseq", "ichorCNA", "BIC-seq2", "CNVkit", "Spectre"]
rows, ylab = [], []
for b in samples:
    for t in tools:
        v = np.zeros(N); g = W[(W["sample"] == b) & (W.level == "full") & (W.tool == t)]
        for c, s_, x in zip(g.chrom, g.start, g.state): v[key[(c, s_)]] = x
        rows.append(v); ylab.append(f"{labels[b]}  {t}" if t == "consensus" else t)
    rows.append(np.full(N, np.nan)); ylab.append("")
M = np.array(rows[:-1]); ylab = ylab[:-1]
ax = fig.add_subplot(gs[1, :]); cmap = ListedColormap([S.LOSS_COLOR, "#f4f4f4", S.GAIN_COLOR]); cmap.set_bad("white")
ax.imshow(np.ma.masked_invalid(M), aspect="auto", cmap=cmap, vmin=-1, vmax=1, interpolation="nearest")
ch = win.chrom.values; edges = [i for i in range(1, N) if ch[i] != ch[i - 1]]
for e in edges: ax.axvline(e - 0.5, color="#bbbbbb", lw=0.4)
mids = [(a + b_) / 2 for a, b_ in zip([0] + edges, edges + [N])]
ax.set_xticks(mids); ax.set_xticklabels([str(i) for i in range(1, 23)], fontsize=5); ax.set_yticks(range(len(ylab))); ax.set_yticklabels(ylab, fontsize=4.3)
for i, l in enumerate(ylab):
    if "consensus" in l: ax.get_yticklabels()[i].set_fontweight("bold")
ax.tick_params(length=0); [sp.set_visible(False) for sp in ax.spines.values()]
ax.set_xlabel("Chromosome (1-Mb windows, GRCh37)"); S.panel(ax, "d", -0.08, 1.01)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=S.GAIN_COLOR, label="gain"), Patch(color=S.LOSS_COLOR, label="loss")], loc="lower right", bbox_to_anchor=(1, 1.0), ncol=2, fontsize=5.5)
src += [{"panel": "d", "sample": b, "tool": t, "chrom": c, "start": s_, "state": x} for b, t, c, s_, x in zip(W["sample"], W.tool, W.chrom, W.start, W.state) if True]

L = pd.read_csv(C.D / "cancer_loo_and_tcga.tsv", sep="\t"); HB = pd.read_csv(C.D / "healthy_tcga_baseline.tsv", sep="\t")
coh = {"IC15": "LUSC", "IC17": "LIHC", "IC20": "LUSC", "IC35": "BRCA", "IC37": "COADREAD"}
t5 = ["QDNAseq", "ichorCNA", "BIC-seq2", "CNVkit", "Spectre"]
ax = fig.add_subplot(gs[2, 0:2]); mk = dict(zip(coh, ["o", "s", "^", "D", "v"]))
for i, t in enumerate(t5):
    hb = HB[(HB.tool == t)].tcga_r.values
    ax.boxplot([hb], positions=[i], widths=0.45, showfliers=True, patch_artist=True, boxprops=dict(facecolor="#e6e6e6", lw=0.5), medianprops=dict(color="#555555"), whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5), flierprops=dict(ms=1.5, mec="#999999"))
    for b, m_ in mk.items():
        v = L[(L["sample"] == b) & (L.tool == t)].tcga_r.values
        ax.scatter([i + 0.32], v, marker=m_, s=14, color=S.color(t), edgecolor="k", lw=0.3, label=labels[b] if i == 0 else None, zorder=3)
    src += [{"panel": "e", "tool": t, "group": "non-cancer", "tcga_r": x} for x in hb]
    src += [{"panel": "e", "tool": t, "group": b, "tcga_r": x} for b, x in zip(L[L.tool == t]["sample"], L[L.tool == t].tcga_r)]
ax.axhline(0, color="#bbbbbb", lw=0.5); ax.set_xticks(range(len(t5))); ax.set_xticklabels(t5); ax.set_ylabel("Correlation with TCGA\ngain/loss frequency")
ax.legend(fontsize=5, ncol=3, loc="lower right"); ax.set_ylim(-0.75, 0.8); ax.set_title("Cancer plasma (symbols) vs 25 non-cancer samples (boxes)", loc="left", fontsize=6.5); S.panel(ax, "e", -0.08)

R_ = pd.read_csv(C.D / "cancer_depth_robustness.tsv", sep="\t")
ax = fig.add_subplot(gs[2, 2])
for i, t in enumerate(t5):
    for k, (lvl, dx) in enumerate([("d1", -0.15), ("d0p1", 0.15)]):
        v = R_[(R_.tool == t) & (R_.level == lvl)].F1.values
        ax.scatter(np.full(len(v), i + dx) + rng.uniform(-.05, .05, len(v)), v, s=8, color=S.color(t), alpha=1 if lvl == "d1" else 0.5, lw=0)
        ax.plot([i + dx - .1, i + dx + .1], [np.median(v)] * 2, color="k", lw=0.9)
        src += [{"panel": "f", "tool": t, "level": lvl, "F1_vs_native": x} for x in v]
ax.set_xticks(range(len(t5))); ax.set_xticklabels(t5, rotation=35, ha="right"); ax.set_ylim(0, 1.05)
ax.set_ylabel("Agreement with own native-\ndepth profile (F1)"); ax.set_title("1× (left) and 0.1× (right)", loc="left", fontsize=6.5); S.panel(ax, "f", -0.35)
S.save(fig, "Figure4_real_plasma"); pd.DataFrame(src).to_csv("output/Figure4_source_data.tsv", sep="\t", index=False)
