#!/usr/bin/env python3
"""
cfCNV Simulation 1 - manuscript figures (vector PDF + 600-dpi PNG + source data).

  python 03_figures.py --tables tables/ --stats stats/ --out figures/

Figure 1  Overall performance: precision-recall plane (all tools x bin sizes, F1
          iso-lines) + per-sample F1 and bin-level MCC at each tool's primary bin.
Figure 2  Sequencing depth within each cfDNA background: F1 (spike-in samples)
          and false calls per Gb (unmodified negative controls).
Figure 3  CNV size: recall by simulated size at the primary bin + bin x size
          recall heatmaps for the bin-scanned tools.
Figure 4  Bin size: F1, precision and recall vs bin size (matched depth).
Figure 5  Negative controls: false calls per Gb and % genome falsely altered.
Figure 6  Runtime per run.
Figure S1 Measured vs nominal depth. S2 Reciprocal-overlap sensitivity.
Figure S3 Configuration sensitivity (ichorCNA, CNVkit, Spectre). S4 Background noise anatomy.
Figure S5 LYCEUM diagnostics.

Colours: one fixed, colour-blind-validated hue + marker per tool. Text 6.5-8 pt,
widths 180 mm (double column) / 88 mm (single column).
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mt
import numpy as np
import pandas as pd
from matplotlib import font_manager

MM = 1 / 25.4
ORDER = ["QDNAseq", "ichorCNA", "BICseq2", "CNV-MEANN", "dudeML", "CNVkit", "Spectre", "LYCEUM"]
LABEL = {"BICseq2": "BIC-seq2", "ichorCNA_default": "ichorCNA (tumour defaults)",
         "CNVkit_flat": "CNVkit (flat reference)", "Spectre_1kb": "Spectre (1-kb windows)"}
COLOR = {"QDNAseq": "#0072B2", "CNVkit": "#E69F00", "ichorCNA": "#009E73", "BICseq2": "#D55E00",
         "Spectre": "#CC79A7", "dudeML": "#56B4E9", "CNV-MEANN": "#882255", "LYCEUM": "#6A3D9A"}
MARK = {"QDNAseq": "o", "CNVkit": "s", "ichorCNA": "^", "BICseq2": "D", "Spectre": "v",
        "dudeML": "P", "CNV-MEANN": "X", "LYCEUM": "*"}
INK, MUTED, GRID = "#1a1a1a", "#5c5c5c", "#e3e3e3"
BGS = ["IH01", "IH02", "IH03"]
SIZES = [50e3, 100e3, 250e3, 500e3, 1e6, 5e6, 10e6, 25e6, 50e6]
SIZE_LAB = ["50 kb", "100 kb", "250 kb", "500 kb", "1 Mb", "5 Mb", "10 Mb", "25 Mb", "50 Mb"]


def style():
    have = {f.name for f in font_manager.fontManager.ttflist}
    fam = next((f for f in ("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans") if f in have), "sans-serif")
    plt.rcParams.update({
        "font.family": fam, "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7.5,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
        "text.color": INK, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
        "axes.axisbelow": True, "lines.linewidth": 1.3, "lines.markersize": 4, "legend.frameon": False,
        "pdf.fonttype": 42, "savefig.dpi": 600, "savefig.bbox": "tight", "savefig.pad_inches": 0.03})


lab = lambda t: LABEL.get(t, t)


def save(fig, out, name, data=None):
    fig.savefig(out / f"{name}.pdf"); fig.savefig(out / f"{name}.png"); plt.close(fig)
    if data is not None:
        data.to_csv(out / f"{name}_source_data.tsv", sep="\t", index=False)
    print("  ", name)


def plab(ax, s, x=-0.16):
    ax.text(x, 1.04, s, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")


def legend(fig, tools, y=-0.01, ncol=None):
    h = [plt.Line2D([], [], color=COLOR[t], marker=MARK[t], lw=1.3, ms=4) for t in tools]
    fig.legend(h, [lab(t) for t in tools], loc="lower center", ncol=ncol or len(tools),
               bbox_to_anchor=(0.5, y), handlelength=2.0, columnspacing=1.3)


def depth_axis(ax, depths):
    ax.set_xscale("log"); ax.set_xticks(depths)
    ax.set_xticklabels([f"{d:g}" for d in depths]); ax.minorticks_off()


def strip(ax, data, tools, col, log=False, jitter=0.18, seed=1, linthresh=1):
    r = np.random.default_rng(seed)
    for i, t in enumerate(tools):
        v = data[data.tool == t][col].dropna().values
        if not len(v):
            continue
        ax.scatter(i + r.uniform(-jitter, jitter, len(v)), v, s=7, color=COLOR[t], alpha=0.55,
                   linewidth=0, zorder=2)
        q1, md, q3 = np.percentile(v, [25, 50, 75])
        ax.plot([i - 0.3, i + 0.3], [md, md], color=INK, lw=1.4, zorder=3, solid_capstyle="butt")
        ax.plot([i, i], [q1, q3], color=INK, lw=0.8, zorder=3)
    ax.set_xticks(range(len(tools))); ax.set_xticklabels([lab(t) for t in tools], rotation=30, ha="right")
    ax.grid(axis="x", visible=False)
    if log:
        ax.set_yscale("symlog", linthresh=linthresh); ax.set_ylim(bottom=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", required=True); ap.add_argument("--stats", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--inputs", default="inputs")
    a = ap.parse_args()
    style()
    T, S, O = Path(a.tables), Path(a.stats), Path(a.out); O.mkdir(parents=True, exist_ok=True)
    pc = pd.read_csv(T / "per_config.tsv", sep="\t")
    td = pd.read_csv(T / "truth_detection.tsv", sep="\t")
    rt = pd.read_csv(T / "runtime.tsv", sep="\t")
    md = pd.read_csv(T / "measured_depth.tsv", sep="\t")
    prim = pd.read_csv(S / "primary_configs.tsv", sep="\t").set_index("tool").bin_bp.to_dict()
    isprim = lambda d: d[[prim.get(t) == b for t, b in zip(d.tool, d.bin_bp)]]
    pos = pc[(pc["class"] == "positive") & pc.tool.isin(ORDER)]
    neg = pc[(pc["class"] == "negative") & pc.tool.isin(ORDER)]
    pm = pos[pos.matched_depth]
    tools = [t for t in ORDER if t in prim]
    active = list(tools)
    depths = sorted(pos.depth.unique())

    # ---------------- Figure 1 ----------------
    fig = plt.figure(figsize=(180 * MM, 66 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1, 1], wspace=0.42)
    ax = fig.add_subplot(gs[0])
    rr = np.linspace(0.01, 1, 200)
    for f in (0.2, 0.4, 0.6, 0.8):
        pr = f * rr / (2 * rr - f); ok = (pr > 0) & (pr <= 1)
        ax.plot(rr[ok], pr[ok], color=GRID, lw=0.8, zorder=1)
        ax.text(0.985, f / (2 - f) + 0.012, f"F1 {f:g}", fontsize=5.5, color=MUTED, ha="right")
    sm = pm.groupby(["tool", "bin_bp"])[["precision", "recall"]].median().reset_index()
    for t in active:
        g = sm[sm.tool == t].sort_values("bin_bp")
        ax.plot(g.recall, g.precision, color=COLOR[t], lw=0.9, alpha=0.7, zorder=2)
        for _, r in g.iterrows():
            big = prim.get(t) == r.bin_bp
            ax.scatter(r.recall, r.precision, s=34 if big else 14, marker=MARK[t], color=COLOR[t],
                       edgecolor=INK if big else "white", linewidth=0.7 if big else 0.4, zorder=3)
        e = g.iloc[[0, -1]] if len(g) > 1 else g
        off = {("BICseq2", 0): (4, -9), ("BICseq2", 1): (2, -10), ("dudeML", 0): (6, 1),
               ("CNVkit", 0): (5, -3), ("QDNAseq", 1): (-4, -9), ("ichorCNA", 0): (5, 3),
               ("Spectre", 0): (5, 2)}
        if t == "LYCEUM":
            continue
        if t == "ichorCNA" and len(e) == 2 and np.allclose(e.recall, e.recall.iloc[0], atol=0.02) and np.allclose(e.precision, e.precision.iloc[0], atol=0.02):
            r = e.iloc[0]
            ax.annotate("500 & 1000", (r.recall, r.precision), xytext=(6, -2), textcoords="offset points",
                        fontsize=5.5, color=COLOR[t], ha="left", va="top")
            continue
        for j, (_, r) in enumerate(e.iterrows()):
            ax.annotate(f"{r.bin_bp/1000:g}", (r.recall, r.precision), xytext=off.get((t, j), (4, 3)),
                        textcoords="offset points", fontsize=5.5, color=COLOR[t],
                        ha="right" if off.get((t, j), (1, 0))[0] < 0 else "left")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
    ax.set_xlabel("Event recall"); ax.set_ylabel("Event precision")
    ax.set_title("Median per tool × bin size (kb labelled)", loc="left"); plab(ax, "a", -0.18)
    ax2 = fig.add_subplot(gs[1]); p1 = isprim(pm)
    strip(ax2, p1, active, "F1"); ax2.set_ylim(-0.03, 1.03); ax2.set_ylabel("Event F1 (RO ≥ 0.5)")
    ax2.set_title("Primary bin, per sample", loc="left"); plab(ax2, "b", -0.3)
    ax3 = fig.add_subplot(gs[2])
    strip(ax3, p1, [t for t in tools if p1[p1.tool == t].bin_MCC.notna().any()], "bin_MCC")
    ax3.set_ylim(-0.03, 1.03); ax3.set_ylabel("Bin-level MCC (10-kb grid)")
    ax3.set_title("Primary bin, per sample", loc="left"); plab(ax3, "c", -0.3)
    legend(fig, active, y=-0.14)
    save(fig, O, "Fig1_overall_performance", pd.concat([sm.assign(panel="a"), p1.assign(panel="b,c")]))

    # ---------------- Figure 2 ----------------
    p2, n2 = isprim(pos), isprim(neg)
    fig, axs = plt.subplots(2, 3, figsize=(180 * MM, 108 * MM), sharex="col", sharey="row")
    for c, bg in enumerate(BGS):
        dd = sorted(pos[pos.background == bg].depth.unique())
        for r, (d, col, yl) in enumerate([(p2, "F1", "Event F1 (spike-in samples)"),
                                          (n2, "false_calls_per_Gb", "False calls per Gb\n(negative controls)")]):
            ax = axs[r][c]
            ax.axvspan(0.08, 1.15, color="#f2f2f2", zorder=0, lw=0)
            for t in active:
                g = d[(d.tool == t) & (d.background == bg)]
                s = g.groupby("depth")[col].agg(["median", lambda x: x.quantile(.25), lambda x: x.quantile(.75)])
                s.columns = ["med", "q1", "q3"]
                if not len(s):
                    continue
                ax.plot(s.index, s.med, color=COLOR[t], marker=MARK[t])
                if col == "F1":
                    ax.fill_between(s.index, s.q1, s.q3, color=COLOR[t], alpha=0.12, lw=0)
            depth_axis(ax, dd)
            if col != "F1":
                ax.set_yscale("symlog", linthresh=1); ax.set_ylim(bottom=0)
                ax.yaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: f"{v:g}"))
                ax.set_xlabel("Sequencing depth (×)")
            else:
                ax.set_ylim(-0.03, 1.03); ax.set_title(f"{bg}", loc="left")
            if c == 0:
                ax.set_ylabel(yl)
            plab(ax, "abcdef"[r * 3 + c], -0.12 if c else -0.22)
    legend(fig, active, y=-0.05); fig.subplots_adjust(hspace=0.18, wspace=0.08)
    save(fig, O, "Fig2_depth_by_background", pd.concat([p2, n2]))

    # ---------------- Figure 3 ----------------
    tp = isprim(td); tp = tp[tp.depth <= 1]
    rs = tp.groupby(["tool", "size_bp"]).detected.mean().reset_index()
    scan = [t for t in ("QDNAseq", "BICseq2", "dudeML") if t in tools]
    fig = plt.figure(figsize=(180 * MM, 70 * MM))
    gs = fig.add_gridspec(1, 1 + len(scan), width_ratios=[1.5] + [1] * len(scan), wspace=0.35)
    ax = fig.add_subplot(gs[0])
    for t in active:
        g = rs[rs.tool == t].set_index("size_bp").reindex(SIZES)
        ax.plot(range(len(SIZES)), g.detected, color=COLOR[t], marker=MARK[t])
    ax.set_xticks(range(len(SIZES))); ax.set_xticklabels(SIZE_LAB, rotation=45, ha="right")
    ax.set_ylim(-0.03, 1.03); ax.set_ylabel("Recall (RO ≥ 0.5)"); ax.set_xlabel("Simulated CNV size")
    ax.set_title("Primary bin, 0.1–1×", loc="left"); plab(ax, "a")
    rb = td[td.depth <= 1].groupby(["tool", "bin_bp", "size_bp"]).detected.mean().reset_index()
    for k, t in enumerate(scan):
        ax = fig.add_subplot(gs[k + 1])
        h = rb[rb.tool == t].pivot(index="bin_bp", columns="size_bp", values="detected").reindex(columns=SIZES)
        im = ax.imshow(h.values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(len(SIZES))); ax.set_xticklabels(SIZE_LAB, rotation=60, ha="right", fontsize=5.5)
        ax.set_yticks(range(len(h.index))); ax.set_yticklabels([f"{b/1000:g}" for b in h.index])
        ax.grid(False); [s.set_visible(False) for s in ax.spines.values()]
        if k == 0:
            ax.set_ylabel("Bin size (kb)")
        ax.set_title(lab(t), loc="left", color=INK); plab(ax, "bcd"[k], -0.22)
    cb = fig.colorbar(im, ax=fig.axes[1:], shrink=0.75, pad=0.02); cb.set_label("Recall"); cb.outline.set_visible(False)
    legend(fig, active, y=-0.2, ncol=5)
    save(fig, O, "Fig3_cnv_size", pd.concat([rs.assign(panel="a"), rb.assign(panel="b-d")]))

    # ---------------- Figure 4 ----------------
    bt = [t for t in ("QDNAseq", "BICseq2", "CNV-MEANN", "dudeML", "CNVkit", "Spectre") if t in tools and pm[pm.tool == t].bin_bp.nunique() > 2]
    fig, axs = plt.subplots(1, 3, figsize=(180 * MM, 60 * MM), sharex=True)
    for k, (col, yl) in enumerate([("F1", "Event F1"), ("precision", "Event precision"), ("recall", "Event recall")]):
        ax = axs[k]
        for t in bt:
            g = pm[pm.tool == t].groupby("bin_bp")[col]
            s = pd.DataFrame({"med": g.median(), "q1": g.quantile(.25), "q3": g.quantile(.75)})
            ax.plot(s.index / 1000, s.med, color=COLOR[t], marker=MARK[t])
            ax.fill_between(s.index / 1000, s.q1, s.q3, color=COLOR[t], alpha=0.08, lw=0)
        ax.set_xscale("log"); ks = [10, 30, 50, 100, 500, 1000]
        ax.set_xticks(ks); ax.set_xticklabels([f"{k:g}" for k in ks]); ax.minorticks_off()
        ax.set_ylim(-0.03, 1.03); ax.set_xlabel("Bin size (kb)"); ax.set_ylabel(yl + " (RO ≥ 0.5)")
        plab(ax, "abc"[k])
    legend(fig, bt, y=-0.12); fig.subplots_adjust(wspace=0.35)
    save(fig, O, "Fig4_bin_size", pm[pm.tool.isin(bt)])

    # ---------------- Figure 5 ----------------
    nm = isprim(neg[neg.matched_depth]).copy(); nm["false_altered_pct"] = 100 * nm.false_altered_fraction
    fig, axs = plt.subplots(1, 2, figsize=(120 * MM, 60 * MM))
    strip(axs[0], nm, active, "false_calls_per_Gb", log=True); axs[0].set_ylabel("False calls per Gb")
    axs[0].yaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: f"{v:g}")); plab(axs[0], "a", -0.28)
    strip(axs[1], nm, active, "false_altered_pct", log=True, linthresh=0.01); axs[1].set_ylabel("Genome falsely altered (%)")
    axs[1].yaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: f"{v:g}")); plab(axs[1], "b", -0.28)
    for ax in axs:
        ax.set_title("Negative controls, 0.1–1×, primary bin", loc="left", fontsize=6.5)
    fig.subplots_adjust(wspace=0.45)
    save(fig, O, "Fig5_negative_controls", nm)

    # ---------------- Figure 6 ----------------
    r6 = rt[rt.state.str.startswith("COMPLETED") & rt.tool.isin(ORDER)].copy(); r6["min"] = r6.wall_s / 60
    rtools = [t for t in ORDER if t in set(r6.tool)]
    fig, ax = plt.subplots(figsize=(88 * MM, 60 * MM))
    strip(ax, r6, rtools, "min", jitter=0.25)
    ax.set_yscale("log"); ax.yaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_ylabel("Wall-clock time per run (min)")
    for i, t in enumerate(rtools):
        nf = int((rt[rt.tool == t].state == "FAILED").sum())
        if nf:
            ax.text(i, ax.get_ylim()[1], f"{nf} failed", ha="center", va="bottom", fontsize=5.5, color=MUTED)
    if "Spectre" in set(rt.tool) and not len(r6[r6.tool == "Spectre"]):
        ax.text(0.99, 0.02, f"Spectre: {int((rt.tool == 'Spectre').sum())}/"
                f"{int((rt.tool == 'Spectre').sum())} runs failed (not shown)",
                transform=ax.transAxes, ha="right", fontsize=5.5, color=MUTED)
    save(fig, O, "Fig6_runtime", r6)

    # ---------------- Figure S1 ----------------
    md["nominal"] = md.sample_id.str.extract(r"_([\d.]+)x$")[0].astype(float)
    md["background"] = md.sample_id.str[:4]
    fig, ax = plt.subplots(figsize=(88 * MM, 66 * MM))
    for bg, m in zip(BGS, ["o", "s", "^"]):
        g = md[md.background == bg]
        ax.scatter(g.nominal, g.measured_depth, s=12, marker=m, facecolor="none", edgecolor=INK, lw=0.6, label=bg)
    ax.plot([0.08, 25], [0.08, 25], color=MUTED, ls="--", lw=0.7)
    ax.set_xscale("log"); ax.set_yscale("log")
    for f in (ax.xaxis, ax.yaxis):
        f.set_major_formatter(mt.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlabel("Nominal depth (×)"); ax.set_ylabel("Measured autosomal depth (×, mosdepth)")
    ax.legend(loc="upper left")
    save(fig, O, "FigS1_measured_depth", md)

    # ---------------- Figure S2 ----------------
    s2 = isprim(pm).groupby("tool")[["recall_10", "recall_30", "recall"]].median().reindex(active)
    fig, ax = plt.subplots(figsize=(88 * MM, 60 * MM))
    for t in active:
        ax.plot([10, 30, 50], s2.loc[t].values, color=COLOR[t], marker=MARK[t], label=lab(t))
    ax.set_xticks([10, 30, 50]); ax.set_xlabel("Reciprocal-overlap threshold (%)")
    ax.set_ylabel("Median event recall"); ax.set_ylim(0, 1); ax.legend(ncol=2, loc="lower left")
    save(fig, O, "FigS2_overlap_sensitivity", s2.reset_index())

    # ---------------- Figure S3: configuration sensitivity ----------------
    pairs = [("ichorCNA", "ichorCNA_default", "ichorCNA, 1 Mb", "pure-sample", "tumour defaults"),
             ("CNVkit", "CNVkit_flat", "CNVkit, 100 kb", "pooled normals", "flat reference"),
             ("Spectre", "Spectre_1kb", "Spectre", "benchmark bin", "1-kb windows")]
    fig, axs = plt.subplots(2, 3, figsize=(180 * MM, 105 * MM))
    src = []
    for k, (m, alt, title, lm, la) in enumerate(pairs):
        for t, ls, mk, lb in ((m, "-", MARK[m], lm), (alt, "--", "o", la)):
            b = prim.get(t, prim.get(m)) if (pc[(pc.tool == t)].bin_bp == prim.get(t, prim.get(m))).any() else pc[pc.tool == t].bin_bp.iloc[0]
            p = pc[(pc.tool == t) & (pc.bin_bp == b) & (pc["class"] == "positive")]
            n = pc[(pc.tool == t) & (pc.bin_bp == b) & (pc["class"] == "negative")]
            g = p.groupby("depth").F1.median(); h = n.groupby("depth").false_altered_fraction.median() * 100
            axs[0, k].plot(g.index, g.values, color=COLOR[m], ls=ls, marker=mk, mfc="white" if t == alt else COLOR[m], label=lb)
            axs[1, k].plot(h.index, h.values, color=COLOR[m], ls=ls, marker=mk, mfc="white" if t == alt else COLOR[m], label=lb)
            src.append(pd.DataFrame({"tool": t, "bin_bp": b, "depth": g.index, "median_F1": g.values,
                                     "median_false_altered_pct": h.reindex(g.index).values}))
        for r in (0, 1):
            depth_axis(axs[r, k], sorted(pc.depth.dropna().unique()))
        axs[0, k].set_ylim(0, 1); axs[0, k].set_title(title, loc="left"); axs[0, k].legend(loc="upper left")
        axs[1, k].set_xlabel("Sequencing depth (×)"); axs[1, k].set_yscale("symlog", linthresh=0.1); axs[1, k].set_ylim(bottom=0)
        plab(axs[0, k], "abc"[k], -0.22)
    axs[0, 0].set_ylabel("Median F1 (spike-ins)"); axs[1, 0].set_ylabel("Genome falsely altered (%)\n(negative controls)")
    fig.subplots_adjust(wspace=0.32, hspace=0.35)
    save(fig, O, "FigS3_configuration_sensitivity", pd.concat(src))

    # ---------------- Figure S4: background noise anatomy ----------------
    nz = pd.read_csv(Path(a.inputs) / "noise_by_sample.tsv", sep="\t")
    nz = nz[nz.bin_bp == nz.bin_bp.min()] if (nz.bin_bp == 50000).sum() == 0 else nz[nz.bin_bp == 50000]
    bg_col = {"IH01": "#882255", "IH02": "#117733", "IH03": "#332288"}
    fig, axs = plt.subplots(1, 3, figsize=(180 * MM, 58 * MM))
    for bg, d in nz.groupby("background"):
        d = d.groupby("depth")[["observed_sd_log2", "poisson_sd_log2", "systematic_sd_log2", "lag1_autocorr"]].median()
        c = bg_col.get(bg, INK)
        axs[0].plot(d.index, d.observed_sd_log2, color=c, marker="o", label=bg)
        axs[1].plot(d.index, d.systematic_sd_log2, color=c, marker="o", label=f"{bg} systematic")
        axs[1].plot(d.index, d.poisson_sd_log2, color=c, ls=":", lw=1.0)
        axs[2].plot(d.index, d.lag1_autocorr, color=c, marker="o", label=bg)
    for ax in axs:
        depth_axis(ax, sorted(nz.depth.unique())); ax.set_xlabel("Sequencing depth (×)")
    axs[0].set_ylabel("SD of bin log2 ratio"); axs[0].set_title("Observed coverage noise", loc="left"); axs[0].legend()
    axs[1].set_ylabel("SD of bin log2 ratio"); axs[1].set_title("Systematic (—) vs counting (···)", loc="left")
    axs[2].set_ylabel("Lag-1 autocorrelation"); axs[2].set_title("Wave-like structure", loc="left"); axs[2].set_ylim(0, 1)
    for k, ax in enumerate(axs):
        plab(ax, "abc"[k], -0.22)
    fig.subplots_adjust(wspace=0.45)
    save(fig, O, "FigS5_background_noise", nz)

    # ---------------- Figure S5: LYCEUM diagnostics ----------------
    ld = pd.read_csv(Path(a.inputs) / "LYCEUM_diagnostics.tsv", sep="\t")
    ld = ld[ld["class"] == "positive"]
    fig, axs = plt.subplots(1, 2, figsize=(150 * MM, 58 * MM))
    for col, mk, lb in (("AUC_del", "v", "Deletions"), ("AUC_dup", "^", "Duplications")):
        g = ld.groupby("depth_x")[col].median()
        axs[0].plot(g.index, g.values, color=COLOR["LYCEUM"], marker=mk, label=lb)
    axs[0].axhline(0.5, color=INK, lw=0.7, ls=":")
    depth_axis(axs[0], sorted(ld.depth_x.unique())); axs[0].set_ylim(0.3, 1); axs[0].set_xlabel("Sequencing depth (×)")
    axs[0].set_ylabel("Exon-level AUC of LYCEUM probability"); axs[0].legend(loc="upper left")
    axs[0].set_title("Model discrimination", loc="left"); plab(axs[0], "a", -0.2)
    for col, mk, lb in (("depth_ratio_loss", "v", "Truth deletions"), ("depth_ratio_gain", "^", "Truth duplications")):
        g = ld[ld.median_rd_neutral >= 0.25].groupby("depth_x")[col].median().dropna()
        axs[1].plot(g.index, g.values, color=INK, marker=mk, label=lb)
    axs[1].axhline(0.5, color="#999999", lw=0.6, ls=":"); axs[1].axhline(1.5, color="#999999", lw=0.6, ls=":")
    depth_axis(axs[1], sorted(ld.depth_x.unique())); axs[1].set_xlabel("Sequencing depth (×)")
    axs[1].set_ylabel("Exon read depth / neutral exons"); axs[1].legend(loc="center right"); axs[1].set_ylim(0, 2.2)
    axs[1].set_title("Signal present in the input", loc="left"); plab(axs[1], "b", -0.2)
    fig.subplots_adjust(wspace=0.35)
    save(fig, O, "FigS6_LYCEUM_diagnostics", ld)

if __name__ == "__main__":
    main()
