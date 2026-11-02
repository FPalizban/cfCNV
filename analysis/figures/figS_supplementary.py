"""Supplementary figures that are drawn directly from the data files in data/.
S1   measured vs nominal depth (Simulation 1)
S9   profile correlation in the 15 pure-normal HG008 mixtures, with detection thresholds
S12  annotation of false calls in non-cancer plasma (enrichment vs size-matched random placement)
S13  Fragle ctDNA fraction vs proportion of fragments < 150 bp
S14  recurrent false-positive windows per caller (1-Mb windows altered in >= 13 of 25 individuals)
Run:  python figS_supplementary.py"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.stats import spearmanr
import style as S, common as C
S.apply()

# S1
md = pd.read_csv(C.D / "sim1_measured_depth.tsv", sep="\t")
md["nominal"] = md.sample_id.str.extract(r"_([\d.]+)x$")[0].astype(float); md["bg"] = md.sample_id.str[:4]
fig, ax = plt.subplots(figsize=(85 * S.MM, 70 * S.MM))
for bg, g in md.groupby("bg"): ax.scatter(g.nominal, g.measured_depth, s=10, color=S.BACKGROUND_COLORS.get(bg, "#333"), label=bg, lw=0)
ax.plot([0.08, 25], [0.08, 25], color="#999999", lw=0.8); ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Nominal depth (×)"); ax.set_ylabel("Measured autosomal depth (×)"); ax.legend(); S.save(fig, "FigureS1_measured_depth"); plt.close(fig)

# S9
pr, tf, mix, prim = C.sim2(); P = C.primary(pr, prim); L = C.primary(C.lod_strict(pr), prim)
fig, ax = plt.subplots(figsize=(120 * S.MM, 65 * S.MM)); tools = [t for t in S.TOOL_ORDER if t in prim.index]; rng = np.random.default_rng(2)
for i, t in enumerate(tools):
    v = P[(P.tool == t) & (P.nominal_fraction == 0)].pearson_r.values
    ax.scatter(i + rng.uniform(-.2, .2, len(v)), v, s=8, color=S.color(t), lw=0)
    th = L[L.tool == t].threshold.iloc[0]; ax.plot([i - .3, i + .3], [th, th], color="k", lw=1.2)
ax.axhline(0, color="#bbbbbb", lw=0.6); ax.set_xticks(range(len(tools))); ax.set_xticklabels(tools, rotation=35, ha="right")
ax.set_ylabel("Profile r in pure-normal mixtures"); S.save(fig, "FigureS9_null_profile_correlation"); plt.close(fig)

# S12
B, TF, INS, CALLS, REC, ANN = C.healthy(); A = ANN[ANN.level == "full"]
cats = [("frac_germline_CNV", "Germline CNV"), ("frac_segdup", "Segmental dup."), ("frac_blacklist", "Blacklist"), ("frac_pericentromeric", "Pericentromeric")]
tools = ["QDNAseq", "Spectre", "BIC-seq2", "CNVkit", "ichorCNA"]
fig, axs = plt.subplots(1, 3, figsize=(S.FULL_WIDTH, 62 * S.MM), gridspec_kw=dict(wspace=0.5))
ax = axs[0]
for j, t in enumerate(tools):
    for i, (m, _) in enumerate(cats):
        r = A[(A.tool == t) & (A.metric == m)]
        if r.empty: continue
        ax.scatter(i + (j - 2) * 0.13, max(r.enrichment.iloc[0], 0.05), color=S.color(t), s=20 if r.p_perm.iloc[0] <= 0.05 else 8, edgecolor="k" if r.p_perm.iloc[0] <= 0.05 else "none", lw=0.4, label=t if i == 0 else None)
ax.axhline(1, color="#999999", ls="--", lw=0.7); ax.set_yscale("log"); ax.set_xticks(range(len(cats))); ax.set_xticklabels([l for _, l in cats], rotation=30, ha="right")
ax.set_ylabel("Enrichment vs random placement"); ax.legend(fontsize=5); S.panel(ax, "a", -0.3)
ax = axs[1]; o = [A[(A.tool == t) & (A.metric == "frac_any")].observed.iloc[0] * 100 for t in tools]; e = [A[(A.tool == t) & (A.metric == "frac_any")].expected.iloc[0] * 100 for t in tools]
ax.bar(range(len(tools)), o, color=[S.color(t) for t in tools]); ax.scatter(range(len(tools)), e, marker="_", s=150, color="k", label="random placement"); ax.set_ylim(0, 80)
ax.set_xticks(range(len(tools))); ax.set_xticklabels(tools, rotation=30, ha="right"); ax.set_ylabel("Called bp in known non-tumour\nor artefact-prone regions (%)"); ax.legend(fontsize=5, loc="upper left"); S.panel(ax, "b", -0.35)
ax = axs[2]
for j, t in enumerate(tools):
    for i, c in enumerate(["LIHC", "LUSC", "COADREAD", "BRCA"]):
        r = A[(A.tool == t) & (A.metric == f"tcga_{c}")]
        ax.scatter(i + (j - 2) * 0.13, r.enrichment.iloc[0], color=S.color(t), s=20 if r.p_perm.iloc[0] <= 0.05 else 8, edgecolor="k" if r.p_perm.iloc[0] <= 0.05 else "none", lw=0.4)
ax.axhline(1, color="#999999", ls="--", lw=0.7); ax.set_xticks(range(4)); ax.set_xticklabels(["LIHC", "LUSC", "COAD/READ", "BRCA"], rotation=30, ha="right")
ax.set_ylabel("TCGA same-direction frequency\n(observed / random)"); S.panel(ax, "c", -0.35)
S.save(fig, "FigureS12_noncancer_call_annotation"); plt.close(fig)

# S13
f = TF[(TF.method == "Fragle") & (TF.level == "full")].merge(INS, on="base")
rho, p = spearmanr(f.frac_lt150, f.TF)
fig, ax = plt.subplots(figsize=(85 * S.MM, 70 * S.MM)); ax.scatter(f.frac_lt150 * 100, f.TF * 100, s=12, color=S.color("Fragle"), lw=0)
ax.set_xlabel("Fragments < 150 bp (%)"); ax.set_ylabel("Fragle ctDNA fraction (%)"); ax.text(0.03, 0.95, f"Spearman ρ = {rho:.2f}, p = {p:.1g}", transform=ax.transAxes, va="top")
S.save(fig, "FigureS13_fragle_fragment_length"); plt.close(fig)

# S14
fig, ax = plt.subplots(figsize=(85 * S.MM, 60 * S.MM)); cf = C.HEALTHY_PRIMARY
n = [REC[(REC.tool == t) & (REC.config == c)].groupby(["chrom", "start"]).n_samples.max().ge(13).sum() for t, c in cf]
ax.barh(range(len(cf)), n, color=[S.color(t) for t, _ in cf]); ax.set_yticks(range(len(cf))); ax.set_yticklabels([f"{t} {c}" for t, c in cf]); ax.invert_yaxis()
for i, v in enumerate(n): ax.text(v + 3, i, str(v), va="center")
ax.set_xlabel("1-Mb windows altered in ≥13 of 25 individuals"); S.save(fig, "FigureS14_recurrent_windows"); plt.close(fig)
