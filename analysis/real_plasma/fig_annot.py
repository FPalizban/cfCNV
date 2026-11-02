import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, '/home/claude/sim2_final'); from s2_analysis import style, COLOR, MARK, save, plab, MM, INK, MUTED
from pathlib import Path
style(); O = Path(sys.argv[1]); FG = Path(sys.argv[2]); FG.mkdir(exist_ok=True)
T = pd.read_csv(O / 'tool_annotation_vs_random.tsv', sep='\t'); T = T[T.level == 'full']
A = pd.read_csv(O / 'calls_annotated.tsv', sep='\t'); A = A[A.level == 'full']
COL = dict(COLOR); COL['BIC-seq2'] = COLOR['BICseq2']
tools = ['QDNAseq', 'Spectre', 'BIC-seq2', 'CNVkit', 'ichorCNA']
fig, axs = plt.subplots(1, 3, figsize=(180 * MM, 62 * MM), gridspec_kw=dict(width_ratios=[1.25, 1, 1.1], wspace=0.55))
# a enrichment by region class
ax = axs[0]; cats = [('frac_germline_CNV', 'Germline\nCNV'), ('frac_segdup', 'Seg.\ndup.'), ('frac_blacklist', 'Black-\nlist'), ('frac_pericentromeric', 'Peri-\ncentrom.')]
for j, t in enumerate(tools):
    for i, (m, l) in enumerate(cats):
        r = T[(T.tool == t) & (T.metric == m)]
        if r.empty: continue
        v = max(r.enrichment.iloc[0], 0.05); sig = r.p_perm.iloc[0] <= 0.05
        ax.scatter(i + (j - 2) * 0.13, v, color=COL[t], marker=MARK.get(t.replace('BIC-seq2', 'BICseq2'), 'o'), s=22 if sig else 12,
                   edgecolor=INK if sig else 'none', lw=0.5, label=t if i == 0 else None, zorder=3)
ax.axhline(1, color=MUTED, lw=0.8, ls='--'); ax.set_yscale('log'); ax.set_ylim(0.04, 80)
ax.set_xticks(range(len(cats))); ax.set_xticklabels([l for _, l in cats], fontsize=6); ax.grid(axis='x', visible=False)
ax.set_ylabel('Enrichment of called bp\n(observed / size-matched random)'); ax.set_title('Where healthy calls fall', loc='left'); plab(ax, 'a', -0.28)
ax.legend(fontsize=5.5, loc='upper center', ncol=3, handletextpad=0.1, columnspacing=0.5, frameon=False)
# b % called bp in known problem regions vs random
ax = axs[1]; x = np.arange(len(tools))
o = [T[(T.tool == t) & (T.metric == 'frac_any')].observed.iloc[0] * 100 for t in tools]
e = [T[(T.tool == t) & (T.metric == 'frac_any')].expected.iloc[0] * 100 for t in tools]
ax.bar(x, o, color=[COL[t] for t in tools], width=0.62); ax.scatter(x, e, marker='_', s=180, color=INK, zorder=3, label='random placement')
for i, v in enumerate(o): ax.text(i, v + 1.5, f'{v:.0f}', ha='center', fontsize=6)
ax.set_xticks(x); ax.set_xticklabels(tools, rotation=35, ha='right', fontsize=6); ax.set_ylim(0, 75); ax.grid(axis='x', visible=False)
ax.set_ylabel('Called bp in known non-tumour\nor artefact-prone regions (%)'); ax.legend(fontsize=5.5, loc='upper right'); ax.set_title('Explained calls', loc='left'); plab(ax, 'b', -0.38)
# c TCGA concordance enrichment (cancer-likeness) ~ 1
ax = axs[2]; coh = ['LIHC', 'LUSC', 'COADREAD', 'BRCA']
for j, t in enumerate(tools):
    for i, c in enumerate(coh):
        r = T[(T.tool == t) & (T.metric == f'tcga_{c}')]
        sig = r.p_perm.iloc[0] <= 0.05
        ax.scatter(i + (j - 2) * 0.13, r.enrichment.iloc[0], color=COL[t], marker=MARK.get(t.replace('BIC-seq2', 'BICseq2'), 'o'),
                   s=22 if sig else 12, edgecolor=INK if sig else 'none', lw=0.5, zorder=3)
ax.axhline(1, color=MUTED, lw=0.8, ls='--'); ax.set_ylim(0.6, 1.5)
ax.set_xticks(range(4)); ax.set_xticklabels(['LIHC', 'LUSC', 'COAD/\nREAD', 'BRCA'], fontsize=6); ax.grid(axis='x', visible=False)
ax.set_ylabel('TCGA same-direction frequency\n(observed / random)'); ax.set_title('Resemblance to cancer CNAs', loc='left'); plab(ax, 'c', -0.38)
save(fig, FG, 'FigS10_noncancer_call_annotation', T)
