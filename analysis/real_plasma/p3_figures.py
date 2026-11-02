import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, '/home/claude/sim2_final'); from s2_analysis import style, COLOR, MARK, save, plab, MM, INK, MUTED
from pathlib import Path
from scipy.stats import spearmanr
style(); FG = Path('figures'); FG.mkdir(exist_ok=True)
B = pd.read_csv('tables/p3_burden.tsv', sep='\t'); T = pd.read_csv('tables/p3_tf.tsv', sep='\t')
I = pd.read_csv('tables/p3_insert_sizes.tsv', sep='\t'); R = pd.read_csv('tables/p3_recurrent_windows.tsv', sep='\t')
COL = dict(COLOR); COL['BIC-seq2'] = COLOR['BICseq2']
CFGS = [('QDNAseq', '100kb'), ('QDNAseq', '50kb'), ('Spectre', '100kb'), ('CNVkit', '30kb'), ('ichorCNA', '1Mb low-TF'), ('ichorCNA', '1Mb defaults')]
CFGS = [c for c in CFGS if ((B.tool == c[0]) & (B.config == c[1])).any()]
if ((B.tool == 'BIC-seq2')).any(): CFGS.insert(3, ('BIC-seq2', '10kb'))
lab = lambda t, c: f"{t}\n{c}"
fig = plt.figure(figsize=(180 * MM, 120 * MM)); gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=0.95, wspace=0.62)
rng = np.random.default_rng(1)
# a FGA full depth
ax = fig.add_subplot(gs[0, :2])
for i, (t, c) in enumerate(CFGS):
    v = B[(B.tool == t) & (B.config == c) & (B.level == 'full')].FGA_pct.values
    ax.scatter(i + rng.uniform(-.18, .18, len(v)), v, s=9, color=COL[t], alpha=.7, lw=0)
    ax.plot([i - .28, i + .28], [np.median(v)] * 2, color=INK, lw=1.3)
ax.set_yscale('symlog', linthresh=0.1); ax.set_ylim(0, 100); ax.set_xticks(range(len(CFGS))); ax.set_xticklabels([lab(*x).replace(' defaults','\ndefaults').replace(' low-TF','\nlow-TF') for x in CFGS], fontsize=5.5)
ax.set_ylabel('Genome called altered (%)'); ax.grid(axis='x', visible=False)
ax.set_title('Non-cancer plasma (n = 25), full depth', loc='left'); plab(ax, 'a', -0.09)
# b depth stability
ax = fig.add_subplot(gs[0, 2]); lv = ['d0p1', 'd1', 'full']; xs = [0.1, 1, 4]
for t, c in CFGS:
    m = [B[(B.tool == t) & (B.config == c) & (B.level == l)].FGA_pct.median() for l in lv]
    ax.plot(xs, m, color=COL[t], marker=MARK.get(t.replace('BIC-seq2', 'BICseq2'), 'o'), ls='--' if 'defaults' in c or c == '50kb' else '-', label=f"{t} {c}")
ax.set_xscale('log'); ax.set_xticks(xs); ax.set_xticklabels(['0.1', '1', 'full\n(3–12)']); ax.minorticks_off()
ax.set_yscale('symlog', linthresh=0.1); ax.set_ylim(0, 100); ax.set_xlabel('Depth (×)'); ax.set_ylabel('Median genome altered (%)')
ax.set_title('Depth', loc='left'); plab(ax, 'b', -0.3); ax.legend(fontsize=5, loc='upper left', bbox_to_anchor=(0, -0.32), ncol=2, handlelength=1.5)
# c tumour fraction
ax = fig.add_subplot(gs[1, 0]); M = [('ichorCNA 1000kb_lowTF', 'ichor\nlow-TF'), ('ichorCNA 1000kb_default', 'ichor\ndefaults'), ('ACE ploidy2 100kb', 'ACE\n100 kb'), ('Fragle', 'Fragle')]
ax.axhspan(0, 0.011, color='#e8e8e8', lw=0); ax.text(3.45, 0.012, 'ichorCNA in HG008\n0% mixtures', fontsize=5, ha='right', va='bottom', color=MUTED)
for i, (m, l) in enumerate(M):
    v = T[(T.method == m) & (T.level == 'full')].TF.values * 100
    ax.scatter(i + rng.uniform(-.18, .18, len(v)), v, s=9, color=['#00441b', '#009E73', '#999999', '#CC79A7'][i], alpha=.75, lw=0)
    ax.plot([i - .28, i + .28], [np.median(v)] * 2, color=INK, lw=1.3)
ax.set_xticks(range(len(M))); ax.set_xticklabels([l for _, l in M], fontsize=6); ax.set_ylabel('Estimated tumour fraction (%)'); ax.grid(axis='x', visible=False)
ax.set_title('Tumour-fraction estimates', loc='left'); plab(ax, 'c', -0.3)
# d Fragle vs short fragments
ax = fig.add_subplot(gs[1, 1]); f = T[(T.method == 'Fragle') & (T.level == 'full')].merge(I, on='base')
ax.scatter(f.frac_lt150 * 100, f.TF * 100, s=12, color='#CC79A7', lw=0); rho, p = spearmanr(f.frac_lt150, f.TF)
ax.set_xlabel('Fragments < 150 bp (%)'); ax.set_ylabel('Fragle ctDNA fraction (%)'); ax.text(0.03, 0.95, f'ρ = {rho:.2f}, p = {p:.1g}', transform=ax.transAxes, va='top', fontsize=6.5)
ax.set_title('Fragle vs short fragments', loc='left'); plab(ax, 'd', -0.3)
# e recurrence track
ax = fig.add_subplot(gs[1, 2]); chroms = [f'chr{i}' for i in range(1, 23)]
for k, (t, c) in enumerate([('QDNAseq', '100kb'), ('Spectre', '100kb'), ('BIC-seq2', '10kb'), ('CNVkit', '30kb'), ('ichorCNA', '1Mb low-TF')]):
    r = R[(R.tool == t) & (R.config == c)].groupby(['chrom', 'start']).n_samples.max().reset_index()
    n_rec = (r.n_samples >= 13).sum()
    ax.barh(k, n_rec, color=COL[t]); ax.text(n_rec + 5, k, str(n_rec), va='center', fontsize=6)
ax.set_yticks(range(5)); ax.set_yticklabels(['QDNAseq', 'Spectre', 'BIC-seq2', 'CNVkit', 'ichorCNA'], fontsize=6); ax.invert_yaxis()
ax.set_xlabel('1-Mb windows in ≥13/25 samples'); ax.grid(axis='y', visible=False); ax.set_title('Recurrent calls', loc='left'); plab(ax, 'e', -0.5)
save(fig, FG, 'Fig10_real_noncancer_plasma', B)
