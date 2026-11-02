# cfCNV figures — regenerate locally

Every figure is drawn from the data files in `data/` by the scripts here. Nothing is hard-coded
except labels: all numbers (medians, limits of detection, summary table) are computed when you run the
scripts, and each main figure also writes its source data to `output/FigureN_source_data.tsv`
(useful for the journal's source-data request and for checking values against the text).

## Install (once)

Python 3.9 or newer:

```bash
pip install pandas numpy matplotlib scipy
```

## Run

```bash
cd cfCNV_figures
python fig2_simulation1.py      # Figure 2  Simulation 1
python fig3_hg008.py            # Figure 3  HG008 tumour-normal mixtures
python fig4_real_plasma.py      # Figure 4  real plasma (a-c non-cancer plasma; d-f stage IV cancer plasma, GSE71378)
python fig5_cost_summary.py     # Figure 5  runtime, memory, summary scorecard
python figS_supplementary.py    # Figures S1, S9, S12, S13, S14
```

Output goes to `output/` as PDF (vector, for submission), SVG (editable in Illustrator/Inkscape) and PNG (600 dpi).
Figure 1 is the BioRender schematic and is not generated here.

## Change the look

Edit **`style.py`** only — every script imports it:

- `TOOL_COLORS` — one colour per caller (used consistently in all figures)
- `HEATMAP_CMAP` — colour map for heatmaps (any matplotlib name, e.g. `"magma"`, `"Blues"`)
- `SCORE_CMAP` — worse → neutral → better colours of the summary table (Figure 5c)
- `BACKGROUND_COLORS` — IH01/IH02/IH03 in Figure S1
- `FONT_FAMILY`, `BASE_SIZE`, `FULL_WIDTH`, `DPI`

Panel layout, axis ranges and labels are in each figure script and are commented by panel letter.

## Data files

| File | Content | Source |
|---|---|---|
| `sim1_per_config.tsv` | Simulation 1 metrics per caller × bin × sample | scoring of all 114 BAMs |
| `sim1_truth_detection.tsv` | detection of every truth event | same |
| `sim1_primary_configs.tsv` | primary bin per caller | highest median F1 at 0.1–1× |
| `sim1_measured_depth.tsv` | mosdepth autosomal depth per BAM | |
| `sim2_profile_correlation.csv` | profile r per mixture × caller × configuration | HG008 scoring |
| `sim2_tf_estimates.csv` | ichorCNA and ACE tumour-fraction estimates | |
| `sim2_mixtures_realised.tsv` | nominal and realised tumour fraction per mixture | |
| `sim2_primary_configs.tsv` | primary configuration per caller | lowest mean LoD |
| `healthy_*.tsv` | PRJNA795275 burden, calls, tumour fraction, fragment sizes, recurrence, annotation | |
| `runtime_memory_summary.tsv` | Slurm wall time and peak memory per caller | |

The limit of detection is recomputed from `sim2_profile_correlation.csv` in `common.py` (`lod_strict`),
with the rule given in the Methods.

## Cancer-plasma data files (Figure 4d-f, Figure 5c)
- `cancer_window_states.tsv` + `grch37_1Mb_windows.tsv`: gain (+1) / neutral (0) / loss (-1) per 1-Mb GRCh37 window, per sample and caller, plus the majority consensus.
- `cancer_loo_and_tcga.tsv`: F1/precision/recall of each caller against the consensus of the other four, and Pearson r with the TCGA frequency of the matching cancer type.
- `healthy_tcga_baseline.tsv`: the same r computed in the 25 non-cancer samples (boxes in 4e).
- `cancer_depth_robustness.tsv`: F1 of the 1x and 0.1x profile against the caller's own native-depth profile (4f).
- `cancer_tumour_fraction.tsv`, `cancer_consensus_arm_events.tsv`, `cancer_gene_direction.tsv`: values quoted in the text.
Gain/loss colours are `GAIN_COLOR` / `LOSS_COLOR` in `style.py`.
