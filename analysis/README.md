# Analysis code used for the cfCNV paper

The Nextflow workflow in the repository root re-runs the benchmark. This folder holds what was used to
turn the runs into the paper's tables and figures, plus the scripts exactly as they were run on the cluster.

| Folder | What it reproduces |
|---|---|
| `figures/` | Figures 2–5 and S1, S9, S12–S14 from the source data in `figures/data/` (`python figN_*.py`; see its README) |
| `simulation1/` | Simulation 1 scoring (event F1, base-pair and bin metrics, background mask), Tables 2 and S2 |
| `simulation2/` | HG008 mixture scoring: profile correlation with the NIST HG008 somatic benchmark, limit of detection, tumour-fraction accuracy |
| `real_plasma/` | Non-cancer plasma (PRJNA795275) burden and concordance (`p3_analysis.py`), cancer plasma (GSE71378) evaluation (`g7_analysis.py`, `g7_eval.py`), annotation figures |
| `ml/` | CNV-MEANN, dudeML and LYCEUM scripts as run for the paper (the workflow versions are in `bin/` and `modules/local/ml.nf`) |
| `cluster_scripts_as_run/` | The original Slurm scripts, unchanged, with the paths of the cluster they ran on |
| `record_tool_versions.sh` | Prints the exact tool versions used on the cluster (for the Methods / supplementary table) |

**Machine-learning callers.** CNV-MEANN, dudeML and LYCEUM are part of the workflow (`--mode simulate` trains
CNV-MEANN and dudeML exactly as for Simulation 1; `--ml_models` applies the trained models in the other modes, as for
Simulation 2). `ml/` keeps the scripts as they were run for the paper, including the MATLAB-based CNV-MEANN training
used for Simulation 1 (MATLAB cannot be distributed in a container).

**Paper result -> workflow output.** The workflow produces the same per-run tables the paper was built on:
`summary/event_scores.tsv` (Simulation 1, `--mode simulate`; same matching rule, reproduces the published F1 values
when given the same calls), `summary/call_burden.tsv` and `summary/tumour_fraction.tsv` (all modes),
`summary/caller_agreement.tsv`, `consensus_*.tsv`, `depth_robustness.tsv` (cancer evaluation, `--mode call`),
and `summary/annotation/` (germline/artefact and TCGA annotation).
