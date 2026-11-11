# cfCNV ML dataset + baseline model

Builds a Hugging Face-style training/test dataset from the BAMs you already have (Simulation 1,
Simulation 2/HG008, PRJNA795275, GSE71378) and trains two simple baselines on it.

* Scripts live in `$P/00_metadata/ml_dataset`
* Everything new is written only to `$P/13_ml_dataset` and `$P/logs/ml_dataset`
* BAMs, truth files, references and the old manifests are only read. Nothing is deleted or overwritten
  (an older dataset folder is renamed, not removed).

| step | file | what it does |
|---|---|---|
| m0 | m0_make_manifest.py | one table of all BAMs + the fixed train/validation/test split |
| A | A_gc_bins.slurm | GC / N / blacklist per 10-kb bin for GRCh37 and GRCh38 |
| B | B_count_bins.slurm | read + fragment-length counts per 10-kb bin, one array task per BAM |
| C | C_assemble.slurm | features + labels -> `$P/13_ml_dataset/hf_dataset` (data + dataset card) |
| D | D_baseline.slurm | threshold baseline + gradient-boosting baseline, all scores, small bundle to send back |

`bash 1_run_all.sh` does m0 and submits A-D with dependencies. `bash 2_check.sh` shows progress.

Splits: train = Sim1 IH01+IH02 replicates 1-4 (+ negatives); validation = Sim1 IH01+IH02 replicate 5;
test_sim1_ih03; test_hg008; test_noncancer (PRJNA795275); test_cancer (GSE71378, weak consensus labels).

The result to send back is `$P/13_ml_dataset/ml_results_for_claude.tar.gz` (tables only, a few MB).
Before publishing the dataset anywhere, check the data-use terms of each source and choose a licence.

## Model v2 (augmented training)

`bash 3_run_v2.sh` submits step E (`E_augmented.slurm`, `e_augmented.py`). It makes 300 training + 60 validation
samples from the Simulation 1 negative controls of the TRAIN split, with random CNV layouts diluted to 3-100 %
tumour fraction and thinned to 0.1x-native depth, retrains, tunes on the augmented validation set only, and scores
every split. Results: `$P/13_ml_dataset/results_v2/`, bundle `$P/13_ml_dataset/ml_results_v2_for_claude.tar.gz`.

## Final run (model v2 + manuscript tables and figure)

`bash 4_final.sh` submits E (model v2) and F (report; waits for E). F scores both models against the 8 callers with
the manuscript definitions (HG008 1-Mb profile r, strict LoD, median r 1-10 %; cancer F1 vs consensus on 1-Mb windows;
non-cancer % genome altered) and writes `$P/13_ml_dataset/manuscript/`:
tables (`table_*.tsv`), per-sample data, `figure_ml_benchmark.pdf/png/svg`, `key_numbers.md`, `methods_notes.md`.
Bundle: `$P/13_ml_dataset/cfCNV_ML_manuscript_results.tar.gz`. The report prints a check that QDNAseq and ichorCNA
reproduce the manuscript LoD and median r exactly.
