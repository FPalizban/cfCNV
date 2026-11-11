---
license: other
license_name: to-be-decided
pretty_name: "cfCNV-bins: bin-level cfDNA low-pass WGS copy-number benchmark"
language:
- en
task_categories:
- tabular-classification
- tabular-regression
tags:
- biology
- genomics
- medical
- cell-free-DNA
- liquid-biopsy
- copy-number-variation
- copy-number-alteration
- whole-genome-sequencing
- benchmark
size_categories:
- 10M<n<100M
configs:
- config_name: bins_100kb
  default: true
  data_files:
  - split: train
    path: data/bins_100kb/train.csv.gz
  - split: validation
    path: data/bins_100kb/validation.csv.gz
  - split: test_sim1_ih03
    path: data/bins_100kb/test_sim1_ih03.csv.gz
  - split: test_hg008
    path: data/bins_100kb/test_hg008.csv.gz
  - split: test_noncancer
    path: data/bins_100kb/test_noncancer.csv.gz
  - split: test_cancer
    path: data/bins_100kb/test_cancer.csv.gz
- config_name: bins_1Mb
  data_files:
  - split: train
    path: data/bins_1Mb/train.csv.gz
  - split: validation
    path: data/bins_1Mb/validation.csv.gz
  - split: test_sim1_ih03
    path: data/bins_1Mb/test_sim1_ih03.csv.gz
  - split: test_hg008
    path: data/bins_1Mb/test_hg008.csv.gz
  - split: test_noncancer
    path: data/bins_1Mb/test_noncancer.csv.gz
  - split: test_cancer
    path: data/bins_1Mb/test_cancer.csv.gz
- config_name: augmented_100kb
  data_files:
  - split: train
    path: data/augmented_100kb/train.csv.gz
  - split: validation
    path: data/augmented_100kb/validation.csv.gz
---

# cfCNV-bins

Bin-level read-depth and fragment-length features from low-coverage whole-genome sequencing of
cell-free DNA (cfDNA), with copy-number labels and fixed train / validation / test splits.
It is the machine-learning companion of the **cfCNV** benchmark of copy-number callers on cfDNA
([code and workflow](https://github.com/FPalizban/cfCNV), [interactive explorer](https://fpalizban.github.io/cfCNV/)).

The dataset is meant for training and testing models that call copy-number gains and losses from
low-pass cfDNA sequencing, and for comparing them with the eight callers of the benchmark on the same samples.

## At a glance

| | |
|---|---|
| Samples | 310 sequencing libraries (BAM files) from four arms |
| Unit | one row = one genomic bin of one sample |
| Resolutions | 100 kb (default) and 1 Mb; autosomes only |
| Labels | loss / neutral / gain per bin (`label_state` = -1 / 0 / +1) |
| Reference builds | GRCh37 (Simulation 1, GSE71378) and GRCh38 (HG008, PRJNA795275); see column `build` |
| Format | gzip-compressed CSV, one file per split |

## Splits

The splits were fixed before any model was trained and are part of the benchmark definition.
All tuning must use `train` and `validation` only.

| split | samples | what it is | labels |
|---|---|---|---|
| `train` | 75 | Simulation 1: CNVs spiked into healthy cfDNA of individuals IH01 and IH02, replicates 1-4, plus negative controls | exact (simulated) |
| `validation` | 15 | Simulation 1: IH01 and IH02, replicate 5 (a CNV layout not present in `train`) | exact (simulated) |
| `test_sim1_ih03` | 24 | Simulation 1: an individual not seen in training (IH03) | exact (simulated) |
| `test_hg008` | 104 | in-silico mixtures of the GIAB HG008 tumour and matched normal, tumour fraction 0-30 %, depth 0.1x-10x | NIST HG008-T draft somatic CNV benchmark |
| `test_noncancer` | 75 | plasma of 25 individuals without cancer (PRJNA795275) at native depth, 1x and 0.1x | assumed neutral everywhere |
| `test_cancer` | 17 | cancer plasma and one healthy control (GSE71378) at several depths | consensus of five callers (weak labels) |

`samples.tsv` lists every sample with its split, depth, tumour fraction, noise level and number of gain / loss bins.

The third config, `augmented_100kb`, holds the synthetic samples used to train the augmented baseline (model v2):
360 samples generated from the Simulation 1 negative controls of the `train` split by thinning counts to lower
depth and adding random copy-number layouts diluted to a tumour fraction of 3-100 %. It contains no test data.

## Columns

**Sample:** `dataset`, `split`, `sample_id`, `base_sample`, `class`, `replicate`, `depth_x`,
`tumour_fraction_nominal`, `tumour_fraction_realised`, `build`.
For `test_cancer`, `tumour_fraction_realised` is the ichorCNA estimate of the full-depth sample, not a ground truth.

**Bin:** `chrom` (1-22, without "chr"), `start`, `end` (0-based, half-open, in the coordinates of `build`).

**Annotation:** `gc` (GC fraction), `map_proxy` (share of fragments with MAPQ >= 20, pooled over all samples of the
build), `frac_usable` (fraction of the bin that passed the filters).

**Counts:** `n_q20` (fragments with MAPQ >= 20), `n_frag` (fragments with known length), `n_short` (<= 150 bp),
`n_long` (> 150 bp).

**Features:** `lr` (GC-corrected log2 ratio, sample median 0), `lr_med3`, `lr_med9`, `lr_med25` (rolling medians
within the chromosome, in bins), `lr_mean9`, `lr_chrom` (chromosome median), `noise` (one value per sample),
`lr_z9`, `log_depth_bin`, `short_frac_lr`, `short_frac_med9`, `flen_dev`.
Fragment-length features are empty for single-end libraries.

**Labels:** `label_state` (-1 loss, 0 neutral, +1 gain, empty = unknown), `label_cn` (copy number where known),
`truth_log2` (expected log2 ratio of the pure altered genome), `event_id`, `event_size` (Simulation 1 only),
`label_source`.

The `augmented_100kb` config has the sample and bin columns, `label_state`, `label_cn`, and the read-depth
features (including `lr_med51`, `lr_med101` and noise-scaled versions `lr_z3` ... `lr_z101`).

## How the data were made

1. **Counting.** One fragment is counted per single-end read or per read 1 of a pair. Unmapped, secondary,
   supplementary, QC-fail and duplicate reads are skipped. Counts are taken in 10-kb bins on autosomes.
2. **Filtering.** A 10-kb bin is dropped when more than 5 % of it is N, when more than 10 % overlaps the ENCODE
   blacklist (GRCh38 only), or when `map_proxy` is below 0.75. The remaining bins are summed to 100 kb and 1 Mb;
   bins with less than half of their length usable are dropped.
3. **Normalisation.** Counts are divided by the median count of bins with similar GC content (20 groups),
   log2-transformed, and centred on the sample median.
4. **Labels.**
   - Simulation 1: a bin is labelled with a spiked-in event when the event covers at least half of the bin.
   - HG008: copy number from the NIST benchmark; gain / loss is relative to the modal copy number of the tumour
     (most of the HG008-T genome is at one copy). Bins with less than half benchmark coverage are unlabelled.
     Mixtures without tumour DNA are labelled neutral. The label describes the tumour genome; the strength of
     the signal in a mixture scales with its tumour fraction.
   - PRJNA795275: neutral everywhere (individuals without cancer).
   - GSE71378: majority consensus of five read-depth callers on 1-Mb windows.

The scripts that build the dataset and the baselines are in `analysis/ml_dataset/` of the [cfCNV repository](https://github.com/FPalizban/cfCNV).

## Baselines

Two simple models are provided as reference points (tables in `baselines/`):
a threshold on the smoothed log2 ratio, and a gradient-boosting classifier (scikit-learn), each trained
once on Simulation 1 only (**v1**) and once with tumour-fraction and depth augmentation (**v2**).
Scores use the definitions of the cfCNV benchmark, so they are directly comparable with the eight callers.

| | HG008 mean limit of detection (tumour fraction) | HG008 median profile r at 1-10 % tumour | non-cancer plasma, % genome called (native depth) | cancer plasma, median F1 vs consensus |
|---|---|---|---|---|
| Gradient boosting v1 | 18.6 % | 0.02 | 0.25 | 0.59 |
| Gradient boosting v2 (augmented) | 4.3 % | 0.17 | 5.4 | 0.46 |
| GC-corrected signal, no model | 1.1 % | 0.30 | - | - |
| QDNAseq (best classic caller) | 1.7 % | 0.51 | 0.0 | 0.86 |

A model trained only on full-strength simulated CNVs does not detect diluted tumour signal; augmentation recovers
most of the sensitivity at a cost in specificity. Beating the read-depth callers on `test_hg008` while keeping
`test_noncancer` clean is the open challenge of this dataset.

## Loading

```python
from datasets import load_dataset

ds = load_dataset("FahimehPalizban/cfCNV-bins", "bins_100kb")        # or "bins_1Mb", "augmented_100kb"
train = ds["train"].to_pandas()
```

or without the `datasets` library:

```python
import pandas as pd
train = pd.read_csv("hf://datasets/FahimehPalizban/cfCNV-bins/data/bins_100kb/train.csv.gz", dtype={"chrom": str})
```

## Limitations

- Real tumour truth comes from a single tumour (HG008), and its benchmark is a draft release.
- Labels of `test_cancer` are the consensus of the callers being benchmarked, not an independent truth, and
  the cohort is small.
- `test_noncancer` is assumed to be copy-number neutral; small germline CNVs are not labelled.
- Simulation 1 CNVs are spiked in at full strength and carry no realistic fragment-length signal; the HG008
  mixtures are made from cell-line DNA, not plasma. Fragment-length features are therefore informative only in
  the two real-plasma test sets.
- Two reference builds are used; coordinates must not be compared across builds.

## Source data

All values are aggregate counts per genomic bin; no reads, genotypes or sequence are included.

- Healthy cfDNA backgrounds and cancer plasma: Snyder et al., *Cell* 2016, NCBI GEO **GSE71378**.
- Non-cancer plasma: NCBI BioProject **PRJNA795275**.
- HG008 tumour / normal and the somatic CNV benchmark: Genome in a Bottle Consortium, NIST (HG008-T, draft benchmark V0.5).

Please cite the original data sources together with this dataset.

## Licence

To be decided by the authors before public release.

## Citation

Palizban F. *et al.* cfCNV: a benchmark of copy-number callers for cell-free DNA whole-genome sequencing.
Manuscript in preparation. (This entry will be updated with the publication and DOI.)

## Contact

Questions and issues: https://github.com/FPalizban/cfCNV/issues
