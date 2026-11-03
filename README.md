# cfCNV — benchmarking copy-number callers on cell-free DNA

[![CI](https://github.com/FPalizban/cfCNV/actions/workflows/ci.yml/badge.svg)](https://github.com/FPalizban/cfCNV/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.04-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)

cfCNV is a benchmark and a reproducible Nextflow workflow for copy-number variant and alteration (CNV/CNA) calling
in low-pass cell-free DNA (cfDNA) whole-genome sequencing. It runs all eight benchmarked callers (five classical, three
machine-learning) and three tumour-fraction estimators in containers, trains the machine-learning callers on simulated
data, harmonises their output, builds a majority consensus, annotates calls against germline,
artefact and cancer references, and scores calls against simulated truth.

**Interactive explorer:** <https://fpalizban.github.io/cfCNV/> — results, a read-depth simulator and in-browser scoring of your own calls.

| Layer | Data | Truth | Question |
|---|---|---|---|
| Simulation 1 | Healthy cfDNA (IH01–IH03, Snyder et al. 2016) with inserted gains and losses, 0.1–1× (up to 20×) | 30 events per replicate, 50 kb–50 Mb, 5 replicates | Can the caller find events of known size and position? |
| Simulation 2 | NIST HG008 tumour mixed into matched normal, 0–30% tumour, 0.1–10× | NIST HG008 somatic copy-number benchmark | How low a tumour fraction can it detect? Is the tumour fraction right? |
| Real plasma | 25 non-cancer WGS plasma samples (PRJNA795275); 5 stage IV cancers + 1 healthy (GSE71378) | none expected / consensus, TCGA PanCanAtlas | How many false alterations in real plasma? Do cancer calls look like the tumour type? |

## Main results

| Caller (primary setting) | Sim 1 event F1 (0.1–1×) | HG008 limit of detection (0.25–1×) | Non-cancer plasma: genome altered | Cancer plasma: F1 vs other callers / TCGA r |
|---|---|---|---|---|
| QDNAseq (50 kb / 100 kb) | 0.71 | 1–2% tumour | 0% | 0.86 / 0.48 |
| ichorCNA (1 Mb, low-TF) | 0.46 | 2% | 26% (tumour fraction ≈1%) | 0.50 / 0.47 |
| BIC-seq2 (30 kb / 10 kb) | 0.45 | 1–2% | 1.0% | 0.65 / 0.43 |
| CNVkit (pooled / flat reference) | 0.34 | 1%* | 6.8% | 0.56 / 0.37 |
| Spectre (500 kb / 100 kb) | 0.23 | 20–30% | 0.4% | 0.16 / 0.29 |
| CNV-MEANN, dudeML (trained on Sim 1) | 0.45, 0.40 | 10–20% | not run | not run |
| LYCEUM | 0 | not detected | not run | not run |

\* narrow null in pure-normal mixtures; CNVkit's profile correlation was the lowest of the classical callers.

---

## Quick start

Requirements: [Nextflow](https://www.nextflow.io/docs/latest/install.html) ≥ 24.04 (Java 17+) and either
Singularity/Apptainer (HPC) or Docker. All tools come in containers that are downloaded automatically.

```bash
# 1. check the installation (no data, no containers, about a minute)
nextflow run FPalizban/cfCNV -r v1.0.0 -profile test -stub

# 2. call CNAs in your BAMs with the paper's real-data settings, plus 1x and 0.1x down-sampling
nextflow run FPalizban/cfCNV -r v1.0.0 -profile singularity,slurm \
    --input samples.csv --fasta /ref/human_g1k_v37.fasta --genome GRCh37 \
    --downsample 1,0.1 --annotation_refs /refs/cfcnv_annotation/prepared \
    --bicseq2_container /images/cfcnv-bicseq2.sif \
    --bicseq2_resources /refs/bicseq2/GRCh37_101mer --bicseq2_map_prefix human_g1k_v37.101mer \
    --outdir results
```

`samples.csv` (header required; `bai`, `group`, `expected_tf`, `tcga_cohort` are optional):

```
sample,bam,bai,group,expected_tf,tcga_cohort
IC17,/data/IC17.bam,/data/IC17.bam.bai,cancer,,LIHC
IC37,/data/IC37.bam,/data/IC37.bam.bai,cancer,,COAD+READ
IH03,/data/IH03.bam,/data/IH03.bam.bai,control,0,
```

`tcga_cohort` is a TCGA study code (or several joined with `+`); it switches on the TCGA concordance score.

## Modes

| `--mode` | What it does | Key options |
|---|---|---|
| `call` (default) | Callers + tumour-fraction estimators on your BAMs (CNV-MEANN/dudeML with `--ml_models`), optional down-sampling, harmonised calls, consensus, caller agreement, depth robustness, annotation | `--input`, `--fasta`, `--genome`, `--downsample`, `--annotation_refs` |
| `simulate` | **Simulation 1.** Inserts 30 single-copy gains/losses (50 kb–50 Mb) per replicate into healthy cfDNA BAMs, down-samples to each depth, adds unspiked negative controls, trains CNV-MEANN and dudeML (leave-one-replicate-out), runs all callers and scores them; trained models are saved to `models/` | `--input` (background BAMs), `--sim_replicates 5`, `--sim_depths 0.1,0.25,0.5,1`, `--truth assets/truth` |
| `mix` | **Simulation 2.** In-silico tumour–normal mixtures at set depths and tumour fractions, then calling (CNV-MEANN/dudeML with the Simulation 1 models, `--ml_models`) and tumour-fraction estimation | `--tumour_bam`, `--normal_bam`, `--tumour_depth`, `--normal_depth`, `--mix_fractions`, `--mix_depths` |
| `score` | Scores one call set (BED/TSV: `chrom start end type`) against a truth set with the paper's matching rule | `--calls`, `--truth`, `--fasta` (for genome size) |

With GRCh37 sizes, the truth generator reproduces the published Simulation 1 truth sets exactly (`assets/truth/`, checked in CI).

## Callers and settings

| Tool | Version | Setting used (default) | Option |
|---|---|---|---|
| QDNAseq | 1.34 (+ QDNAseq.hg19 / .hg38) | 100-kb bins, CBS, `callBins` | `--qdnaseq_bins_kb 100` (comma list allowed) |
| ichorCNA | 0.2.0 | 1-Mb windows, low-tumour-fraction settings (normal 0.95–0.999, ploidy 2, no subclones) | `--ichor_modes lowTF,default` |
| BIC-seq2 | norm 0.2.4 / seg 0.7.2 | 10-kb bins, λ = 2, call when \|log2\| ≥ 0.2 | `--bicseq2_bin`, `--bicseq2_lambda`, `--bicseq2_log2_threshold` |
| CNVkit | 0.9.10 | WGS mode, 30-kb bins, flat reference (pooled if normals are given) | `--cnvkit_bin`, `--cnvkit_normals 'n1.bam,n2.bam'` |
| Spectre | 0.2.1 | mosdepth 100-kb coverage, `--cancer`, min. CNV length max(bin, 50 kb) | `--spectre_bin`, `--blacklist` |
| CNV-MEANN | Python port of the published network | 500-kb bins; 4 features → 15 tanh units → 3 classes; trained by the workflow | `--cnvmeann_bin` |
| dudeML | commit a4c7c0b | 1-Mb windows, 5 flanking windows, random forest (RFC100) per depth; trained by the workflow | `--dudeml_bin` |
| LYCEUM | commit a5a8b07, released weights | hg38 exon targets (lifted one-to-one to GRCh37 when needed), exon calls merged into segments | – |
| ACE | Bioconductor | cellularity at ploidy 2 and free ploidy, from QDNAseq segments | `--tf_methods` |
| Fragle | commit 3212220 | raw WGS mode (paired-end BAMs) | `--tf_methods` |

**Machine-learning callers.** LYCEUM is pretrained and runs in every mode. CNV-MEANN and dudeML need labelled
training data: `--mode simulate` trains them (each CNV-positive sample is predicted by models that never saw its
replicate; negative controls by all models) and writes the models to `results/models/`. Pass that folder as
`--ml_models` to apply them in `call` or `mix` mode, as was done for Simulation 2 in the paper (dudeML uses the models
trained at the depth closest to each sample). Without `--ml_models` these two callers are skipped with a warning. In
the paper, CNV-MEANN for Simulation 1 was trained with the authors' MATLAB code, which needs a MATLAB licence; the
workflow uses the Python implementation of the same network that was used for Simulation 2.

Simulation 1 primary settings were QDNAseq 50 kb, BIC-seq2 30 kb, CNVkit 100 kb pooled and Spectre 500 kb; pass them
with the options above. BIC-seq2 needs per-chromosome FASTA and mappability files for the read length; without
`--bicseq2_resources` the workflow builds them once (GEM, ~1–2 days) and stores them under `results/resources`.

## Output

```
results/
├── callers/<tool>/<sample>/        native output of every caller (VCF, .cns, .seg.txt, .cnv, plots)
├── summary/
│   ├── all_calls.tsv              every call, harmonised: sample_id base level group tool config chrom start end type
│   ├── call_burden.tsv            calls, gains, losses, % genome altered, calls >= 10 Mb per caller run
│   ├── tumour_fraction.tsv        ichorCNA, ACE (ploidy 2 / free), Fragle estimates (+ expected TF in mix mode)
│   ├── consensus_calls.tsv        majority consensus (>= 3 callers) in 1-Mb windows
│   ├── caller_agreement.tsv       each caller vs the other callers (precision, recall, F1) and TCGA r
│   ├── depth_robustness.tsv       F1 of each down-sampled profile vs the same caller at native depth
│   ├── window_states.tsv.gz       gain/neutral/loss per 1-Mb window, caller and sample
│   ├── event_scores.tsv           (simulate/score) TP, FP, FN, precision, recall, F1, false calls/Gb, bp metrics
│   ├── truth_detection.tsv        (simulate/score) detected yes/no per truth event
├── models/                        (simulate) trained CNV-MEANN and dudeML models, for --ml_models
│   └── annotation/                calls vs gnomAD-SV, segmental duplications, blacklist, OncoKB, TCGA; recurrent hotspot mask
└── pipeline_info/                 run report, timeline, trace (runtime and memory per task)
```

## Annotation references

```bash
bash bin/get_annotation_refs.sh /refs/cfcnv_annotation      # downloads public data, ~15 min
# -> pass /refs/cfcnv_annotation/prepared as --annotation_refs
```

TCGA PanCanAtlas arm calls and SNP6 segments, GENCODE v19, OncoKB cancer genes, gnomAD-SV v2.1 (AF ≥ 1%), UCSC
segmental duplications, ENCODE blacklist v2 and the hg38→hg19 chain (GRCh38 calls are lifted over piecewise).

## Running on an HPC

`-profile singularity,slurm` works on most Slurm clusters; `conf/chop.config` is an example institutional profile
(`-profile chop`). Set `NXF_SINGULARITY_CACHEDIR` to a large shared folder so images are downloaded once, start the
run inside `screen`/`tmux` or as a small Slurm job, and add `-resume` to continue after an interruption.
Resource limits: `--max_cpus`, `--max_memory`, `--max_time`.

## Containers

Built from `containers/*/Dockerfile` by GitHub Actions and published at `ghcr.io/fpalizban/cfcnv-{r,cnvkit,spectre,fragle,ml,lyceum,tools}:1.0.0`. Two tools are not inside any published image because their terms do not clearly allow redistribution: **BIC-seq2** (build your own image, `assets/hpc/README_bicseq2.md`, and pass `--bicseq2_container`; otherwise BIC-seq2 is skipped) and **dudeML** (downloaded from the authors' repository when the workflow runs).
To use your own image for a step, add a config file with e.g.
`process { withName: 'QDNASEQ' { container = '/path/my_image.sif' } }` and pass it with `-c my.config`.

## Reproducing the paper

`analysis/` contains the scripts that turned the runs into the tables and figures (`analysis/figures` regenerates
Figures 2–5 from source data), the machine-learning caller scripts, and the cluster scripts exactly as run. See
`analysis/README.md`.

## Citation

If you use cfCNV, please cite the paper (see `CITATION.cff`) and the tools it runs: QDNAseq (Scheinin et al. 2014),
ichorCNA (Adalsteinsson et al. 2017), BIC-seq2 (Xi et al. 2016), CNVkit (Talevich et al. 2016), Spectre
(github.com/fritzsedlazeck/Spectre), CNV-MEANN, dudeML, LYCEUM, ACE (Poell et al. 2019), Fragle (Zhu et al., bioRxiv 2023, doi:10.1101/2023.07.28.550922), mosdepth (Pedersen & Quinlan 2018),
SAMtools (Danecek et al. 2021), Nextflow (Di Tommaso et al. 2017).

## Licence

**No licence has been chosen yet for the cfCNV code.** Until one is added, all rights are reserved: the code can be
read here but not reused or redistributed. A licence will be added after institutional review.

The tools that cfCNV runs keep their own licences, and some restrict use: BIC-seq2 is free for academic use only,
Fragle is CC BY-NC-SA 4.0 (non-commercial), and dudeML has no licence file. See
[THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md) for every tool, its licence and how cfCNV obtains it.
