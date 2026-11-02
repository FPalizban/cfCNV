# Changelog

## v1.0.0 (November 2026)
First release, accompanying the manuscript.
- Nextflow DSL2 workflow with four modes: `call`, `simulate` (Simulation 1), `mix` (Simulation 2), `score`
- All eight benchmarked callers: QDNAseq, ichorCNA, BIC-seq2, CNVkit, Spectre, CNV-MEANN, dudeML, LYCEUM; plus ACE and Fragle
- CNV-MEANN and dudeML trained in `--mode simulate` (leave-one-replicate-out) and applied elsewhere with `--ml_models`
- Eight containers (GHCR)
- Harmonised call table, call burden, tumour fraction, majority consensus, caller agreement, depth robustness
- Annotation against gnomAD-SV, segmental duplications, ENCODE blacklist, OncoKB and TCGA PanCanAtlas
- Simulation 1 truth sets and the scoring rule of the paper; analysis scripts for every figure and table
