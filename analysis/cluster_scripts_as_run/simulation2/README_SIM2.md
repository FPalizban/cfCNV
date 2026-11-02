# cfCNV – Simulation 2 (HG008 tumour–normal mixtures, GRCh38): full run

**Input:** the 104 calibrated mixtures already built in `05_simulation/simulation2_HG008/mixtures_v2`.
The design is 0–30 % tumour (including 2 %) at 0.1–10× depth, with 3 replicates at ≤ 1×.
The only mixture missing is 10× / 30 %, which needs more tumour reads than the pool holds.

**Truth:** NIST HG008-T V0.5 draft somatic CNV benchmark (PASS regions).

## What runs

| Step | Tool | Settings (same as Simulation 1 unless noted) | Tasks |
|---|---|---|---|
| A1 | QDNAseq | QDNAseq.hg38 bins, 6 bin sizes | 624 |
| A2 | CNVkit | flat GRCh38 reference, 6 bins (a pooled normal is impossible: one donor, and its normal is inside every mixture) | 624 |
| A3 | ichorCNA | 1 Mb tumour defaults, 1 Mb low-tumour-fraction settings, 500 kb defaults | 312 |
| A4 | Spectre | coverage at the 6 bin sizes, blacklist, `--cancer` | 624 |
| A5 | ACE | tumour fraction from the QDNAseq segmentation, 100 kb / 500 kb / 1 Mb, ploidy 2 and free ploidy | 312 |
| A6 | dudeML | Simulation 1 models of the same depth and bin size (no retraining) | 104 + 624 |
| A7 | CNV-MEANN v2 | features → normalisation → network trained on all Simulation 1 positives | 624 + 6 |
| B | BIC-seq2 | new GRCh38 150-mer mappability (one long job, up to ~1–2 days), then 6 bins | 1 + 104 + 624 |
| C | LYCEUM | **native** GRCh38 run with its own exon targets | 104 |
| D | scoring | runs twice: after A + C, and again after BIC-seq2 | 2 |

Fragle is not run here. It reads cfDNA fragment-length patterns, and these mixtures are tissue and cell-line DNA. It belongs to the real-cfDNA phase.

**How it is scored (because HG008-T has a mean ploidy of ~1.44):**
1. Tumour-fraction estimates (ichorCNA, ACE) compared with the measured fraction.
2. Correlation of each tool's genome profile with the tumour's true copy-number pattern in 1-Mb windows. The 0 % mixtures give the background, from which a detection limit is derived per depth.
3. Events relative to the tumour's modal copy number, with a two-copy baseline also reported.

New disk use: about 150 GB, mostly BIC-seq2 read-position files and LYCEUM arrays.

## Steps

**1. Copy the zip from your Mac to the cluster** (Mac Terminal):
```bash
scp ~/Downloads/cfCNV_sim2_full.zip palizbanf@palizbanf-hpc9.research.chop.edu:~/
```
**2. Log in and unzip:**
```bash
ssh palizbanf@palizbanf-hpc9.research.chop.edu
cd /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/simulation2
unzip -o ~/cfCNV_sim2_full.zip          # creates the folder  simulation2/full
cd full
```
**3. Setup (login node, needs internet, ~15 min, once):**
```bash
bash 0_setup_sim2.sh
```
It must end with **`SETUP_SIM2_OK`**. If it stops with `MISSING`, send me the output.

**4. Submit everything:**
```bash
bash 1_submit_sim2.sh
```
**5. After ~30 min, check that nothing fails early.** Paste me any errors.
```bash
bash 2_check_sim2.sh
```
Check again whenever you like. Most tools finish within about a day. BIC-seq2 waits for its mappability job, so it finishes last.

**6. First results:** once `2_check_sim2.sh` shows `scores: ready`, collect and send:
```bash
bash 3_collect_sim2.sh                  # -> ~/cfcnv_sim2_for_claude.tar.gz
```
On your Mac, then attach the file here:
```bash
scp palizbanf@palizbanf-hpc9.research.chop.edu:~/cfcnv_sim2_for_claude.tar.gz ~/Documents/cfCNV/4Claude/
```
**7. When BIC-seq2 has finished** (`squeue -u $USER` is empty), the second scoring runs by itself. Run step 6 again and send the new file.

Re-submitting any step is safe, because finished tasks are skipped. Example:
`sbatch --array=0-623%40 A1_qdnaseq.slurm`. To re-score at any time: `sbatch D_score_sim2.slurm`.
