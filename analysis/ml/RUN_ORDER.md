# cfCNV — remaining tools: run order

## 0. Transfer and set up (once)

On your Mac:
```bash
scp cfCNV_ML_scripts.zip palizbanf@<HPC_HOST>:/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/
```
On the HPC:
```bash
export PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project
cd $PROJECT/00_metadata
unzip -j -o cfCNV_ML_scripts.zip && chmod +x *.sh
source config.sh && source 00_env.sh
mkdir -p $PROJECT/logs/{tf,hg008} $PROJECT/00_metadata/manifests/{tf,hg008}
export R_MODULE=<the R module you used for QDNAseq>
```

---

## 1. dudeML (already running)

```bash
bash 80_check_ml_jobs.sh dude_                       # status
scontrol update JobId=25046738 ArrayTaskThrottle=12  # optional: speed up features
```
When finished, inspect one prediction file and **send me the output**:
```bash
python 81_harmonize_dudeml.py --inspect <one prediction file>
```

---

## 2. CNV-MEANN

```bash
# 2a. environment
module load Python/3.10.8-GCCcore-12.2.0
python -m venv $PROJECT/06_tools/envs/cnvmeann
source $PROJECT/06_tools/envs/cnvmeann/bin/activate && pip install numpy pandas pysam && deactivate

# 2b. per-sample manifest (expect 115 lines)
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; print "sample_id\tbam\tbackground\ttarget_depth\treplicate\tclass"; next}
  !s[$h["sample_id"]]++{print $h["sample_id"]"\t"$h["bam"]"\t"$h["background"]"\t"$h["target_depth"]"\t"$h["replicate"]"\t"$h["class"]}' \
  $SIM1_MAN/bicseq2_tasks.tsv > $SIM1_MAN/sim1_samples.tsv
wc -l $SIM1_MAN/sim1_samples.tsv

# 2c. features: test one, then all
sbatch --array=0 91_cnvmeann_features.slurm
sbatch --array=1-113%10 91_cnvmeann_features.slurm

# 2d. MATLAB code dump -> send me cnvmeann_matlab_dump.txt
cd $PROJECT/06_tools/ML_based && git clone https://github.com/huangtihao/CNV-MEANN.git; cd CNV-MEANN
for f in main.m bpcnv_train.m bpcnv_test_simulation.m; do echo "===== $f"; iconv -f GBK -t UTF-8 "$f" 2>/dev/null || cat "$f"; done > ~/cnvmeann_matlab_dump.txt
cd $PROJECT/00_metadata
```

---

## 3. Tumour-fraction tools: ACE + Fragle

```bash
# 3a. install (login node)
bash 99_install_tf_tools.sh

# 3b. manifest: sample_id bam dataset build true_tf target_depth
#     start with the 19 Sim-1 negatives (true TF = 0)
TF=$PROJECT/00_metadata/manifests/tf/tf_samples.tsv
awk -F'\t' 'NR==1{print "sample_id\tbam\tdataset\tbuild\ttrue_tf\ttarget_depth"; next}
  $6!="positive"{print $1"\t"$2"\tSIM1\thg19\t0\t"$4}' $SIM1_MAN/sim1_samples.tsv > $TF
#     later append: HG008 mixtures (dataset HG008, build hg38, true_tf 0..0.30)
#                   PRJNA795275 BAMs (true_tf NA)
N=$(( $(wc -l < $TF) - 2 ))

# 3c. run
sbatch --array=0-$N%10 98_ace.slurm
sbatch --array=0-$N%10 99_fragle.slurm      # skips HG008 rows by design
bash 80_check_ml_jobs.sh tf_

# 3d. collect + error table
python 100_collect_tf.py --manifest $TF --tfroot $PROJECT/07_tool_outputs/tf \
  --outdir $PROJECT/09_metrics/tf
```

---

## 4. LYCEUM (after HG008 mixtures exist)

```bash
bash 95_lyceum_setup.sh                    # install + toy test; send me the label values it prints
# manifest: sample_id bam tumor_fraction target_depth -> manifests/hg008/hg008_mixtures.tsv
sbatch --array=0-<N-1>%8 96_lyceum_hg008.slurm
python 97_harmonize_lyceum.py --root $PROJECT/07_tool_outputs/hg008/LYCEUM \
  --manifest $PROJECT/00_metadata/manifests/hg008/hg008_mixtures.tsv \
  --outdir $PROJECT/08_harmonized_calls/hg008/LYCEUM
```

---

## Send me

1. dudeML `--inspect` output
2. `~/cnvmeann_matlab_dump.txt`
3. Log of the first CNV-MEANN, ACE and Fragle task (`$PROJECT/logs/...`)
