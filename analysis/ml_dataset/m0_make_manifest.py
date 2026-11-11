#!/usr/bin/env python
"""m0  One table of every BAM that goes into the ML dataset (reads the four existing manifests).

Output: $MLO/manifest.tsv   (one row per BAM; row number = Slurm array index)
Splits (decided BEFORE any model is trained, never changed afterwards):
  train            Simulation 1, backgrounds IH01 + IH02, replicates 1-4 and their negative controls
  validation       Simulation 1, backgrounds IH01 + IH02, replicate 5  (new CNV layout, used to tune)
  test_sim1_ih03   Simulation 1, background IH03 (unseen individual)
  test_hg008       Simulation 2, HG008 tumour/normal mixtures (real tumour CNVs, GRCh38)
  test_noncancer   PRJNA795275 non-cancer plasma (no CNVs expected)
  test_cancer      GSE71378 cancer plasma (weak labels = consensus of the callers)
"""
import os, re, sys
from pathlib import Path
import pandas as pd

E = os.environ
OUT = Path(E["MLO"]) / "manifest.tsv"
COLS = ["dataset", "sample_id", "base_sample", "class", "replicate", "depth_x",
        "tumour_fraction_nominal", "tumour_fraction_realised", "build", "layout", "split", "bam", "truth"]


def num(x):
    m = re.search(r"[0-9]*\.?[0-9]+", str(x))
    return float(m.group(0)) if m else float("nan")


def resolve(p, base):
    p = Path(str(p))
    return p if p.is_absolute() else (base / p)


rows, problems = [], []

# ---------------------------------------------------------------- Simulation 1 (GRCh37)
m = Path(E["SIM1_MANIFEST"])
if m.exists():
    d = pd.read_csv(m, sep="\t", dtype=str).fillna("NA")
    for _, r in d.iterrows():
        bg = r.get("background", r["sample_id"][:4])
        neg = str(r.get("class", "")).lower().startswith("neg")
        rep = "NA" if neg else str(int(num(r.get("replicate", "nan"))))
        if bg == "IH03":
            split = "test_sim1_ih03"
        else:
            split = "validation" if rep == "5" else "train"
        truth = "NA"
        if not neg:
            t = str(r.get("truth", "NA"))
            truth = str(resolve(t, m.parent.parent)) if t not in ("NA", "", "nan") else str(Path(E["SIM1_TRUTH_DIR"]) / f"rep{rep}.truth.tsv")
        rows.append(dict(dataset="sim1", sample_id=r["sample_id"], base_sample=bg, **{"class": "negative" if neg else "spiked"},
                         replicate=rep, depth_x=num(r.get("depth", "nan")), tumour_fraction_nominal=float("nan"),
                         tumour_fraction_realised=float("nan"), build="GRCh37", layout="auto", split=split,
                         bam=str(resolve(r["bam"], m.parent.parent)), truth=truth))
else:
    problems.append(f"missing {m}")

# ---------------------------------------------------------------- Simulation 2, HG008 (GRCh38)
m = Path(E["SIM2_MANIFEST"])
if m.exists():
    d = pd.read_csv(m, sep="\t", dtype=str)
    for _, r in d.iterrows():
        tf = num(r.get("tumor_fraction", r.get("tumour_fraction", "nan")))
        rows.append(dict(dataset="hg008", sample_id=r["sample_id"], base_sample="HG008", **{"class": "mixture"},
                         replicate=str(r.get("replicate", "NA")), depth_x=num(r.get("depth_x", "nan")),
                         tumour_fraction_nominal=tf, tumour_fraction_realised=num(r.get("realised_fraction", "nan")),
                         build="GRCh38", layout="PAIRED", split="test_hg008", bam=r["bam"], truth=E["TRUTH38"]))
else:
    problems.append(f"missing {m}")

# ---------------------------------------------------------------- real plasma
tfs = {}
if Path(E["GSE_TF"]).exists():
    t = pd.read_csv(E["GSE_TF"], sep="\t", dtype=str)
    t = t[(t["method"] == "ichorCNA 1000kb_lowTF") & (t["level"] == "full")]
    tfs = {b: float(v) for b, v in zip(t["base"], t["TF"])}
for ds, key, build, split in (("prjna795275", "PRJNA_MANIFEST", "GRCh38", "test_noncancer"),
                              ("gse71378", "GSE_MANIFEST", "GRCh37", "test_cancer")):
    m = Path(E[key])
    if not m.exists():
        problems.append(f"missing {m}")
        continue
    d = pd.read_csv(m, sep="\t", dtype=str)
    for _, r in d.iterrows():
        base = r["base_sample"]
        cls = "noncancer" if ds == "prjna795275" else ("healthy" if base.startswith("IH") else "cancer")
        rows.append(dict(dataset=ds, sample_id=r["sample_id"], base_sample=base, **{"class": cls}, replicate="NA",
                         depth_x=num(r.get("depth_x", "nan")), tumour_fraction_nominal=float("nan"),
                         tumour_fraction_realised=tfs.get(base, float("nan")) if ds == "gse71378" else float("nan"),
                         build=build, layout=r.get("layout", "auto"), split=split,
                         bam=str(resolve(r["bam"], m.parent)), truth="NA"))

man = pd.DataFrame(rows, columns=COLS)
ok = [Path(b).exists() and (Path(b + ".bai").exists() or Path(b[:-4] + ".bai").exists()) for b in man.bam]
miss = man[[not x for x in ok]]
man = man[ok].drop_duplicates("bam").reset_index(drop=True)
OUT.parent.mkdir(parents=True, exist_ok=True)
man.to_csv(OUT, sep="\t", index=False)
if len(miss):
    miss.to_csv(OUT.with_name("manifest_missing_bams.tsv"), sep="\t", index=False)

print(f"manifest: {OUT}   ({len(man)} BAMs)")
print(man.groupby(["split", "dataset"]).size().to_string())
if len(miss):
    print(f"\nWARNING {len(miss)} BAMs not found or not indexed -> {OUT.with_name('manifest_missing_bams.tsv')}")
for p in problems:
    print("WARNING", p)
if len(man) == 0:
    sys.exit("No BAMs found - nothing to do.")
