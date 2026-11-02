#!/usr/bin/env python3
"""
cfCNV Simulation 1 - FINAL tables: merges the corrected scoring (scorer v2) with the
re-runs and configuration checks.

Main configuration per tool (used in Figures 1-6, Table 1):
  QDNAseq, BIC-seq2, dudeML          as before (scorer v2)
  ichorCNA   pure-sample settings    (check K1)      default run -> "ichorCNA_default"
  CNVkit     pooled normals (LOBO)   (check K2)      flat run    -> "CNVkit_flat"
  CNV-MEANN  v2, depth-normalised    (rescue B)
  Spectre    germline, coverage at the benchmark bin sizes (rescue A)
             original 1-kb --cancer run -> "Spectre_1kb"
  LYCEUM     exon calls merged into events, argmax (rescue C)

  python 06_merge_final.py --v2 ../sim1_v2/tables --inputs inputs --out tables
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

MAP = {"truth_events": "n_truth", "called_events": "n_calls", "grid10k_MCC": "bin_MCC",
       "grid10k_balanced_accuracy": "bin_balanced_accuracy", "grid10k_specificity": "bin_specificity",
       "grid10k_recall": "bin_sensitivity", "grid10k_precision": "bin_precision", "grid10k_f1": "bin_F1",
       "bp_jaccard": "jaccard", "median_breakpoint_error_bp": "breakpoint_error_bp",
       "detected_10pct": "detected_10", "detected_30pct": "detected_30", "sample": "sample_id",
       "depth_x": "depth", "sample_class": "class"}


def conv(d, tool):
    d = d.rename(columns=MAP).copy(); d["tool"] = tool; d["level"] = "segment"
    pos = d["class"] == "positive"
    d.loc[pos & d.F1.isna() & d.status.str.startswith("OK"), "F1"] = 0.0
    d.loc[pos & d.recall.isna() & d.status.str.startswith("OK"), "recall"] = 0.0
    d["recall_10"] = d.detected_10 / d.n_truth; d["recall_30"] = d.detected_30 / d.n_truth
    d["matched_depth"] = d.depth <= 1
    return d


def tconv(t, tool):
    t = t.rename(columns={"sample": "sample_id", "depth_x": "depth", "det50": "detected"}).copy()
    t["tool"] = tool; t["type"] = t["type"].str.upper(); return t


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--v2", required=True); ap.add_argument("--inputs", required=True)
    ap.add_argument("--out", required=True); a = ap.parse_args()
    V, I, O = Path(a.v2), Path(a.inputs), Path(a.out); O.mkdir(parents=True, exist_ok=True)
    pc = pd.read_csv(V / "per_config.tsv", sep="\t"); td = pd.read_csv(V / "truth_detection.tsv", sep="\t")
    rename = {"CNVkit": "CNVkit_flat", "ichorCNA": "ichorCNA_default", "Spectre": "Spectre_1kb"}
    pc = pc[~pc.tool.isin(["CNV-MEANN", "LYCEUM"])].replace({"tool": rename})
    td = td[~td.tool.isin(["CNV-MEANN", "LYCEUM"])].replace({"tool": rename})
    add, addt = [], []
    ck = pd.read_csv(I / "checks_per_config.csv"); ckt = pd.read_csv(I / "checks_truth_detection.csv")
    for tool, src in (("ichorCNA", "ichorCNA"), ("CNVkit", "CNVkit")):
        add.append(conv(ck[ck.tool == src], tool)); addt.append(tconv(ckt[ckt.tool == src], tool))
    for name, tool, cfg in (("CNV-MEANN", "CNV-MEANN", "v2_normalised"), ("Spectre", "Spectre", "germline"),
                            ("LYCEUM", "LYCEUM", "argmax")):
        d = pd.read_csv(I / f"{name}_per_config.csv"); t = pd.read_csv(I / f"{name}_truth_detection.csv")
        add.append(conv(d[d.config == cfg], tool)); addt.append(tconv(t[t.config == cfg], tool))
    cols = list(pc.columns)
    pc = pd.concat([pc] + [x.reindex(columns=cols) for x in add], ignore_index=True)
    td = pd.concat([td] + [x.reindex(columns=list(td.columns)) for x in addt], ignore_index=True)
    pc.to_csv(O / "per_config.tsv", sep="\t", index=False); td.to_csv(O / "truth_detection.tsv", sep="\t", index=False)
    for f in ("measured_depth.tsv", "runtime.tsv"):
        x = pd.read_csv(V / f, sep="\t")
        if "tool" in x: x = x.replace({"tool": {"Spectre": "Spectre"}})
        x.to_csv(O / f, sep="\t", index=False)
    pd.DataFrame([
        {"tool": "ichorCNA", "note": "main = pure-sample settings (normal 0.001-0.1, ploidy 2, no subclonal states); default tumour settings reported as ichorCNA_default"},
        {"tool": "CNVkit", "note": "main = pooled-normal reference from the other two backgrounds' negatives; flat reference reported as CNVkit_flat"},
        {"tool": "CNV-MEANN", "note": "re-implemented network (4-15-3), depth/GC-normalised features, leave-one-replicate-out"},
        {"tool": "Spectre", "note": "germline mode, mosdepth coverage at the benchmark bin size, ENCODE blacklist; original 1-kb run = Spectre_1kb"},
        {"tool": "LYCEUM", "note": "exon calls (argmax) merged into events (gap <= 2 exons, >= 2 exons); bin_bp 0 = exon level"},
    ]).to_csv(O / "configurations.tsv", sep="\t", index=False)
    print(pc.groupby(["tool", "class"]).size().unstack(fill_value=0).to_string())


if __name__ == "__main__":
    main()
