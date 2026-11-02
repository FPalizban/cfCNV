#!/usr/bin/env python3
"""
cfCNV Simulation 1 - build unified analysis tables from the metrics bundle.

  python 01_build_tables.py --data <unzipped cfCNV_Sim1_METRICS_ONLY> --out tables/

Writes
  per_config.tsv        one row per tool x bin x sample (event, bin-level, FP metrics)
  truth_detection.tsv   one row per tool x bin x sample x truth event (detected at RO 0.5)
  measured_depth.tsv    mosdepth autosomal mean depth per simulated BAM
  runtime.tsv           wall-clock seconds per run (Slurm accounting)
  exclusions.tsv        configurations removed from analysis, with the reason

Event metrics follow the bundle's scorers (one-to-one, same chromosome and
direction, >= 50% reciprocal overlap). Precision is NA when a tool made no
calls; F1 is then set to 0 for positive samples (recall is 0).
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

AUTO_BP = 2_881_033_286  # autosomal genome size used by the bundle's scorer


def f1(p, r):
    return np.where((p + r) > 0, 2 * p * r / (p + r), 0.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    D = Path(a.data) / "09_metrics" / "simulation1"
    O = Path(a.out); O.mkdir(parents=True, exist_ok=True)
    excl, rows = [], []

    # ---------------- QDNAseq / CNVkit / ichorCNA / Spectre (master) ----------------
    m = pd.read_csv(D / "simulation1_overlap_master.csv")
    rs = pd.read_csv(D / "run_status.csv")
    rs["out_dir_bin"] = rs.output_file.str.extract(r"/(\d+)/[^/]+$")[0].astype(float)
    m = m.merge(rs[["tool", "sample", "analysis_bin_bp", "out_dir_bin"]],
                on=["tool", "sample", "analysis_bin_bp"], how="left")
    bad = m.out_dir_bin.notna() & (m.out_dir_bin != m.analysis_bin_bp)
    for (t, b), g in m[bad].groupby(["tool", "analysis_bin_bp"]):
        excl.append({"tool": t, "bin_bp": b, "n_configurations": len(g),
                     "reason": f"scored from the {int(g.out_dir_bin.iloc[0])}-bp output file "
                               f"(path error in run_status); results duplicate that bin"})
    m = m[~bad]
    sp = m[m.tool == "Spectre"]
    if len(sp) and (sp.called_events.fillna(0) == 0).all():
        excl.append({"tool": "Spectre", "bin_bp": 1000, "n_configurations": len(sp),
                     "reason": "zero calls in every sample incl. 20x (OK_ZERO_CALLS); "
                               "kept in tables, reported as non-functional on this input"})
    mm = pd.DataFrame({
        "tool": m.tool, "bin_bp": m.analysis_bin_bp.astype(int), "sample_id": m["sample"],
        "background": m.background, "depth": m.depth_x, "replicate": m.replicate,
        "class": m.sample_class, "n_truth": m.truth_events, "n_calls": m.called_events,
        "TP": m.TP_events_50pct, "FP": m.FP_events_50pct, "FN": m.FN_events_50pct,
        "precision": m.precision_50pct, "recall": m.recall_50pct,
        "bin_MCC": m.grid10k_MCC, "bin_balanced_accuracy": m.grid10k_balanced_accuracy,
        "bin_specificity": m.grid10k_specificity, "bin_sensitivity": m.grid10k_recall,
        "jaccard": m.bp_jaccard, "breakpoint_error_bp": m.median_breakpoint_error_bp,
        "false_calls_per_Gb": m.false_calls_per_Gb, "false_altered_fraction": m.false_altered_fraction,
        "detected_10": m.overlap_events_10pct, "detected_30": m.overlap_events_30pct,
        "level": "segment"})
    rows.append(mm)

    # ---------------- BIC-seq2 ----------------
    b = pd.read_csv(D / "BICseq2" / "BICseq2_configuration_metrics.tsv", sep="\t")
    bn = pd.read_csv(D / "BICseq2" / "BICseq2_negative_controls.tsv", sep="\t")
    bb = pd.DataFrame({
        "tool": "BICseq2", "bin_bp": b.bin_bp, "sample_id": b.sample_id, "background": b.background,
        "depth": b.target_depth, "replicate": b.replicate, "class": b["class"],
        "n_truth": b.truth_events, "n_calls": b.called_events, "TP": b.TP_50, "FP": b.FP_50, "FN": b.FN_50,
        "precision": b.precision_50, "recall": b.recall_50,
        "breakpoint_error_bp": b.mean_breakpoint_error_bp_50,
        "detected_10": b.TP_10, "detected_30": b.TP_30, "level": "segment"})
    bn = bn.rename(columns={"sample_id": "sample_id"})
    bn["false_calls_per_Gb"] = bn.false_calls / (AUTO_BP / 1e9)
    bn["false_altered_fraction"] = bn.false_altered_bp / AUTO_BP
    bb = bb.merge(bn[["sample_id", "bin_bp", "false_calls_per_Gb", "false_altered_fraction"]],
                  on=["sample_id", "bin_bp"], how="left")
    pos = bb["class"] == "positive"
    bb.loc[pos, "false_calls_per_Gb"] = bb.loc[pos, "FP"] / (AUTO_BP / 1e9)
    rows.append(bb)

    # ---------------- dudeML ----------------
    dm = pd.read_csv(D / "dudeML" / "plot_ready" / "dudeml_threshold_sensitivity.csv")
    d5 = dm[dm.overlap_threshold == 0.5].copy()
    det = {t: dm[dm.overlap_threshold == t].set_index(["sample_id", "bin_bp"]).TP for t in (0.1, 0.3)}
    dd = pd.DataFrame({
        "tool": "dudeML", "bin_bp": d5.bin_bp, "sample_id": d5.sample_id, "background": d5.background,
        "depth": d5.depth, "replicate": d5.replicate, "class": "positive",
        "n_truth": d5.TP + d5.FN, "n_calls": d5.n_calls, "TP": d5.TP, "FP": d5.FP, "FN": d5.FN,
        "precision": d5.precision, "recall": d5.recall, "bin_MCC": d5.bin_MCC,
        "bin_balanced_accuracy": d5.bin_balanced_accuracy, "bin_specificity": d5.bin_specificity,
        "jaccard": d5.jaccard,
        "breakpoint_error_bp": (d5.median_start_bp.astype(float) + d5.median_end_bp.astype(float)) / 2,
        "false_calls_per_Gb": d5.FP / (AUTO_BP / 1e9), "level": "segment (merged windows)"})
    idx = list(zip(d5.sample_id, d5.bin_bp))
    dd["detected_10"] = [det[0.1].get(k, np.nan) for k in idx]
    dd["detected_30"] = [det[0.3].get(k, np.nan) for k in idx]
    dn = pd.read_csv(D / "dudeML" / "plot_ready" / "dudeml_negative_fp_metrics.csv")
    dn2 = pd.DataFrame({
        "tool": "dudeML", "bin_bp": dn.bin_bp, "sample_id": dn.sample_id, "background": dn.background,
        "depth": dn.depth, "replicate": 0, "class": "negative", "n_truth": 0,
        "n_calls": dn.n_apparent_fp_events, "FP": dn.n_apparent_fp_events,
        "bin_specificity": dn.apparent_specificity, "false_calls_per_Gb": dn.apparent_fp_events_per_gb,
        "false_altered_fraction": dn.apparent_fp_fraction, "level": "segment (merged windows)"})
    rows += [dd, dn2]

    # ---------------- CNV-MEANN (bin-level only) ----------------
    for f in sorted((D / "CNV-MEANN" / "plot_ready").glob("meann_bin_*_metrics.csv")):
        if "pilot" in f.name:
            continue
        c = pd.read_csv(f)
        rows.append(pd.DataFrame({
            "tool": "CNV-MEANN", "bin_bp": c.bin_bp, "sample_id": c.sample_id, "background": c.background,
            "depth": c.depth, "replicate": c.replicate, "class": c["class"],
            "bin_MCC": c.MCC_binary, "bin_specificity": c.specificity_binary,
            "bin_sensitivity": c.recall_binary, "bin_precision": c.precision_binary,
            "bin_F1": c.F1_binary, "predicted_cnv_fraction": c.predicted_cnv_fraction,
            "n_calls": c.total_event_runs, "level": "bin only"}))
    excl.append({"tool": "CNV-MEANN", "bin_bp": "all", "n_configurations": 684,
                 "reason": "bin-level classifier output only (no event intervals); class-weighted "
                           "model (prior 20:1:1) predicts ~0 CNV bins - event metrics not computable"})

    # ---------------- LYCEUM (exon level, build mismatch) ----------------
    ly = pd.read_csv(D / "LYCEUM" / "adapted_GRCh37" / "LYCEUM_per_sample_exon_metrics.csv")
    rows.append(pd.DataFrame({
        "tool": "LYCEUM", "bin_bp": 0, "sample_id": ly.sample_id, "background": ly.background,
        "depth": ly.depth, "replicate": ly.replicate, "class": ly["class"],
        "bin_precision": ly.binary_precision, "bin_sensitivity": ly.binary_recall, "bin_F1": ly.binary_F1,
        "n_calls": ly.n_called_exons, "level": "exon (hg38 model on GRCh37 truth)"}))
    excl.append({"tool": "LYCEUM", "bin_bp": "exon", "n_configurations": 114,
                 "reason": "released model/targets are hg38, Simulation 1 is GRCh37; scope table marks "
                           "'excluded_incompatible'; adapted exon scores reported descriptively only"})

    pc = pd.concat(rows, ignore_index=True)
    pc["F1"] = np.where(pc["class"] == "positive",
                        np.where(pc.precision.notna(), f1(pc.precision.fillna(0), pc.recall.fillna(0)),
                                 np.where(pc.recall.notna(), 0.0, np.nan)), np.nan)
    pc.loc[pc.level.isin(["bin only"]) | pc.level.str.startswith("exon"), "F1"] = np.nan
    for t in (10, 30):
        pc[f"recall_{t}"] = pc[f"detected_{t}"] / pc.n_truth.replace(0, np.nan)
    pc["matched_depth"] = pc.depth <= 1
    pc.to_csv(O / "per_config.tsv", sep="\t", index=False)

    # ---------------- truth detection ----------------
    td = []
    t = pd.read_csv(D / "truth_event_detection.csv")
    t = t.merge(m[["tool", "sample", "analysis_bin_bp"]].drop_duplicates(),
                on=["tool", "sample", "analysis_bin_bp"])  # drops excluded bins
    td.append(pd.DataFrame({"tool": t.tool, "bin_bp": t.analysis_bin_bp, "sample_id": t["sample"],
                            "background": t.background, "depth": t.depth_x, "replicate": t.replicate,
                            "size_bp": t.truth_size_bp, "type": t.truth_type.str.upper(),
                            "detected": t.detected_50pct}))
    t = pd.read_csv(D / "BICseq2" / "BICseq2_truth_event_detection.tsv", sep="\t")
    td.append(pd.DataFrame({"tool": "BICseq2", "bin_bp": t.bin_bp, "sample_id": t.sample_id,
                            "background": t.background, "depth": t.target_depth, "replicate": t.replicate,
                            "size_bp": t.cnv_size_bp, "type": t.cnv_type.str.upper(), "detected": t.detected_50}))
    t = pd.read_csv(D / "dudeML" / "truth_detection.tsv", sep="\t")
    td.append(pd.DataFrame({"tool": "dudeML", "bin_bp": t.bin_bp, "sample_id": t.sample_id,
                            "background": t.background, "depth": t.depth, "replicate": t.replicate,
                            "size_bp": t.cnv_size_bp, "type": t.alteration.str.upper(), "detected": t.detected}))
    pd.concat(td, ignore_index=True).to_csv(O / "truth_detection.tsv", sep="\t", index=False)

    # ---------------- measured depth (mosdepth) ----------------
    md = []
    for f in (Path(a.data) / "07_tool_outputs" / "simulation1" / "Spectre").glob("*/mosdepth/*.summary.txt"):
        s = pd.read_csv(f, sep="\t")
        s = s[s.chrom.astype(str).isin([str(i) for i in range(1, 23)])]
        md.append({"sample_id": f.parent.parent.name,
                   "measured_depth": (s.bases.sum() / s.length.sum())})
    if not md:   # v2 bundle: use samtools-coverage QC tables (03_qc/simulation1_depths)
        for f in (Path(a.data) / "03_qc" / "simulation1_depths").glob("*.coverage.tsv"):
            s = pd.read_csv(f, sep="\t")
            s = s[s["#rname"].astype(str).isin([str(i) for i in range(1, 23)])]
            L = s.endpos - s.startpos + 1
            md.append({"sample_id": f.name.replace(".coverage.tsv", ""),
                       "measured_depth": (s.meandepth * L).sum() / L.sum()})
    md = pd.DataFrame(md)
    md.to_csv(O / "measured_depth.tsv", sep="\t", index=False)

    # ---------------- runtime ----------------
    rt = []
    s = pd.read_csv(D / "slurm_resource_usage.csv")
    s = s[s.component != "reference"]
    st = s.State.where(~((s.component == "Spectre") & (s.ExitCode == "5:0")),
                       "COMPLETED (Spectre finished; wrapper exit 5)")
    rt.append(pd.DataFrame({"tool": s.component, "wall_s": s.ElapsedRaw, "state": st}))
    s = pd.read_csv(D / "BICseq2" / "BICseq2_slurm_resources.tsv", sep="|")
    s = s[s.JobName.str.startswith("bic_normseg") & ~s.JobID.astype(str).str.contains(r"\.")]
    rt.append(pd.DataFrame({"tool": "BICseq2", "wall_s": s.ElapsedRaw, "state": s.State}))
    s = pd.read_csv(D / "CNV-MEANN" / "plot_ready" / "CNV_MEANN_runtime_by_bin.csv")
    rt.append(pd.DataFrame({"tool": "CNV-MEANN", "wall_s": s.batch_walltime_seconds / s.n_configurations,
                            "state": "COMPLETED (batch / n)"}))
    pd.concat(rt, ignore_index=True).to_csv(O / "runtime.tsv", sep="\t", index=False)

    pd.DataFrame(excl).to_csv(O / "exclusions.tsv", sep="\t", index=False)
    print(pd.DataFrame(excl).to_string(index=False))
    print("\nconfigurations per tool (after exclusions):")
    print(pc.groupby(["tool", "class"]).size().unstack(fill_value=0).to_string())


if __name__ == "__main__":
    main()
