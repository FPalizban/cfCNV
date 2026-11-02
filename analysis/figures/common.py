"""Data loading and harmonised names shared by the figure scripts (no plotting style here)."""
import numpy as np, pandas as pd
from pathlib import Path
D = Path(__file__).resolve().parent / "data"
NAME = {"BICseq2": "BIC-seq2"}
nm = lambda t: NAME.get(t, t)

def sim1():
    pc = pd.read_csv(D / "sim1_per_config.tsv", sep="\t"); pc["tool"] = pc.tool.map(nm)
    td = pd.read_csv(D / "sim1_truth_detection.tsv", sep="\t"); td["tool"] = td.tool.map(nm)
    prim = pd.read_csv(D / "sim1_primary_configs.tsv", sep="\t"); prim["tool"] = prim.tool.map(nm)
    return pc, td, prim.set_index("tool").bin_bp.to_dict()

def sim2():
    pr = pd.read_csv(D / "sim2_profile_correlation.csv"); pr["tool"] = pr.tool.map(nm)
    tf = pd.read_csv(D / "sim2_tf_estimates.csv")
    mix = pd.read_csv(D / "sim2_mixtures_realised.tsv", sep="\t")
    prim = pd.read_csv(D / "sim2_primary_configs.tsv", sep="\t"); prim["tool"] = prim.tool.map(nm)
    return pr, tf, mix, prim.set_index("tool")

def lod_strict(pr):
    """Strict limit of detection: null = pure-normal mixtures; threshold = max(null max, null mean + 3 SD);
    LoD = lowest tumour fraction at which every replicate exceeds the threshold."""
    rows = []
    for (t, c, b), g in pr.groupby(["tool", "config", "bin_bp"]):
        null = g[g.nominal_fraction == 0].pearson_r; thr = max(null.max(), null.mean() + 3 * null.std())
        for d, h in g[g.nominal_fraction > 0].groupby("depth_x"):
            mn = h.groupby("nominal_fraction").pearson_r.min(); ok = mn[mn > thr]
            rows.append({"tool": t, "config": c, "bin_bp": b, "depth_x": d, "threshold": thr,
                         "LoD": float(ok.index.min()) if len(ok) else np.nan})
    return pd.DataFrame(rows)

def primary(df, prim):
    """keep rows of each tool's primary Simulation 2 configuration"""
    k = df.merge(prim[["config", "bin_bp"]].reset_index(), on=["tool", "config", "bin_bp"])
    return k

def healthy():
    b = pd.read_csv(D / "healthy_burden.tsv", sep="\t")
    tf = pd.read_csv(D / "healthy_tumour_fraction.tsv", sep="\t")
    ins = pd.read_csv(D / "healthy_insert_sizes.tsv", sep="\t")
    calls = pd.read_csv(D / "healthy_calls.tsv", sep="\t")
    rec = pd.read_csv(D / "healthy_recurrent_windows.tsv", sep="\t")
    ann = pd.read_csv(D / "healthy_annotation_vs_random.tsv", sep="\t")
    return b, tf, ins, calls, rec, ann

HEALTHY_PRIMARY = [("QDNAseq", "100kb"), ("Spectre", "100kb"), ("BIC-seq2", "10kb"), ("CNVkit", "30kb"), ("ichorCNA", "1Mb low-TF")]
