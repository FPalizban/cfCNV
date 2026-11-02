#!/usr/bin/env python3
"""Harmonise the native output of every caller into one table and summarise call burden and
tumour-fraction estimates (same parsers and thresholds as the cfCNV paper).

usage: collect_calls.py --manifest manifest.tsv --fai ref.fa.fai --outdir . [--bic-threshold 0.2]
manifest columns (tab): sample_id base level group expected_tf tool config kind file
   kind = calls | tf
outputs: all_calls.tsv  call_burden.tsv  tumour_fraction.tsv
all_calls.tsv columns: sample_id base level group tool config chrom start end type   (chrom without 'chr')
"""
import argparse, gzip, re
from pathlib import Path
import numpy as np, pandas as pd

AUTO = {str(i) for i in range(1, 23)}
nc = lambda c: str(c)[3:] if str(c).startswith("chr") else str(c)


def qdnaseq(f):
    out = []
    for ln in open(f):
        if ln.startswith("#") or ln.count("\t") < 7:
            continue
        p = ln.split("\t"); info = dict(x.split("=", 1) for x in p[7].split(";") if "=" in x)
        if nc(p[0]) in AUTO and "END" in info:
            out.append((nc(p[0]), int(p[1]), int(info["END"]), "gain" if info.get("SVTYPE") == "DUP" else "loss"))
    return out


def cnvkit(f):
    d = pd.read_csv(f, sep="\t")
    return [(nc(c), int(s), int(e), "gain" if n > 2 else "loss")
            for c, s, e, n in zip(d.chromosome, d.start, d.end, d.cn) if nc(c) in AUTO and n != 2]


def ichor(f):
    d = pd.read_csv(f, sep="\t")
    col = "Corrected_Copy_Number" if "Corrected_Copy_Number" in d.columns else "copy.number"
    return [(nc(c), int(s), int(e), "gain" if n > 2 else "loss")
            for c, s, e, n in zip(d.chrom, d.start, d.end, d[col]) if nc(c) in AUTO and pd.notna(n) and n != 2]


def spectre(f):
    out = []
    op = gzip.open if str(f).endswith(".gz") else open
    with op(f, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            p = ln.rstrip().split("\t")
            if len(p) >= 4 and nc(p[0]) in AUTO:
                out.append((nc(p[0]), int(p[1]), int(p[2]), "gain" if p[3].upper().startswith("DUP") else "loss"))
    return out


def bicseq2(f, thr):
    d = pd.read_csv(f, sep="\t")
    lc = [c for c in d.columns if "log2" in c.lower()][0]
    return [(nc(c), int(s), int(e), "gain" if v >= thr else "loss")
            for c, s, e, v in zip(d.iloc[:, 0], d.iloc[:, 1], d.iloc[:, 2], d[lc]) if nc(c) in AUTO and abs(v) >= thr]


def dudeml(f):
    """dudeML per-window predictions (chr, start 1-based, end, ID, coverage, CNV, ...): adjacent windows of the
    same class are merged, as in the paper's Simulation 1 scoring"""
    op = gzip.open if str(f).endswith(".gz") else open
    rows = []
    with op(f, "rt") as fh:
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 6 or nc(p[0]) not in AUTO: continue
            k = p[5].strip().lower()
            t = "loss" if k.startswith("del") else "gain" if k.startswith("dup") else None
            rows.append((nc(p[0]), int(float(p[1])) - 1, int(float(p[2])), t))
    rows.sort(key=lambda r: (int(r[0]), r[1])); out, cur = [], None
    for c, s, e, t in rows:
        if cur and t and cur[0] == c and cur[3] == t and cur[2] == s: cur[2] = e; continue
        if cur: out.append(tuple(cur)); cur = None
        if t: cur = [c, s, e, t]
    if cur: out.append(tuple(cur))
    return out


def generic(f):
    """already harmonised: chrom start end type (header optional)"""
    out = []
    for ln in open(f):
        p = ln.rstrip().split("\t")
        if len(p) < 4 or not p[1].isdigit():
            continue
        t = p[3].lower()
        t = "loss" if t.startswith(("loss", "del")) else "gain" if t.startswith(("gain", "dup", "amp")) else None
        if t and nc(p[0]) in AUTO:
            out.append((nc(p[0]), int(p[1]), int(p[2]), t))
    return out


def tf_rows(tool, f):
    f = Path(f)
    if tool == "ichorCNA":
        t = f.read_text()
        m = re.search(r"Tumor Fraction:\s*([0-9.eE-]+)", t); p = re.search(r"Ploidy:\s*([0-9.eE-]+)", t)
        return [("ichorCNA", float(m.group(1)), float(p.group(1)) if p else np.nan)]
    if tool == "ACE":
        a = pd.read_csv(f, sep="\t").iloc[0]
        return [("ACE ploidy 2", a.ace_cellularity_ploidy2, 2.0),
                ("ACE free ploidy", a.ace_cellularity_free, a.ace_ploidy_free)]
    if tool == "Fragle":
        d = pd.read_csv(f)
        return [("Fragle", float(d.iloc[0, 1]), np.nan)]
    return []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True); ap.add_argument("--fai", required=True)
    ap.add_argument("--outdir", default="."); ap.add_argument("--bic-threshold", type=float, default=0.2)
    a = ap.parse_args(); O = Path(a.outdir); O.mkdir(parents=True, exist_ok=True)
    genome = sum(int(l.split()[1]) for l in open(a.fai) if nc(l.split()[0]) in AUTO)
    M = pd.read_csv(a.manifest, sep="\t", dtype=str).fillna("")
    parsers = {"QDNAseq": qdnaseq, "CNVkit": cnvkit, "ichorCNA": ichor, "Spectre": spectre,
               "BIC-seq2": lambda f: bicseq2(f, a.bic_threshold), "dudeML": dudeml}
    calls, burden, tf = [], [], []
    for r in M.itertuples():
        key = dict(sample_id=r.sample_id, base=r.base, level=r.level, group=r.group)
        if r.kind == "tf":
            for meth, v, pl in tf_rows(r.tool, r.file):
                tf.append({**key, "expected_tf": r.expected_tf, "method": f"{meth} {r.config}".strip(), "TF": v, "ploidy": pl})
            continue
        c = [x for x in parsers.get(r.tool, generic)(r.file) if x[2] > x[1]]
        calls += [{**key, "tool": r.tool, "config": r.config, "chrom": ch, "start": s, "end": e, "type": t} for ch, s, e, t in c]
        L = np.array([e - s for _, s, e, _ in c]) if c else np.array([])
        burden.append({**key, "tool": r.tool, "config": r.config, "n_calls": len(c),
                       "n_gain": sum(t == "gain" for *_, t in c), "n_loss": sum(t == "loss" for *_, t in c),
                       "altered_Mb": L.sum() / 1e6, "genome_altered_pct": 100 * L.sum() / genome,
                       "median_call_kb": float(np.median(L) / 1e3) if len(L) else np.nan,
                       "calls_ge_10Mb": int((L >= 1e7).sum()) if len(L) else 0})
    cols = ["sample_id", "base", "level", "group", "tool", "config", "chrom", "start", "end", "type"]
    pd.DataFrame(calls, columns=cols).to_csv(O / "all_calls.tsv", sep="\t", index=False)
    pd.DataFrame(burden).to_csv(O / "call_burden.tsv", sep="\t", index=False)
    pd.DataFrame(tf, columns=["sample_id", "base", "level", "group", "expected_tf", "method", "TF", "ploidy"]).to_csv(
        O / "tumour_fraction.tsv", sep="\t", index=False)
    print(f"{len(calls)} calls from {len(burden)} caller runs; {len(tf)} tumour-fraction estimates")


if __name__ == "__main__":
    main()
