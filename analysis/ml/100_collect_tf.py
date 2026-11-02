#!/usr/bin/env python3
"""
cfCNV - collect tumour-fraction estimates (ACE, Fragle, optional ichorCNA)
and compare with the known fraction (HG008 mixtures; Sim-1 healthy = 0).

Outputs (in --outdir):
  tf_estimates_long.tsv   one row per sample x tool
  tf_error_summary.tsv    n, MAE, signed bias, RMSE per tool x dataset x true_tf x depth
                          (only rows with a known true_tf)

Usage:
  python 100_collect_tf.py --manifest .../tf_samples.tsv \
     --tfroot $PROJECT/07_tool_outputs/tf --ace-bin 1000 \
     [--ichor-glob "$PROJECT/07_tool_outputs/hg008/ichorCNA/*/*.params.txt"] \
     --outdir $PROJECT/09_metrics/tf
"""
import argparse
import csv
import glob
import math
import re
import sys
from collections import defaultdict
from pathlib import Path


def fnum(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def read_ace(path):
    with open(path) as fh:
        r = next(csv.DictReader(fh, delimiter="\t"))
    return fnum(r["ace_bestfit"])


def read_fragle(path, sample):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        return None
    cols = list(rows[0].keys())
    tf_col = next((c for c in cols if re.search(r"ctdna|fraction|\btf\b", c, re.I)), None)
    if tf_col is None:
        sys.exit(f"ERROR: cannot find the ctDNA-fraction column in {path}; columns: {cols}")
    row = next((r for r in rows if sample in " ".join(map(str, r.values()))), rows[0])
    v = fnum(row[tf_col])
    # Fragle may report percent; harmonise to a 0-1 fraction
    return v / 100 if v is not None and v > 1 else v


def read_ichor(path):
    lines = Path(path).read_text().splitlines()
    for i, l in enumerate(lines):
        m = re.match(r"\s*Tumor Fraction:\s*([0-9.eE+-]+)", l)
        if m:
            return fnum(m.group(1))
        if l.startswith("Sample") and "Tumor Fraction" in l and i + 1 < len(lines):
            h = l.split("\t"); v = lines[i + 1].split("\t")
            return fnum(v[h.index("Tumor Fraction")])
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--tfroot", required=True)
    ap.add_argument("--ace-bin", default="1000")
    ap.add_argument("--ichor-glob")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    man = list(csv.DictReader(open(a.manifest), delimiter="\t"))
    root = Path(a.tfroot)

    ichor = {}
    if a.ichor_glob:
        for p in glob.glob(a.ichor_glob):
            sid = Path(p).name.replace(".params.txt", "")
            ichor[sid] = read_ichor(p)

    long_rows = []
    for m in man:
        s, ds = m["sample_id"], m["dataset"]
        true_tf = fnum(m.get("true_tf"))
        ests = {}
        p = root / "ACE" / f"bin{a.ace_bin}kb" / ds / s / f"{s}.ace.tsv"
        if p.exists():
            ests["ACE"] = read_ace(p)
        p = root / "Fragle" / ds / s / "Fragle.csv"
        if p.exists():
            ests["Fragle"] = read_fragle(p, s)
        if s in ichor:
            ests["ichorCNA"] = ichor[s]
        for tool, est in ests.items():
            long_rows.append({
                "sample_id": s, "dataset": ds, "build": m["build"],
                "target_depth": m.get("target_depth", ""),
                "true_tf": "" if true_tf is None else true_tf,
                "tool": tool, "estimated_tf": "" if est is None else est,
                "error": "" if (est is None or true_tf is None) else est - true_tf,
            })

    with open(out / "tf_estimates_long.tsv", "w", newline="") as fh:
        f = ["sample_id", "dataset", "build", "target_depth", "true_tf",
             "tool", "estimated_tf", "error"]
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(long_rows)

    groups = defaultdict(list)
    for r in long_rows:
        if r["error"] != "":
            groups[(r["tool"], r["dataset"], r["true_tf"], r["target_depth"])].append(r["error"])
    summ = []
    for (tool, ds, tf, dp), errs in sorted(groups.items(), key=lambda x: str(x[0])):
        n = len(errs)
        summ.append({"tool": tool, "dataset": ds, "true_tf": tf, "target_depth": dp,
                     "n": n, "MAE": sum(abs(e) for e in errs) / n,
                     "bias": sum(errs) / n,
                     "RMSE": math.sqrt(sum(e * e for e in errs) / n)})
    with open(out / "tf_error_summary.tsv", "w", newline="") as fh:
        f = ["tool", "dataset", "true_tf", "target_depth", "n", "MAE", "bias", "RMSE"]
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(summ)

    per_tool = defaultdict(int)
    for r in long_rows:
        per_tool[r["tool"]] += 1
    print("estimates per tool:", dict(per_tool))
    print("summary groups:", len(summ))
    print("output:", out)


if __name__ == "__main__":
    main()
