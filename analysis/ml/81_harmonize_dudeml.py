#!/usr/bin/env python3
"""
cfCNV - harmonize dudeML window-level predictions into CNV intervals.

dudeML `predict` writes one row per window. This script
  1. reads each prediction file,
  2. maps the predicted class to GAIN / LOSS / NEUTRAL,
  3. merges adjacent windows with the same state into segments,
  4. writes the same harmonized layout used for BIC-seq2
     (one row per GAIN/LOSS segment, plus a status table).

The exact column layout of `predict` output must be confirmed first:

    python 81_harmonize_dudeml.py --inspect /path/to/one_prediction_file

then pass the column indices (0-based) with --chrom-col/--start-col/
--end-col/--class-col. Nothing is guessed silently.

Example:
    python 81_harmonize_dudeml.py \
        --tasks  $SIM1_MAN/bicseq2_tasks.tsv \
        --pattern "$SIM1_OUT/dudeML/{sample}/w{bin}/{sample}.w{bin}.pred.bed" \
        --chrom-col 0 --start-col 1 --end-col 2 --class-col 3 \
        --outdir $PROJECT/08_harmonized_calls/simulation1/dudeML

`--tasks` can be any TSV with columns sample_id, background, target_depth,
replicate, class, bin_bp (the BIC-seq2 task manifest has all of them,
and dudeML used the same 114 samples x 6 window sizes = 684 configs).
"""

import argparse
import csv
import sys
from pathlib import Path


def inspect(path):
    print(f"FILE: {path}\n")
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i >= 8:
                break
            fields = line.rstrip("\n").split("\t")
            print(f"line {i + 1}: {len(fields)} tab-separated fields")
            for j, f in enumerate(fields[:12]):
                print(f"   [{j}] {f}")
            if len(fields) > 12:
                print(f"   ... {len(fields) - 12} more")
    print("\nPick the column indices for chrom/start/end/predicted class "
          "and re-run with --chrom-col --start-col --end-col --class-col.")


def to_state(label):
    s = str(label).strip().lower()
    if s in ("del", "deletion", "loss") or s.startswith("del"):
        return "LOSS"
    if s in ("dup", "duplication", "gain") or s.startswith("dup"):
        return "GAIN"
    if s in ("n", "normal", "neutral", "0", "none"):
        return "NEUTRAL"
    return None


def merge_windows(rows, max_gap):
    """rows: list of (chrom, start, end, state) sorted by chrom/start."""
    segs = []
    cur = None
    for chrom, start, end, state in rows:
        if state == "NEUTRAL":
            if cur:
                segs.append(cur)
                cur = None
            continue
        if (cur and cur["chrom"] == chrom and cur["state"] == state
                and start - cur["end"] <= max_gap):
            cur["end"] = max(cur["end"], end)
            cur["n_windows"] += 1
        else:
            if cur:
                segs.append(cur)
            cur = {"chrom": chrom, "start": start, "end": end,
                   "state": state, "n_windows": 1}
    if cur:
        segs.append(cur)
    return segs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inspect")
    ap.add_argument("--tasks")
    ap.add_argument("--pattern",
                    help="prediction path with {sample} and {bin} placeholders")
    ap.add_argument("--chrom-col", type=int)
    ap.add_argument("--start-col", type=int)
    ap.add_argument("--end-col", type=int)
    ap.add_argument("--class-col", type=int)
    ap.add_argument("--max-gap-windows", type=float, default=1.0,
                    help="merge same-state windows separated by <= this many "
                         "window lengths (default 1 = only touching/adjacent)")
    ap.add_argument("--min-windows", type=int, default=1,
                    help="drop segments supported by fewer windows (default 1 = keep all)")
    ap.add_argument("--outdir")
    a = ap.parse_args()

    if a.inspect:
        inspect(a.inspect)
        return

    need = [a.tasks, a.pattern, a.outdir, a.chrom_col, a.start_col,
            a.end_col, a.class_col]
    if any(x is None for x in need):
        sys.exit("ERROR: run --inspect first, then give --tasks --pattern "
                 "--outdir and all four column indices.")

    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(a.tasks) as fh:
        tasks = list(csv.DictReader(fh, delimiter="\t"))

    calls, status = [], []
    unknown_labels = set()

    for t in tasks:
        sample, bin_bp = t["sample_id"], t["bin_bp"]
        pred = Path(a.pattern.format(sample=sample, bin=bin_bp))
        base = {"task_id": t.get("task_id", ""), "sample_id": sample,
                "bin_bp": bin_bp}

        if not pred.exists():
            status.append({**base, "status": "MISSING", "message": str(pred)})
            continue

        rows = []
        with pred.open() as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) <= max(a.chrom_col, a.start_col, a.end_col, a.class_col):
                    continue
                try:
                    start = int(float(f[a.start_col]))
                    end = int(float(f[a.end_col]))
                except ValueError:
                    continue  # header line
                state = to_state(f[a.class_col])
                if state is None:
                    unknown_labels.add(f[a.class_col])
                    continue
                rows.append((f[a.chrom_col].replace("chr", ""), start, end, state))

        rows.sort(key=lambda r: (r[0], r[1]))
        win = int(bin_bp)
        segs = merge_windows(rows, max_gap=a.max_gap_windows * win)
        segs = [s for s in segs if s["n_windows"] >= a.min_windows]

        for s in segs:
            calls.append({
                "tool": "dudeML",
                "task_id": t.get("task_id", ""),
                "sample_id": sample,
                "background": t.get("background", ""),
                "target_depth": t.get("target_depth", ""),
                "replicate": t.get("replicate", ""),
                "class": t.get("class", ""),
                "bin_bp": bin_bp,
                "chrom": s["chrom"],
                "start_0based": s["start"],
                "end": s["end"],
                "state": s["state"],
                "n_windows": s["n_windows"],
                "source_file": str(pred),
            })
        status.append({**base, "status": "OK",
                       "message": f"windows={len(rows)};segments={len(segs)}"})

    if unknown_labels:
        sys.exit("ERROR: unrecognised class labels "
                 f"{sorted(unknown_labels)[:10]} - check --class-col, or add "
                 "the label to to_state().")

    def write(path, rows, fields):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                               lineterminator="\n")
            w.writeheader()
            w.writerows(rows)

    call_fields = ["tool", "task_id", "sample_id", "background", "target_depth",
                   "replicate", "class", "bin_bp", "chrom", "start_0based",
                   "end", "state", "n_windows", "source_file"]
    write(outdir / "dudeML_calls.tsv", calls, call_fields)
    write(outdir / "dudeML_harmonization_status.tsv", status,
          ["task_id", "sample_id", "bin_bp", "status", "message"])

    n_ok = sum(s["status"] == "OK" for s in status)
    print(f"Configurations: {len(tasks)}  OK: {n_ok}  "
          f"MISSING: {len(tasks) - n_ok}")
    print(f"GAIN/LOSS segments: {len(calls)}")
    print(f"Output: {outdir}")


if __name__ == "__main__":
    main()
