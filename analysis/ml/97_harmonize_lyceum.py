#!/usr/bin/env python3
"""
cfCNV - harmonize LYCEUM exon-level calls (HG008 arm).

LYCEUM exon-level output columns (README): sample, chromosome, exon start,
exon end, prediction. Calls stay at exon resolution (Methods 2.7.2: the
segment vs exon/gene distinction is retained); an optional merged-segment
table is also written for interval-overlap scoring against NIST truth.

Prediction labels are checked, not assumed: any label that is not
recognised stops the script and lists what it found.

Usage:
  python 97_harmonize_lyceum.py --root $PROJECT/07_tool_outputs/hg008/LYCEUM \
      --manifest $PROJECT/00_metadata/manifests/hg008/hg008_mixtures.tsv \
      --outdir  $PROJECT/08_harmonized_calls/hg008/LYCEUM
"""
import argparse
import csv
import sys
from pathlib import Path


def to_state(x):
    s = str(x).strip().upper()
    if s.startswith("DEL") or s in ("LOSS", "1"):
        return "LOSS"
    if s.startswith("DUP") or s in ("GAIN", "2"):
        return "GAIN"
    if s in ("NO-CALL", "NOCALL", "NEUTRAL", "NORMAL", "N", "0", "NO_CALL"):
        return "NEUTRAL"
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--merge-gap", type=int, default=100_000,
                    help="merge same-state exons closer than this (bp) into one segment")
    a = ap.parse_args()

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    meta = {r["sample_id"]: r for r in csv.DictReader(open(a.manifest), delimiter="\t")}

    exon_rows, seg_rows, unknown, status = [], [], set(), []
    for sample, m in meta.items():
        files = sorted((Path(a.root) / sample / "calls").glob("*"))
        if not files:
            status.append({"sample_id": sample, "status": "MISSING"})
            continue
        rows = []
        for f in files:
            for line in open(f):
                p = line.rstrip("\n").split("\t")
                if len(p) < 5:
                    continue
                try:
                    st, en = int(float(p[2])), int(float(p[3]))
                except ValueError:
                    continue
                state = to_state(p[4])
                if state is None:
                    unknown.add(p[4])
                    continue
                rows.append((p[1].replace("chr", ""), st, en, state, p[4]))
        rows.sort(key=lambda r: (r[0], r[1]))
        for chrom, st, en, state, raw in rows:
            if state != "NEUTRAL":
                exon_rows.append({"tool": "LYCEUM", "sample_id": sample,
                                  "tumor_fraction": m.get("tumor_fraction", ""),
                                  "target_depth": m.get("target_depth", ""),
                                  "chrom": chrom, "start_0based": st, "end": en,
                                  "state": state, "raw_label": raw})
        cur = None
        for chrom, st, en, state, _ in rows:
            if state == "NEUTRAL":
                # an exon called neutral ends the current segment
                if cur:
                    seg_rows.append(cur)
                    cur = None
                continue
            if cur and cur["chrom"] == chrom and cur["state"] == state and st - cur["end"] <= a.merge_gap:
                cur["end"], cur["n_exons"] = en, cur["n_exons"] + 1
            else:
                if cur:
                    seg_rows.append(cur)
                cur = {"tool": "LYCEUM", "sample_id": sample,
                       "tumor_fraction": m.get("tumor_fraction", ""),
                       "target_depth": m.get("target_depth", ""),
                       "chrom": chrom, "start_0based": st, "end": en,
                       "state": state, "n_exons": 1}
        if cur:
            seg_rows.append(cur)
        status.append({"sample_id": sample, "status": "OK"})

    if unknown:
        sys.exit(f"ERROR: unrecognised LYCEUM labels: {sorted(unknown)} - "
                 "add them to to_state() after checking the README/toy output.")

    def write(path, rows):
        if not rows:
            Path(path).write_text("")
            return
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                               lineterminator="\n")
            w.writeheader()
            w.writerows(rows)

    write(out / "LYCEUM_exon_calls.tsv", exon_rows)
    write(out / "LYCEUM_merged_segments.tsv", seg_rows)
    write(out / "LYCEUM_status.tsv", status)
    print(f"samples OK {sum(s['status']=='OK' for s in status)}/{len(status)}; "
          f"exon calls {len(exon_rows)}; merged segments {len(seg_rows)}")


if __name__ == "__main__":
    main()
