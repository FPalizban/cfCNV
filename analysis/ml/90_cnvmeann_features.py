#!/usr/bin/env python3
"""
cfCNV - CNV-MEANN feature extraction, genome-wide (ADAPTED input layer).

The released CNV-MEANN preprocessing (main_simulation.py) is research code
for chr21 only: hard-coded paths, a per-base read-count matrix for the whole
genome (not feasible for 22 autosomes), and labels assigned by the ROW
POSITION of events in a GroundTruthCNV file. It cannot run on cfCNV data
unchanged. This script re-implements the same five features per bin and
writes them in the SAME 6-column layout main_simulation.py writes
(bin, rd, gc, rel, qua, label; tab-separated; no header), one file per
chromosome, so the released MATLAB training/test code can read them.
The model architecture and decision logic are not touched.

Features per bin (as in main_simulation.py):
  rd    mean per-base read depth in the bin
  gc    GC fraction of the reference bin
  rel   |rd - mean(rd of 5 bins on each side)|  (one side at chromosome ends)
  qua   mean MAPQ of reads starting in the bin, divided by the max over bins
  label 0 normal, 1 gain, 2 hemi_loss, 3 homo_loss

Documented deviations from the released code (record in the Supplement):
  D1  all autosomes processed, one chromosome at a time (original: chr21).
  D2  bin size is a parameter (original: 1,000 bp fixed).
  D3  rel uses true 5-bin neighbour means; the original slices 4 bins and
      divides by 5.
  D4  MAPQ uses all reads starting in the bin; the original keeps only the
      first read at each start position.
  D5  labels come from truth coordinates (GAIN -> 1, LOSS CN1 -> 2,
      CN0 -> 3) instead of event row order.
  D6  zero-depth bins: the original DROPS every bin with rd == 0. At
      0.1-0.5x most 1-kb bins have no reads, so this removes most of the
      genome. Default here keeps them (--drop-zero reproduces the original).
Bins containing any N in the reference are dropped, as in the original.

Usage:
  python 90_cnvmeann_features.py --bam X.bam --ref human_g1k_v37.fasta \
      --sample X --outdir OUT --bin 1000 [--truth rep1.truth.tsv]
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

AUTOSOMES = [str(i) for i in range(1, 23)]


def read_truth(path):
    """Read a truth table robustly; returns DataFrame chrom,start,end,label.

    Handles headers that start with '#', tab/comma/space separators,
    'chr' prefixes, and type columns given as GAIN/LOSS/DUP/DEL or CN.
    """
    raw = Path(path).read_text().splitlines()
    raw = [l for l in raw if l.strip()]
    if not raw:
        sys.exit(f"ERROR: empty truth file {path}")

    sep = "\t" if "\t" in raw[0] else ("," if "," in raw[0] else r"\s+")
    first = re.split(sep, raw[0].lstrip("#").strip())
    has_header = not re.fullmatch(r"\d+", first[1] if len(first) > 1 else "")

    if has_header:
        cols = [re.sub(r"[^a-z0-9]", "", c.lower()) for c in first]
        df = pd.DataFrame([re.split(sep, l.strip()) for l in raw[1:]
                           if not l.startswith("#")])
        df = df.iloc[:, :len(cols)]
        df.columns = cols[:df.shape[1]]
    else:
        df = pd.DataFrame([re.split(sep, l.strip()) for l in raw
                           if not l.startswith("#")])
        cols = [f"c{i}" for i in range(df.shape[1])]
        df.columns = cols

    def pick(*names):
        for n in names:
            if n in df.columns:
                return n
        return None

    c = pick("chrom", "chr", "chromosome") or df.columns[0]
    s = pick("start", "startpos", "begin") or df.columns[1]
    e = pick("end", "endpos", "stop") or df.columns[2]
    t = pick("type", "state", "cnvtype", "eventtype", "alteration", "svtype")
    cn = pick("cn", "copynumber", "expectedcn", "copy", "expectedcopynumber")

    out = pd.DataFrame({
        "chrom": df[c].astype(str).str.replace("chr", "", regex=False),
        "start": pd.to_numeric(df[s], errors="coerce"),
        "end": pd.to_numeric(df[e], errors="coerce"),
    })

    def label(row_t, row_cn):
        txt = str(row_t).upper() if row_t is not None else ""
        if row_cn is not None:
            try:
                v = float(row_cn)
                if v > 2:
                    return 1
                if v == 1:
                    return 2
                if v == 0:
                    return 3
            except ValueError:
                pass
        if "GAIN" in txt or "DUP" in txt or "CN3" in txt:
            return 1
        if "CN0" in txt or "HOMO" in txt:
            return 3
        if "LOSS" in txt or "DEL" in txt or "CN1" in txt:
            return 2
        return None

    other = [x for x in df.columns if x not in (c, s, e)]

    def row_label(i):
        lab = label(df[t].iloc[i] if t else None,
                    df[cn].iloc[i] if cn else None)
        if lab is None:  # no named type column: scan the remaining columns
            for col in other:
                lab = label(df[col].iloc[i], None)
                if lab is not None:
                    break
        return lab

    out["label"] = [row_label(i) for i in range(len(df))]
    bad = out["label"].isna() | out["start"].isna() | out["end"].isna()
    if bad.any():
        sys.exit("ERROR: could not parse these truth rows:\n"
                 f"{df[bad].head().to_string()}\n"
                 f"columns seen: {list(df.columns)}")
    out[["start", "end", "label"]] = out[["start", "end", "label"]].astype(int)
    return out


def chrom_name(bam, wanted):
    names = set(bam.references)
    if wanted in names:
        return wanted
    if "chr" + wanted in names:
        return "chr" + wanted
    return None


def process_chrom(bam, fasta, chrom_bam, chrom_ref, binsize,
                  include_dups, min_mapq):
    L = bam.get_reference_length(chrom_bam)
    nbins = L // binsize + 1

    diff = np.zeros(L + 1, dtype=np.int32)
    mapq_sum = np.zeros(nbins, dtype=np.float64)
    mapq_n = np.zeros(nbins, dtype=np.int64)

    starts, ends = [], []
    for r in bam.fetch(chrom_bam):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        if r.is_duplicate and not include_dups:
            continue
        if r.mapping_quality < min_mapq:
            continue
        s, e = r.reference_start, r.reference_end
        if e is None:
            continue
        starts.append(s)
        ends.append(e)
        b = s // binsize
        mapq_sum[b] += r.mapping_quality
        mapq_n[b] += 1
        if len(starts) >= 2_000_000:
            np.add.at(diff, np.array(starts), 1)
            np.add.at(diff, np.array(ends), -1)
            starts, ends = [], []
    if starts:
        np.add.at(diff, np.array(starts), 1)
        np.add.at(diff, np.array(ends), -1)

    depth = np.cumsum(diff[:L], dtype=np.int64)
    del diff
    pad = nbins * binsize - L
    depth = np.concatenate([depth, np.zeros(pad, dtype=np.int64)])
    rd = depth.reshape(nbins, binsize).mean(axis=1)
    del depth

    seq = fasta.fetch(chrom_ref).upper()
    gc = np.zeros(nbins)
    hasN = np.zeros(nbins, dtype=bool)
    for i in range(nbins):
        chunk = seq[i * binsize:(i + 1) * binsize]
        if not chunk:
            hasN[i] = True
            continue
        hasN[i] = "N" in chunk
        gc[i] = (chunk.count("G") + chunk.count("C")) / len(chunk)

    with np.errstate(invalid="ignore", divide="ignore"):
        mapq = np.where(mapq_n > 0, mapq_sum / np.maximum(mapq_n, 1), np.nan)
    return rd, gc, hasN, mapq


def neighbour_rel(rd):
    n = len(rd)
    rel = np.zeros(n)
    csum = np.concatenate([[0.0], np.cumsum(rd)])
    for i in range(n):
        lo, hi = max(0, i - 5), min(n, i + 6)
        left = (csum[i] - csum[lo]) / (i - lo) if i > lo else None
        right = (csum[hi] - csum[i + 1]) / (hi - i - 1) if hi > i + 1 else None
        if left is not None and right is not None:
            rel[i] = abs(((rd[i] - left) + (rd[i] - right)) / 2)
        elif left is not None:
            rel[i] = abs(rd[i] - left)
        elif right is not None:
            rel[i] = abs(rd[i] - right)
    return rel


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--bin", type=int, default=1000)
    ap.add_argument("--truth", help="truth table for CNV-positive samples")
    ap.add_argument("--drop-zero", action="store_true",
                    help="drop rd==0 bins exactly like the released code")
    ap.add_argument("--include-duplicates", action="store_true", default=True)
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--chroms", default=",".join(AUTOSOMES))
    a = ap.parse_args()

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    truth = read_truth(a.truth) if a.truth else None
    bam = pysam.AlignmentFile(a.bam, "rb")
    fasta = pysam.FastaFile(a.ref)

    genome_rows = []
    for chrom in a.chroms.split(","):
        cb = chrom_name(bam, chrom)
        cr = chrom if chrom in fasta.references else (
            "chr" + chrom if "chr" + chrom in fasta.references else None)
        if cb is None or cr is None:
            sys.exit(f"ERROR: chromosome {chrom} missing from BAM or FASTA")

        rd, gc, hasN, mapq = process_chrom(bam, fasta, cb, cr, a.bin,
                                           a.include_duplicates, a.min_mapq)
        keep = ~hasN
        if a.drop_zero:
            keep &= rd > 0

        idx = np.nonzero(keep)[0]
        rd_k, gc_k = rd[idx], gc[idx]
        mq = mapq[idx]
        mq = np.where(np.isnan(mq), 60.0, mq)  # released code sets empty bins to 60
        qua = mq / mq.max() if mq.max() > 0 else mq
        rel = neighbour_rel(rd_k)

        label = np.zeros(len(idx), dtype=int)
        if truth is not None:
            starts = idx * a.bin
            ends = starts + a.bin
            for ev in truth[truth.chrom == chrom].itertuples():
                ov = (starts < ev.end) & (ends > ev.start)
                # bin is labelled if its midpoint lies inside the event
                mid = starts + a.bin // 2
                label[ov & (mid >= ev.start) & (mid < ev.end)] = ev.label

        binno = idx + 1  # 1-based, as in the released code
        df = pd.DataFrame({"bin": binno, "rd": rd_k, "gc": gc_k,
                           "rel": rel, "qua": qua, "label": label})
        df.to_csv(out / f"{a.sample}.chr{chrom}.bin{a.bin}.txt",
                  sep="\t", header=False, index=False, float_format="%.6g")

        g = df.copy()
        g.insert(0, "chrom", chrom)
        g["start"] = idx * a.bin
        g["end"] = g["start"] + a.bin
        genome_rows.append(g)
        print(f"chr{chrom}: bins kept {len(idx)}, labelled "
              f"gain={int((label == 1).sum())} hemi_loss={int((label == 2).sum())}",
              flush=True)

    allg = pd.concat(genome_rows, ignore_index=True)
    allg.to_csv(out / f"{a.sample}.genome.bin{a.bin}.tsv.gz", sep="\t",
                index=False, float_format="%.6g")
    (out / f".done.bin{a.bin}").touch()
    print("DONE", a.sample)


if __name__ == "__main__":
    main()
