#!/usr/bin/env python
"""B  Read and fragment counts in 10-kb bins for ONE BAM (row N of $MLO/manifest.tsv).

Counted unit = one fragment: single-end reads, or read 1 of a pair (read 2 is skipped).
Always skipped: unmapped, secondary, supplementary, QC-fail, duplicate.
Columns written (autosomes 1..22, 'chr' removed):
  n_all     fragments with any MAPQ           (with n_q20 gives a mappability proxy)
  n_q20     fragments with MAPQ >= 20         (the read-depth signal)
  n_frag    MAPQ>=20 proper pairs with 0 < |insert| <= 1000 bp (fragment length known)
  n_short   ... of these, fragment <= 150 bp
  n_long    ... of these, fragment  > 150 bp
  sum_flen  sum of their fragment lengths      (mean = sum_flen / n_frag)
Output: $MLO/counts_10kb/<dataset>/<sample_id>.bins10kb.tsv.gz, .flen_hist.tsv, .done
The BAM is opened read-only.
"""
import os, re, sys, json, time
from array import array
from pathlib import Path
import numpy as np, pandas as pd, pysam

BIN = 10000
E = os.environ
AUTO = [str(i) for i in range(1, 23)]
nc = lambda c: re.sub(r"^chr", "", str(c))
SKIP = 0x4 | 0x100 | 0x200 | 0x400 | 0x800


def count_one(r):
    od = Path(E["MLO"]) / "counts_10kb" / r["dataset"]; od.mkdir(parents=True, exist_ok=True)
    s = r["sample_id"]
    done = od / f"{s}.done"
    if done.exists() and (od / f"{s}.bins10kb.tsv.gz").exists():
        print("SKIP (done)", s); return
    t0 = time.time()
    bam = pysam.AlignmentFile(r["bam"], "rb", threads=2)
    hist = np.zeros(1001, dtype=np.int64)
    parts, n_pair, n_tot = [], 0, 0
    for contig, L in zip(bam.references, bam.lengths):
        c = nc(contig)
        if c not in AUTO:
            continue
        pa, pq, pf, fl = array("q"), array("q"), array("q"), array("q")
        for a in bam.fetch(contig):
            f = a.flag
            if f & SKIP:
                continue
            if f & 1:
                if f & 0x80:                       # read 2 of a pair
                    continue
                n_pair += 1
            p = a.reference_start
            pa.append(p)
            if a.mapping_quality >= 20:
                pq.append(p)
                if f & 2:
                    t = abs(a.template_length)
                    if 0 < t <= 1000:
                        pf.append(p); fl.append(t)
        nb = (L + BIN - 1) // BIN
        b_all = np.frombuffer(pa, dtype=np.int64) // BIN if len(pa) else np.zeros(0, np.int64)
        b_q = np.frombuffer(pq, dtype=np.int64) // BIN if len(pq) else np.zeros(0, np.int64)
        b_f = np.frombuffer(pf, dtype=np.int64) // BIN if len(pf) else np.zeros(0, np.int64)
        flen = np.frombuffer(fl, dtype=np.int64) if len(fl) else np.zeros(0, np.int64)
        n_tot += len(pa)
        hist += np.bincount(flen, minlength=1001)[:1001]
        start = np.arange(nb) * BIN
        parts.append(pd.DataFrame({
            "chrom": c, "start": start, "end": np.minimum(start + BIN, L),
            "n_all": np.bincount(b_all, minlength=nb)[:nb], "n_q20": np.bincount(b_q, minlength=nb)[:nb],
            "n_frag": np.bincount(b_f, minlength=nb)[:nb],
            "n_short": np.bincount(b_f[flen <= 150], minlength=nb)[:nb],
            "n_long": np.bincount(b_f[flen > 150], minlength=nb)[:nb],
            "sum_flen": np.bincount(b_f, weights=flen, minlength=nb)[:nb].astype(np.int64)}))
        print(s, c, len(pq), "MAPQ>=20 fragments", flush=True)
    bam.close()
    if not parts:
        sys.exit(f"{s}: no autosomes 1..22 found in the BAM header")
    d = pd.concat(parts, ignore_index=True)
    d["chrom"] = pd.Categorical(d["chrom"], AUTO, ordered=True)
    d = d.sort_values(["chrom", "start"]).astype({"chrom": str})
    tmp = od / f"{s}.bins10kb.tsv.gz.tmp"
    d.to_csv(tmp, sep="\t", index=False, compression="gzip"); tmp.rename(od / f"{s}.bins10kb.tsv.gz")
    pd.DataFrame({"flen": np.arange(1001), "n": hist}).query("n > 0").to_csv(od / f"{s}.flen_hist.tsv", sep="\t", index=False)
    nf = int(hist.sum())
    med = int(np.searchsorted(np.cumsum(hist), nf / 2)) if nf else None
    info = {"sample_id": s, "fragments_counted": int(n_tot), "fragments_mapq20": int(d.n_q20.sum()),
            "paired_fraction": round(n_pair / n_tot, 4) if n_tot else None, "fragments_with_length": nf,
            "median_fragment_length": med, "seconds": round(time.time() - t0)}
    (od / f"{s}.info.json").write_text(json.dumps(info) + "\n")
    done.write_text("ok\n")
    print("DONE", json.dumps(info), flush=True)


if __name__ == "__main__":
    man = pd.read_csv(Path(E["MLO"]) / "manifest.tsv", sep="\t", dtype=str)
    task, ntask = int(sys.argv[1]), int(sys.argv[2])      # this task handles rows task, task+ntask, ...
    for i in range(task, len(man), ntask):
        count_one(man.iloc[i])
