#!/usr/bin/env python3
"""
Stream FASTQ reads for BWA, optionally subsampled.
  fq_stream.py <fraction> <seed> PAIRED r1.fq.gz r2.fq.gz [r1b r2b ...]   -> interleaved FASTQ on stdout
  fq_stream.py <fraction> <seed> SINGLE r.fq.gz [r2.fq.gz ...]           -> FASTQ on stdout
Pairs are kept or dropped together (one random draw per pair), so mates stay in sync.
"""
import gzip, random, sys


def recs(fh):
    while True:
        h = fh.readline()
        if not h:
            return
        yield h + fh.readline() + fh.readline() + fh.readline()


def main():
    frac, seed, layout, files = float(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4:]
    rng = random.Random(seed); out = sys.stdout.buffer; keep = frac >= 0.9999
    if layout == "PAIRED":
        for a, b in zip(files[0::2], files[1::2]):
            with gzip.open(a, "rb") as fa, gzip.open(b, "rb") as fb:
                for r1, r2 in zip(recs(fa), recs(fb)):
                    if keep or rng.random() < frac:
                        out.write(r1); out.write(r2)
    else:
        for a in files:
            with gzip.open(a, "rb") as fa:
                for r1 in recs(fa):
                    if keep or rng.random() < frac:
                        out.write(r1)


if __name__ == "__main__":
    main()
