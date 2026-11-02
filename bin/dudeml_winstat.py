#!/usr/bin/env python3

import argparse
import math
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


def kth(hist, k):
    """0-based kth order statistic from value->count histogram."""

    total = 0

    for value in sorted(hist):

        total += hist[value]

        if total > k:
            return float(value)

    return float(max(hist))


def quantile(hist, q):

    n = sum(hist.values())

    if n == 0:
        return 0.0

    pos = (n - 1) * q

    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))

    vlo = kth(hist, lo)
    vhi = kth(hist, hi)

    if lo == hi:
        return vlo

    return (
        vlo
        + (pos - lo)
        * (vhi - vlo)
    )


def stats(hist):

    n = sum(hist.values())

    if n == 0:
        return 0.0, 0.0, 0.0, 0.0

    mean = (
        sum(v * c for v, c in hist.items())
        / n
    )

    mean2 = (
        sum((v * v) * c for v, c in hist.items())
        / n
    )

    variance = max(
        0.0,
        mean2 - mean * mean
    )

    std = math.sqrt(variance)

    med = quantile(hist, 0.50)

    q25 = quantile(hist, 0.25)
    q75 = quantile(hist, 0.75)

    iqr = q75 - q25

    return med, iqr, mean, std


def main():

    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--bam",
        required=True
    )

    ap.add_argument(
        "--sample",
        required=True
    )

    ap.add_argument(
        "--outdir",
        required=True
    )

    ap.add_argument(
        "--chroms",
        required=True
    )

    ap.add_argument(
        "--bins",
        nargs="+",
        type=int,
        required=True
    )

    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    chroms = [
        x.strip()
        for x in open(args.chroms)
        if x.strip()
    ]

    allowed = set(chroms)

    bins = sorted(args.bins)

    handles = {}

    for b in bins:

        p = outdir / f"{args.sample}.{b}.winStat.bed"

        handles[b] = p.open("w")


    current = None

    chr_nonzero = Counter()

    # bin_hists[bin_size][bin_number] = Counter(depth -> bases)
    bin_hists = {
        b: defaultdict(Counter)
        for b in bins
    }

    chr_end = 0


    def reset():

        nonlocal chr_nonzero
        nonlocal bin_hists
        nonlocal chr_end

        chr_nonzero = Counter()

        bin_hists = {
            b: defaultdict(Counter)
            for b in bins
        }

        chr_end = 0


    def finish_chrom(chrom):

        if chrom is None:
            return

        if chrom not in allowed:
            return

        if not chr_nonzero:

            print(
                f"WARNING: no covered bases on chr {chrom}"
            )

            chrom_median = 1.0

        else:

            chrom_median = quantile(
                chr_nonzero,
                0.50
            )

            if chrom_median <= 0:
                chrom_median = 1.0


        for b in bins:

            n_full = chr_end // b

            out = handles[b]

            for idx in range(n_full):

                hist = bin_hists[b].get(idx)

                if hist is None:

                    hist = Counter({
                        0: b
                    })

                # Ensure exactly b bases represented.
                observed = sum(hist.values())

                if observed < b:
                    hist[0] += (
                        b - observed
                    )

                med, iqr, mean, std = stats(hist)

                med /= chrom_median
                iqr /= chrom_median
                mean /= chrom_median
                std /= chrom_median

                # Coordinates intentionally match dudeML's
                # fvecTrain key convention:
                #
                # 1..b, b+1..2b, ...
                start = idx * b + 1
                end = (idx + 1) * b

                out.write(
                    "\t".join([
                        str(chrom),
                        str(start),
                        str(end),
                        f"{med:.4f}",
                        f"{iqr:.4f}",
                        f"{mean:.4f}",
                        f"{std:.4f}",
                    ])
                    + "\n"
                )

        print(
            f"Finished chr {chrom}; "
            f"nonzero median={chrom_median}"
        )


    cmd = [
        "bedtools",
        "genomecov",
        "-bga",
        "-ibam",
        args.bam,
    ]

    print("Running:", " ".join(cmd))

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1024 * 1024,
    )

    for line in proc.stdout:

        fields = line.rstrip().split("\t")

        if len(fields) < 4:
            continue

        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        depth = int(float(fields[3]))

        if chrom != current:

            if current is not None:
                finish_chrom(current)

            reset()
            current = chrom

        if chrom not in allowed:
            continue

        chr_end = max(
            chr_end,
            end
        )

        length = end - start

        if depth > 0:
            chr_nonzero[depth] += length


        # Split this run across each requested bin size.
        for b in bins:

            pos = start

            while pos < end:

                idx = pos // b

                stop = min(
                    end,
                    (idx + 1) * b
                )

                n = stop - pos

                bin_hists[b][idx][depth] += n

                pos = stop


    if current is not None:
        finish_chrom(current)

    rc = proc.wait()

    if rc != 0:
        raise SystemExit(
            f"bedtools genomecov failed: {rc}"
        )

    for out in handles.values():
        out.close()

    print()
    print("WINSTAT COMPLETE")

    for b in bins:

        p = outdir / f"{args.sample}.{b}.winStat.bed"

        print(
            b,
            p,
            p.stat().st_size
        )


if __name__ == "__main__":
    main()
