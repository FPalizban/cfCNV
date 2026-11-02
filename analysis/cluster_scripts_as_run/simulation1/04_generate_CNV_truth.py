#!/usr/bin/env python3

import argparse
import random

p = argparse.ArgumentParser()
p.add_argument("--chrom-sizes", required=True)
p.add_argument("--output", required=True)
p.add_argument("--seed", required=True, type=int)
args = p.parse_args()

random.seed(args.seed)

# CNV size -> number of events in each simulation replicate
DESIGN = {
    50_000: 4,
    100_000: 4,
    250_000: 4,
    500_000: 4,
    1_000_000: 4,
    5_000_000: 4,
    10_000_000: 2,
    25_000_000: 2,
    50_000_000: 2
}

chromosomes = []

with open(args.chrom_sizes) as f:
    for line in f:
        chrom, length = line.rstrip().split()[:2]

        # Primary truth events only on autosomes.
        if chrom.isdigit() and 1 <= int(chrom) <= 22:
            chromosomes.append((chrom, int(length)))

events = []

def overlaps_existing(chrom, start, end):
    # Maintain at least 1 Mb separation between independent events.
    padding = 1_000_000

    for e in events:
        ec, es, ee = e[0], e[1], e[2]

        if ec != chrom:
            continue

        if start < ee + padding and end > es - padding:
            return True

    return False

event_number = 1

for cnv_size, n_events in DESIGN.items():

    for i in range(n_events):

        for attempt in range(10000):

            chrom, chrom_length = random.choice(chromosomes)

            margin = max(5_000_000, int(cnv_size * 0.2))

            if chrom_length <= cnv_size + 2 * margin:
                continue

            start = random.randint(
                margin,
                chrom_length - cnv_size - margin
            )

            end = start + cnv_size

            if overlaps_existing(chrom, start, end):
                continue

            cnv_type = "loss" if event_number % 2 else "gain"

            copy_number = 1 if cnv_type == "loss" else 3

            size_class = (
                "small" if cnv_size <= 250_000
                else "medium" if cnv_size <= 5_000_000
                else "large"
            )

            events.append(
                (
                    chrom,
                    start,
                    end,
                    cnv_type,
                    copy_number,
                    size_class,
                    f"CNV_{event_number}",
                    cnv_size
                )
            )

            event_number += 1
            break
        else:
            raise RuntimeError(
                f"Unable to place CNV of size {cnv_size}"
            )

# Sort numerically by chromosome and coordinate.
events.sort(key=lambda x: (int(x[0]), x[1]))

with open(args.output, "w") as out:

    out.write(
        "#chrom\tstart\tend\tcnv_type\tcopy_number\t"
        "size_class\tcnv_id\tcnv_size\n"
    )

    for e in events:
        out.write("\t".join(map(str, e)) + "\n")
