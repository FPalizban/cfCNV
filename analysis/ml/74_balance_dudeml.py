#!/usr/bin/env python3

import argparse
import random
from pathlib import Path


ap = argparse.ArgumentParser()

ap.add_argument(
    "--files-from",
    required=True
)

ap.add_argument(
    "--output",
    required=True
)

ap.add_argument(
    "--per-class",
    type=int,
    default=10000
)

ap.add_argument(
    "--seed",
    type=int,
    default=2026
)

args = ap.parse_args()

random.seed(args.seed)

classes = ["N", "Dup", "Del"]

reservoir = {
    x: []
    for x in classes
}

seen = {
    x: 0
    for x in classes
}

K = args.per_class


def consider(label, line):

    if label not in reservoir:
        return

    seen[label] += 1

    arr = reservoir[label]

    if len(arr) < K:
        arr.append(line)
        return

    j = random.randint(
        1,
        seen[label]
    )

    if j <= K:
        arr[j - 1] = line


files = [
    Path(x.strip())
    for x in open(args.files_from)
    if x.strip()
]

if not files:
    raise SystemExit(
        "ERROR: no training files"
    )

for f in files:

    if not f.exists():
        raise SystemExit(
            f"ERROR missing {f}"
        )

    with f.open() as fh:

        for line in fh:

            fields = line.rstrip().split("\t")

            if len(fields) < 6:
                continue

            consider(
                fields[3],
                line.rstrip()
            )


for label in classes:

    if seen[label] == 0:
        raise SystemExit(
            f"ERROR: zero {label} examples"
        )


output = []

for label in classes:

    arr = reservoir[label]

    # If fewer than K examples exist,
    # resample with replacement.
    if len(arr) < K:

        original = list(arr)

        while len(arr) < K:
            arr.append(
                random.choice(original)
            )

    output.extend(arr[:K])


random.shuffle(output)

with open(args.output, "w") as out:

    for line in output:
        out.write(line + "\n")


print("Input files:", len(files))

for label in classes:
    print(
        label,
        "observed=",
        seen[label],
        "written=",
        K
    )

print("TOTAL =", len(output))
print("OUTPUT =", args.output)
