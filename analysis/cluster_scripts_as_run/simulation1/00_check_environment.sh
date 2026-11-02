#!/bin/bash

echo "Checking basic commands..."

for cmd in wget git awk sed samtools bedtools fastqc prefetch fasterq-dump bwa java Rscript python3; do
  echo -n "$cmd: "
  which $cmd || echo "NOT FOUND"
done

echo
echo "Loaded modules:"
module list 2>&1 || true
