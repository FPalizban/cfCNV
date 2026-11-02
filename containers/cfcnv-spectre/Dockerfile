# Spectre 0.2.1 (CNVCaller run on mosdepth bin coverage)
FROM python:3.11-slim
RUN apt-get update && apt-get install -y --no-install-recommends git procps && rm -rf /var/lib/apt/lists/* \
 && pip install --no-cache-dir "pysam>=0.22" "numpy>=1.24.3,<2" "pandas>=2.0.1" "matplotlib>=3.7.1" "scipy>=1.10.1" \
 && pip install --no-cache-dir "git+https://github.com/fritzsedlazeck/Spectre.git@0.2.1" \
 && spectre --help > /dev/null 2>&1 || spectre -h | head -5
LABEL org.opencontainers.image.source="https://github.com/FPalizban/cfCNV" \
      org.opencontainers.image.description="cfCNV: Spectre"
