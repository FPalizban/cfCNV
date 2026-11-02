# Fragle (skandlab/FRAGLE, commit pinned below) with CPU PyTorch; dependencies as listed in its README
FROM python:3.9-slim
ARG FRAGLE_COMMIT=3212220151e4250cd8eeb5f26f50c17b878a55da
RUN apt-get update && apt-get install -y --no-install-recommends git samtools procps && rm -rf /var/lib/apt/lists/* \
 && pip install --no-cache-dir torch==2.2.2 --index-url https://download.pytorch.org/whl/cpu \
 && pip install --no-cache-dir scikit-learn==1.2.2 pandas==1.5.3 numpy==1.24.4 pysam==0.22.0 \
 && git clone https://github.com/skandlab/FRAGLE.git /opt/FRAGLE && cd /opt/FRAGLE && git checkout ${FRAGLE_COMMIT} \
 && chmod -R a+rX /opt/FRAGLE
ENV FRAGLE_DIR=/opt/FRAGLE
LABEL org.opencontainers.image.source="https://github.com/FPalizban/cfCNV" \
      org.opencontainers.image.description="cfCNV: Fragle"
