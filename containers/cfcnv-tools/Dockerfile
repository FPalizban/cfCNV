# cfCNV helper tools: samtools, bedtools, mosdepth, python (pandas/numpy) for spike-in, down-sampling,
# mixtures, harmonisation, consensus, annotation and scoring.
FROM mambaorg/micromamba:1.5.10
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH
LABEL org.opencontainers.image.source="https://github.com/FPalizban/cfCNV" \
      org.opencontainers.image.description="cfCNV helper tools"
