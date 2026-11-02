# CNVkit 0.9.10 (version used in the paper)
FROM mambaorg/micromamba:1.5.10
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH
RUN cnvkit.py version
LABEL org.opencontainers.image.source="https://github.com/FPalizban/cfCNV" \
      org.opencontainers.image.description="cfCNV: CNVkit"
