# Environment for the machine-learning callers trained by the workflow: the Python CNV-MEANN network
# (bin/cnvmeann_*.py, scikit-learn) and dudeML's dependencies. dudeML itself has no licence file and is NOT
# included: bin/get_dudeml.sh downloads it from the authors' repository when the workflow runs.
FROM mambaorg/micromamba:1.5.10
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH
RUN python -c "import sklearn, pysam, pandas, joblib; print('sklearn', sklearn.__version__, 'pandas', pandas.__version__)"
LABEL org.opencontainers.image.source="https://github.com/FPalizban/cfCNV" \
      org.opencontainers.image.description="cfCNV: environment for CNV-MEANN and dudeML (dudeML downloaded at run time)"
