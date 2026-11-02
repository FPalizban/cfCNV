# Fragle (skandlab/FRAGLE, pinned commit) with CPU PyTorch and the dependencies listed in its README
FROM mambaorg/micromamba:1.5.10
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && micromamba clean -a -y
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG FRAGLE_COMMIT=3212220151e4250cd8eeb5f26f50c17b878a55da
USER root
ENV PATH=/opt/conda/bin:$PATH FRAGLE_DIR=/opt/FRAGLE
RUN git clone https://github.com/skandlab/FRAGLE.git /opt/FRAGLE \
 && cd /opt/FRAGLE && git checkout ${FRAGLE_COMMIT} && rm -rf .git \
 && chmod -R a+rwX /opt/FRAGLE \
 && python -c "import torch, sklearn, pysam, pandas, joblib; print('torch', torch.__version__, 'sklearn', sklearn.__version__)" \
 && python -c "import joblib; joblib.load('/opt/FRAGLE/Models/select_model.pkl'); print('Fragle model selector loads')"
USER $MAMBA_USER
LABEL org.opencontainers.image.source="https://github.com/FPalizban/cfCNV" \
      org.opencontainers.image.description="cfCNV: Fragle"
