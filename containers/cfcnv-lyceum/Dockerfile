# LYCEUM (ciceklab/LYCEUM, pinned commit; released pretrained weights, unchanged), CPU versions of the
# package versions in the authors' environment file (PyTorch 1.7, TensorFlow 2.2, einops 0.3, local-attention 1.1.1).
FROM mambaorg/micromamba:1.5.10
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && micromamba clean -a -y
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG LYCEUM_COMMIT=6e843329f0792e10275b0b00b561edaf0ebbc5c5
USER root
ENV PATH=/opt/conda/bin:$PATH LYCEUM_DIR=/opt/LYCEUM
RUN pip install --no-cache-dir torch==1.7.0+cpu -f https://download.pytorch.org/whl/torch_stable.html \
 && pip install --no-cache-dir numpy==1.18.5 pandas==1.0.5 scikit-learn==0.23.1 tensorflow-cpu==2.2.0 \
        protobuf==3.20.3 einops==0.3.0 local-attention==1.1.1 tqdm==4.66.2 \
 && python -c "import urllib.request, tarfile, io, sys; c=sys.argv[1]; d=urllib.request.urlopen('https://github.com/ciceklab/LYCEUM/archive/'+c+'.tar.gz').read(); tarfile.open(fileobj=io.BytesIO(d)).extractall('/opt')" ${LYCEUM_COMMIT} \
 && mv /opt/LYCEUM-${LYCEUM_COMMIT} /opt/LYCEUM && test -s /opt/LYCEUM/models/lyceum_model.pt \
 && chmod -R a+rX /opt/LYCEUM \
 && python -c "import torch, tensorflow, einops, local_attention, sklearn; print('torch', torch.__version__, 'tf', tensorflow.__version__)" \
 && cd /opt/LYCEUM/scripts && python -c "import performer_pytorch; print('LYCEUM modules import')"
USER $MAMBA_USER
LABEL org.opencontainers.image.source="https://github.com/FPalizban/cfCNV" \
      org.opencontainers.image.description="cfCNV: LYCEUM"
