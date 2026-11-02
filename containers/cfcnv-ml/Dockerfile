# Machine-learning callers trained by the workflow: dudeML (tomh1lll/dudeML, pinned commit, with the two
# compatibility edits used for the paper: joblib import and pandas row-wise comparison) and the Python
# CNV-MEANN network (bin/cnvmeann_*.py; scikit-learn).
FROM mambaorg/micromamba:1.5.10
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && micromamba clean -a -y
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG DUDEML_COMMIT=a4c7c0bc58966f9eda024ecab1729d865fa39819
USER root
ENV PATH=/opt/conda/bin:$PATH DUDEML_DIR=/opt/dudeML
RUN git clone https://github.com/tomh1lll/dudeML.git /opt/dudeML && cd /opt/dudeML && git checkout ${DUDEML_COMMIT} && rm -rf .git \
 && sed -i 's/\r$//' dudeML.py \
 && sed -i 's/from sklearn.externals import joblib/import joblib/' dudeML.py \
 && sed -i 's/out_bs_1\[out_bs_1 == bs_1\]/out_bs_1[out_bs_1.eq(bs_1, axis=0)]/; s/out_bs_2\[out_bs_2 == bs_2\]/out_bs_2[out_bs_2.eq(bs_2, axis=0)]/' dudeML.py \
 && grep -c "import joblib" dudeML.py && grep -q "eq(bs_1, axis=0)" dudeML.py \
 && chmod -R a+rX /opt/dudeML \
 && python /opt/dudeML/dudeML.py classify -h > /dev/null \
 && python -c "import sklearn, pysam, pandas; print('sklearn', sklearn.__version__, 'pandas', pandas.__version__)"
USER $MAMBA_USER
LABEL org.opencontainers.image.source="https://github.com/FPalizban/cfCNV" \
      org.opencontainers.image.description="cfCNV: dudeML and CNV-MEANN"
