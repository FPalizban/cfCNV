# cfCNV R image: QDNAseq (+ hg19/hg38 bin annotations), ACE, ichorCNA 0.3.2 (GitHub master, commit 5bfc03e, as used for the paper) (+ readCounter)
FROM mambaorg/micromamba:1.5.10
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && micromamba clean -a -y
ARG MAMBA_DOCKERFILE_ACTIVATE=1
USER root
ENV PATH=/opt/conda/bin:$PATH ICHOR_DIR=/opt/ichorCNA
# QDNAseq bin annotations for GRCh38 (not on Bioconductor) and ichorCNA (pinned commit used in the paper; the v0.2.0 tag lacks --genomeBuild)
RUN Rscript -e 'remotes::install_github("asntech/QDNAseq.hg38@v1.0.0", upgrade = "never", dependencies = FALSE)' \
 && git clone https://github.com/broadinstitute/ichorCNA.git /opt/ichorCNA \
 && git -C /opt/ichorCNA checkout 5bfc03ed854f0e93fe5b624c97c1290fa0053837 \
 && R CMD INSTALL /opt/ichorCNA \
 && Rscript -e 'suppressPackageStartupMessages({library(QDNAseq); library(QDNAseq.hg19); library(QDNAseq.hg38); library(ACE); library(ichorCNA)}); cat("R packages OK\n")' \
 && grep -q -- '"--genomeBuild"' /opt/ichorCNA/scripts/runIchorCNA.R \
 && which readCounter
USER $MAMBA_USER
LABEL org.opencontainers.image.source="https://github.com/FPalizban/cfCNV" \
      org.opencontainers.image.description="cfCNV: QDNAseq, ACE, ichorCNA"
