# RECIPE ONLY - cfCNV does not publish this image. BIC-seq2 is free for academic use; non-academic use needs a
# licence from the Harvard Office of Technology Development, and its redistribution terms are not stated. Build it
# yourself (see assets/hpc/README_bicseq2.md) and pass it with --bicseq2_container.
# Contents: NBICseq-norm 0.2.4 + NBICseq-seg 0.7.2 (Bioconda), the Perl "unique" read filter from
# samtools-0.1.7a_getUnique-0.1.3 (downloaded from the BIC-seq site), and GEM + UCSC tools for mappability resources.
FROM mambaorg/micromamba:1.5.10
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && micromamba clean -a -y
ARG MAMBA_DOCKERFILE_ACTIVATE=1
USER root
ENV PATH=/opt/conda/bin:$PATH
WORKDIR /
# getUnique "unique" read filter (only its Perl script is used), from the BIC-seq site
RUN wget -q --tries=5 http://compbio.med.harvard.edu/BIC-seq/samtools-0.1.7a_getUnique-0.1.3.tar.gz \
 && tar xzf samtools-0.1.7a_getUnique-0.1.3.tar.gz && rm -f samtools-0.1.7a_getUnique-0.1.3.tar.gz \
 && test -s /samtools-0.1.7a_getUnique-0.1.3/misc/samtools.pl && chmod -R a+rX /samtools-0.1.7a_getUnique-0.1.3
# GEM mappability tools (2013 binary release that includes gem-2-wig)
RUN wget -q --tries=5 -O gem.tbz2 "https://downloads.sourceforge.net/project/gemlibrary/gem-library/Binary%20pre-release%203/GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2" \
 && tar xjf gem.tbz2 && cp GEM-binaries-*/bin/* /usr/local/bin/ && rm -rf gem.tbz2 GEM-binaries-* \
 || echo "WARNING: GEM download failed - only needed when --bicseq2_resources is not given"
ENV BICSEQ2_NORM=/opt/conda/bin/NBICseq-norm.pl \
    BICSEQ2_SEG=/opt/conda/bin/NBICseq-seg.pl \
    BICSEQ2_GETUNIQUE=/samtools-0.1.7a_getUnique-0.1.3/misc/samtools.pl
RUN NBICseq-norm.pl 2>&1 | head -3; NBICseq-seg.pl 2>&1 | head -3; samtools --version | head -1
USER $MAMBA_USER
LABEL org.opencontainers.image.source="https://github.com/FPalizban/cfCNV" \
      org.opencontainers.image.description="cfCNV: BIC-seq2"
