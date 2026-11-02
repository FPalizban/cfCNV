# BIC-seq2 (NBICseq-norm 0.2.4, NBICseq-seg 0.7.2, samtools-0.1.7a getUnique 0.1.3 - only its Perl filter is used) + GEM and UCSC tools to
# build mappability resources. Programs are installed at the same paths as in the image used for the paper.
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends build-essential perl r-base-core r-cran-mgcv \
      samtools zlib1g-dev libncurses-dev wget ca-certificates bzip2 procps \
 && rm -rf /var/lib/apt/lists/*
WORKDIR /
RUN wget -q http://compbio.med.harvard.edu/BIC-seq/NBICseq-norm_v0.2.4.tar.gz \
 && wget -q http://compbio.med.harvard.edu/BIC-seq/NBICseq-seg_v0.7.2.tar.gz \
 && wget -q http://compbio.med.harvard.edu/BIC-seq/samtools-0.1.7a_getUnique-0.1.3.tar.gz \
 && for f in *.tar.gz; do tar xzf "$f"; done && rm -f *.tar.gz \
 && (cd /NBICseq-norm_v0.2.4 && make clean >/dev/null 2>&1; make) \
 && (cd /NBICseq-seg_v0.7.2 && make clean >/dev/null 2>&1; make) \
 && test -f /samtools-0.1.7a_getUnique-0.1.3/misc/samtools.pl
# mappability tools (only needed for --bicseq2_resources not given)
RUN wget -q -O gem.tbz2 "https://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%203/GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2/download" \
 && tar xjf gem.tbz2 && cp GEM-binaries-*/bin/* /usr/local/bin/ && rm -rf gem.tbz2 GEM-binaries-* \
 && for t in wigToBigWig bigWigToBedGraph; do wget -q -O /usr/local/bin/$t http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/$t; chmod +x /usr/local/bin/$t; done
ENV \
    BICSEQ2_NORM=/NBICseq-norm_v0.2.4/NBICseq-norm.pl \
    BICSEQ2_SEG=/NBICseq-seg_v0.7.2/NBICseq-seg.pl \
    BICSEQ2_GETUNIQUE=/samtools-0.1.7a_getUnique-0.1.3/misc/samtools.pl
LABEL org.opencontainers.image.source="https://github.com/FPalizban/cfCNV" \
      org.opencontainers.image.description="cfCNV: BIC-seq2"
