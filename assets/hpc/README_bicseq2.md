# BIC-seq2 image (build it yourself)

BIC-seq2 (Xi et al. 2016, Park lab) is free for academic use; non-academic use requires a licence from the
Harvard Office of Technology Development, and the BIC-seq site does not state whether the software may be
redistributed. cfCNV therefore distributes only the recipe (`containers/cfcnv-bicseq2/`), not an image. Each
user builds it from the original sources (Bioconda packages and the BIC-seq site) under the BIC-seq2 terms.

**With Docker (any computer):**
```bash
docker build -t cfcnv-bicseq2:local containers/cfcnv-bicseq2
# use with -profile docker:      --bicseq2_container cfcnv-bicseq2:local
# or convert for an HPC:         singularity build cfcnv-bicseq2.sif docker-daemon://cfcnv-bicseq2:local
```

**On an HPC with Singularity/Apptainer (fakeroot enabled):**
```bash
singularity build --fakeroot cfcnv-bicseq2.sif assets/hpc/cfcnv-bicseq2.def
```

Then run cfCNV with `--bicseq2_container /path/to/cfcnv-bicseq2.sif`. Without it, BIC-seq2 is skipped with a warning.
