#!/bin/bash
# =====================================================================
# 0_setup_sim2.sh  - LOGIN node, run ONCE (needs internet, ~15 min).
#   bash 0_setup_sim2.sh
# 1. checks inputs   2. installs QDNAseq.hg38 bins + ACE into 06_tools/Rlib
# 3. autosome FASTA + ENCODE hg38 blacklist   4. LYCEUM keep-probabilities copy (if missing)
# 5. task lists for all tools.  Ends with "SETUP_SIM2_OK".
# =====================================================================
set -euo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/simulation2/full/s2_env.sh
cd "$S2"; chmod +x ./*.sh 2>/dev/null || true
module load SAMtools/1.16.1-GCC-11.3.0

echo "== 1. inputs"
ok=1
for f in "$REF38" "$REF38.fai" "$TRUTH38" "$S2_META/HG008_mixtures_v2.tsv" "$S2_OUT/CNVkit/reference/flat_GRCh38_50000.cnn" \
         "$SPECTRE" "$MOSDEPTH" "$READCOUNTER" "$PYD" "$PYL" "$TOOLS/bin/run_bicseq2" "$TOOLS/bin/run_dudeml" \
         "$P/00_metadata/72_dudeml_winstat_bga.py" "$P/00_metadata/cfcnv_ml/90_cnvmeann_features.py" \
         "$SIM1_DUDE_MODELS/depth_1/bin_50000" "$SIM1_MEANN_V2/bin_500000"; do
  if [[ -e "$f" ]]; then echo "  ok      $f"; else echo "  MISSING $f"; ok=0; fi
done
echo "  mixtures built: $(ls "$S2_MIX"/*.done 2>/dev/null | wc -l) (expected 104)"
[[ $ok == 1 ]] || { echo "STOP: something is missing - send Claude this output"; exit 1; }

echo; echo "== 2. R packages (QDNAseq.hg38, ACE) into $R_LIBS_USER"
r_modules
Rscript -e '
lib <- Sys.getenv("R_LIBS_USER"); .libPaths(c(lib, .libPaths()))
rep <- "https://cloud.r-project.org"
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", lib=lib, repos=rep)
if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", lib=lib, repos=rep)
if (!requireNamespace("QDNAseq.hg38", quietly=TRUE)) remotes::install_github("asntech/QDNAseq.hg38@main", lib=lib, upgrade="never")
if (!requireNamespace("ACE", quietly=TRUE)) BiocManager::install("ACE", lib=lib, update=FALSE, ask=FALSE)
for (p in c("QDNAseq","QDNAseq.hg38","ACE")) cat(p, as.character(packageVersion(p)), "\n")'
module purge >/dev/null 2>&1 || true
module load SAMtools/1.16.1-GCC-11.3.0

echo; echo "== 3. reference helpers"
[[ -s "$REF38_AUTO" ]] || samtools faidx "$REF38" $(printf 'chr%s ' {1..22}) > "$REF38_AUTO"
[[ -s "$REF38_AUTO.fai" ]] || samtools faidx "$REF38_AUTO"
if [[ ! -s "$BLACKLIST38" ]]; then
  curl -fsSL https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz | gzip -dc | cut -f1-3 > "$BLACKLIST38"
fi
printf 'chr%s\n' {1..22} > "$S2_META/chromosomes_hg38.txt"
echo "  autosome FASTA, blacklist ($(wc -l < "$BLACKLIST38") regions), chromosome list ready"

echo; echo "== 4. LYCEUM copy that keeps probabilities"
if [[ ! -s "$LYC_REPO/scripts/LYCEUM_call_keepconf.py" ]]; then
  bash "$P/00_metadata/rescue/C_lyceum/C0_prepare_lyceum.sh"
  sed -i 's/os\.rmdir(/(lambda *a: None)(/' "$LYC_REPO/scripts/LYCEUM_call_keepconf.py"
fi
echo "  $LYC_REPO/scripts/LYCEUM_call_keepconf.py"

echo; echo "== 5. Spectre options (blacklist support)"
"$SPECTRE" CNVCaller --help > "$S2/spectre_help.txt" 2>&1 || true
grep -q -- "--blacklist" "$S2/spectre_help.txt" && echo "  --blacklist supported"

echo; echo "== 6. task lists"
pyd "$S2/s0_make_tasks.py"
echo; echo "SETUP_SIM2_OK"
