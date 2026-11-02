#!/bin/bash
# 0_setup_g7.sh - run ONCE on the LOGIN node (needs internet), ~1 min:   bash 0_setup_g7.sh
# Reuses everything built for Phase 3 (BWA index, modules, Fragle, CNVkit/ichorCNA references)
# and lists the GSE71378 samples in samples.txt. Must end with SETUP_G7_OK.
set -uo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/gse71378/full/p3_env.sh
cd "$S2"; bad=0; P3M="$P/00_metadata/phase3"
for f in p3_modules.sh p3_sra_module.sh; do [[ -s "$P3M/$f" ]] && cp "$P3M/$f" "$S2_META/"; done
[[ -s "$S2_META/p3_modules.sh" ]] && echo "modules: OK (from Phase 3)" || { echo "modules: MISSING"; bad=1; }
[[ -s "$REF38.bwt" ]] && echo "BWA index: OK" || { echo "BWA index: MISSING"; bad=1; }
[[ -s "$FRAGLE_SRC/main.py" ]] && echo "Fragle: OK" || echo "Fragle: not found (callers still run)"
ls -d "$TOOLS"/BICseq2_resources/GRCh38_GIABv3_*mer 2>/dev/null | sed 's/^/BIC-seq2 mappability available: /'
echo; echo "== GSE71378 inventory (ENA) =="
"$SYSPY" "$S2/p1_inventory.py" || { echo "INVENTORY FAILED"; bad=1; }
echo; [[ $bad -eq 0 ]] && echo SETUP_G7_OK || echo "SETUP_G7: fix the MISSING lines (or send me this output)"
