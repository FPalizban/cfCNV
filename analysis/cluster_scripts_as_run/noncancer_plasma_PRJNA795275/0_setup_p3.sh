#!/bin/bash
# =====================================================================
# 0_setup_p3.sh - run ONCE on the LOGIN node (needs internet), ~10-20 min:
#   bash 0_setup_p3.sh
# 1. finds the BWA and SAMtools modules and writes p3_modules.sh
# 2. checks the GRCh38 BWA index (if missing, 1_submit_stage1.sh builds it first)
# 3. installs Fragle (environment from Zenodo + code from GitHub) into $TOOLS
# 4. lists PRJNA795275 from ENA and selects the WGS samples (p1_inventory.py)
# It must end with SETUP_P3_OK.
# =====================================================================
set -uo pipefail
source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/phase3/full/p3_env.sh
cd "$S2"; bad=0
say() { printf '%-34s %s\n' "$1" "$2"; }

# ---- 1. modules -----------------------------------------------------
pick() { module -t avail "$1" 2>&1 | grep -iE "^$1/" | sort -V | tail -1; }
BWA_MOD=$(pick BWA); SAM_MOD="SAMtools/1.16.1-GCC-11.3.0"
if [[ -z "$BWA_MOD" ]] && ! command -v bwa >/dev/null; then say "BWA" "MISSING (no module, not in PATH)"; bad=1
else say "BWA" "${BWA_MOD:-$(command -v bwa)}"; fi
{ echo "module load $SAM_MOD"; [[ -n "$BWA_MOD" ]] && echo "module load $BWA_MOD"; } > "$S2_META/p3_modules.sh"
( source "$S2_META/p3_modules.sh" >/dev/null 2>&1; command -v bwa >/dev/null && command -v samtools >/dev/null ) \
  && say "BWA + SAMtools load together" OK || { say "BWA + SAMtools load together" "PROBLEM - send me this output"; bad=1; }
[[ -x "$MOSDEPTH" ]] && say mosdepth OK || { say mosdepth MISSING; bad=1; }

# ---- 2. reference index ----------------------------------------------
[[ -s "$REF38" ]] && say "GRCh38 FASTA" OK || { say "GRCh38 FASTA" "MISSING $REF38"; bad=1; }
if [[ -s "$REF38.bwt" ]]; then say "BWA index" "present"; else say "BWA index" "absent - will be built (job X0, ~2-3 h)"; fi
[[ -s "$P/00_metadata/simulation2/full/spectre_help.txt" ]] && cp "$P/00_metadata/simulation2/full/spectre_help.txt" "$S2/"
for b in 30000 50000 100000; do
  [[ -s "$P/07_tool_outputs/simulation2/CNVkit/reference/flat_GRCh38_${b}.cnn" ]] || { say "CNVkit flat reference ${b}" MISSING; bad=1; }
done

# ---- 3. Fragle ------------------------------------------------------------
if [[ ! -s "$FRAGLE_SRC/main.py" ]]; then
  ( cd "$TOOLS" && wget -q -O FRAGLE-main.zip https://github.com/skandlab/FRAGLE/archive/refs/heads/main.zip \
    && unzip -oq FRAGLE-main.zip && rm -f FRAGLE-main.zip ) || say "Fragle code" "DOWNLOAD FAILED"
fi
if [[ ! -s "$FRAGLE_ENV/bin/activate" ]]; then
  mkdir -p "$FRAGLE_ENV"
  ( cd "$TOOLS" && wget -q -O Fragle.tar.gz "https://zenodo.org/records/13968359/files/Fragle.tar.gz?download=1" \
    && tar -xzf Fragle.tar.gz -C "$FRAGLE_ENV" && rm -f Fragle.tar.gz ) || say "Fragle environment" "DOWNLOAD FAILED"
fi
if [[ ! -f "$FRAGLE_ENV/.unpacked" && -s "$FRAGLE_ENV/bin/activate" ]]; then
  ( set +u; source "$FRAGLE_ENV/bin/activate" && conda-unpack && touch "$FRAGLE_ENV/.unpacked" ) >/dev/null 2>&1 || say "Fragle conda-unpack" FAILED
fi
if [[ -s "$FRAGLE_SRC/main.py" && -s "$FRAGLE_ENV/bin/activate" ]] && \
   ( set +u; source "$FRAGLE_ENV/bin/activate" && python -c "import torch, sklearn, pysam, pandas" ) >/dev/null 2>&1; then
  say Fragle OK
else say Fragle "NOT READY (the CNV callers still run; send me this output)"; fi

# ---- SRA toolkit (only needed for runs without an ENA FASTQ copy)
SRA_MOD=$(pick SRA-Toolkit); say "SRA toolkit" "${SRA_MOD:-not available}"
[[ -n "$SRA_MOD" ]] && echo "module load $SRA_MOD" > "$S2_META/p3_sra_module.sh"

# ---- 4. inventory -----------------------------------------------------------
echo; echo "== PRJNA795275 inventory (ENA) =="
"$SYSPY" "$S2/p1_inventory.py" || { echo "INVENTORY FAILED"; bad=1; }
echo
[[ $bad -eq 0 ]] && echo "SETUP_P3_OK" || echo "SETUP_P3: fix the MISSING / FAILED lines above (or send me this output)"
