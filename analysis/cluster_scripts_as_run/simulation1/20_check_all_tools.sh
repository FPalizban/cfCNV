#!/bin/bash

source /mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project/00_metadata/config.sh

echo "============================================"
echo "cfCNV TOOL INSTALLATION AUDIT"
echo "============================================"

echo
echo "===== 1. QDNAseq / R packages ====="

module load R/4.2.2-foss-2022b
module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2

export R_LIBS_USER=$TOOLS/Rlib

Rscript - <<'RS'
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

pkgs <- c(
    "QDNAseq",
    "QDNAseq.hg19",
    "HMMcopy",
    "DNAcopy",
    "ichorCNA"
)

for (p in pkgs) {
    cat(sprintf("%-18s %s\n",
                p,
                requireNamespace(p, quietly=TRUE)))
}
RS

echo
echo "===== 2. CNVkit ====="

module load CNVkit/0.9.10-foss-2022b-R-4.2.2

if command -v cnvkit.py >/dev/null 2>&1
then
    echo -n "CNVkit: READY - "
    cnvkit.py version
else
    echo "CNVkit: MISSING"
fi

echo
echo "===== 3. ichorCNA readCounter ====="

if command -v readCounter >/dev/null 2>&1
then
    echo "readCounter: READY"
    command -v readCounter
else
    RC=$(find "$TOOLS" -type f -name readCounter -perm -111 2>/dev/null | head -1)

    if [[ -n "$RC" ]]
    then
        echo "readCounter: FOUND"
        echo "$RC"
    else
        echo "readCounter: MISSING"
    fi
fi

echo
echo "===== 4. BIC-seq2 ====="

if [[ -d "$TOOLS/classical/BICSEQ2" ]]
then
    echo "BIC-seq2 repository: FOUND"

    find "$TOOLS/classical/BICSEQ2" \
        -type f \( \
        -name 'NBICseq-norm.pl' \
        -o -name 'NBICseq-seg.pl' \
        \) -print
else
    echo "BIC-seq2: MISSING"
fi

echo
echo "===== 5. Spectre ====="

if [[ -x "$TOOLS/envs/spectre/bin/spectre" ]]
then
    echo "Spectre: READY"
    "$TOOLS/envs/spectre/bin/spectre" --version 2>&1 || true

elif command -v spectre >/dev/null 2>&1
then
    echo "Spectre: READY"
    command -v spectre

else
    echo "Spectre: MISSING"
fi

echo
echo "mosdepth:"

if [[ -x "$TOOLS/bin/mosdepth" ]]
then
    "$TOOLS/bin/mosdepth" --version
elif command -v mosdepth >/dev/null 2>&1
then
    mosdepth --version
else
    echo "MISSING"
fi

echo
echo "===== 6. dudeML ====="

if [[ -f "$TOOLS/ML_based/dudeML/dudeML.py" ]]
then
    echo "dudeML repository: FOUND"

    if [[ -x "$TOOLS/envs/dudeml/bin/python" ]]
    then

        "$TOOLS/envs/dudeml/bin/python" - <<'PY'
mods = ["pandas","numpy","sklearn","Bio"]

for m in mods:
    try:
        __import__(m)
        print(m, "OK")
    except Exception as e:
        print(m, "MISSING", e)
PY

    else
        echo "dudeML Python environment: MISSING"
    fi

else
    echo "dudeML: MISSING"
fi

echo
echo "===== 7. CNV-MEANN ====="

if [[ -d "$TOOLS/ML_based/CNV-MEANN" ]]
then
    echo "CNV-MEANN repository: FOUND"

    echo "Python files:"
    find "$TOOLS/ML_based/CNV-MEANN" \
        -maxdepth 2 -type f -name '*.py' | head -10

    echo "MATLAB files:"
    find "$TOOLS/ML_based/CNV-MEANN" \
        -maxdepth 2 -type f -name '*.m' | head -10
else
    echo "CNV-MEANN: MISSING"
fi

echo
echo "MATLAB availability:"
command -v matlab || echo "MATLAB executable not currently in PATH"

echo
echo "===== 8. LYCEUM ====="

if [[ -d "$TOOLS/ML_based/LYCEUM" ]]
then
    echo "LYCEUM repository: FOUND"

    find "$TOOLS/ML_based/LYCEUM" \
        -maxdepth 3 \
        -type f \
        \( -name '*LYCEUM*.py' \
        -o -name '*.yml' \
        -o -name '*.yaml' \
        -o -name 'requirements*.txt' \
        \) | head -30
else
    echo "LYCEUM: MISSING"
fi

echo
echo "============================================"
echo "AUDIT FINISHED"
echo "============================================"
