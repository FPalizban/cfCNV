#!/bin/bash
set -euo pipefail

PROJECT=/mnt/isilon/cag_inf_op/palizbanf/cfCNV/cfCNV_project

source "$PROJECT/00_metadata/config.sh"
source "$PROJECT/00_metadata/00_env.sh"

echo "=========================================="
echo "cfCNV SIMULATION-1 MASTER SUBMISSION"
echo "=========================================="

mkdir -p \
"$SIM1_MAN" \
"$SIM1_OUT" \
"$SIM1_HARM" \
"$SIM1_MET" \
"$LOGS/simulation1"

# ------------------------------------------------------------
# 1. Syntax checks
# ------------------------------------------------------------

echo
echo "[1/8] Checking scripts..."

for x in \
02_prepare_reference.sh \
10_qdnaseq_array.sh \
11_cnvkit_array.sh \
12_ichorcna_array.sh \
13_spectre_array.sh \
31_score.sh
do
    bash -n "$PROJECT/00_metadata/$x"
done

python3 -m py_compile \
"$PROJECT/00_metadata/01_prepare_manifest.py" \
"$PROJECT/00_metadata/20_collect_calls.py" \
"$PROJECT/00_metadata/21_score_sim1.py" \
"$PROJECT/00_metadata/22_bin_metrics.py"

echo "Syntax checks: OK"

# ------------------------------------------------------------
# 2. Build manifests
# ------------------------------------------------------------

echo
echo "[2/8] Building Simulation-1 manifests..."

python3 "$PROJECT/00_metadata/01_prepare_manifest.py"

MASTER="$SIM1_MAN/master.tsv"

N=$(($(wc -l < "$MASTER") - 1))

if [[ "$N" -ne 114 ]]; then
    echo "ERROR: expected 114 samples; found $N"
    exit 10
fi

# ------------------------------------------------------------
# 3. Normalize class labels
# ------------------------------------------------------------

echo
echo "[3/8] Normalizing class labels..."

python3 - <<'PY'
import csv
import os
from pathlib import Path

d = Path(os.environ["SIM1_MAN"])

files = [
    "master.tsv",
    "qdnaseq_tasks.tsv",
    "cnvkit_tasks.tsv",
    "ichor_tasks.tsv",
    "spectre_tasks.tsv"
]

for name in files:
    p = d / name

    with open(p, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        rows = list(r)
        fields = r.fieldnames

    for row in rows:
        x = row["class"].strip().lower()

        if (
            "neg" in x or
            "control" in x or
            x == "normal"
        ):
            row["class"] = "negative"
        else:
            row["class"] = "positive"

    tmp = p.with_suffix(".tmp")

    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n"
        )
        w.writeheader()
        w.writerows(rows)

    tmp.replace(p)
PY

POS=$(awk -F'\t' 'NR>1 && $6=="positive"{n++} END{print n+0}' "$MASTER")
NEG=$(awk -F'\t' 'NR>1 && $6=="negative"{n++} END{print n+0}' "$MASTER")

echo "Positive = $POS"
echo "Negative = $NEG"

if [[ "$POS" -ne 95 || "$NEG" -ne 19 ]]; then
    echo "ERROR: expected 95 positive + 19 negative."
    exit 11
fi

# ------------------------------------------------------------
# 4. Confirm BAMs
# ------------------------------------------------------------

echo
echo "[4/8] Checking all 114 BAM files..."

missing=0

while IFS=$'\t' read -r sample bam background depth rep class
do
    [[ "$sample" == "sample_id" ]] && continue

    if [[ ! -s "$bam" ]]; then
        echo "MISSING BAM: $bam"
        missing=$((missing+1))
    fi

    if [[ ! -s "${bam}.bai" && ! -s "${bam%.bam}.bai" ]]; then
        echo "WARNING: BAM index not found for $bam"
    fi

done < "$MASTER"

if [[ "$missing" -ne 0 ]]; then
    echo "ERROR: $missing BAMs missing."
    exit 12
fi

echo "All BAMs present."

# ------------------------------------------------------------
# 5. Check required caller resources
# ------------------------------------------------------------

echo
echo "[5/8] Checking caller resources..."

[[ -x "$READCOUNTER" ]] || {
    echo "ERROR: readCounter missing"
    exit 20
}

[[ -x "$MOSDEPTH" ]] || {
    echo "ERROR: mosdepth missing"
    exit 21
}

[[ -x "$SPECTRE" ]] || {
    echo "ERROR: Spectre wrapper missing"
    exit 22
}

for f in \
"$ICHOR_EXT/gc_hg19_1000kb.wig" \
"$ICHOR_EXT/map_hg19_1000kb.wig" \
"$ICHOR_EXT/GRCh37.p13_centromere_UCSC-gapTable.txt"
do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing ichorCNA resource:"
        echo "$f"
        exit 23
    fi
done

echo "Caller resources: OK"

# ------------------------------------------------------------
# 6. Fix CNVkit flat-reference invocation
# ------------------------------------------------------------

echo
echo "[6/8] Checking CNVkit WGS command..."

CNVK="$PROJECT/00_metadata/11_cnvkit_array.sh"

# Current script has: -m wgs -f
# For tumor-only WGS, create flat reference with -n.
if grep -q -- '-m wgs -f' "$CNVK"; then
    sed -i 's/-m wgs -f/-m wgs -n -f/' "$CNVK"
fi

bash -n "$CNVK"

grep 'cnvkit.py batch' "$CNVK" || true

echo "CNVkit command prepared."

# ------------------------------------------------------------
# 7. Submit reference + caller arrays
# ------------------------------------------------------------

echo
echo "[7/8] Submitting jobs..."

NQ=$(($(wc -l < "$SIM1_MAN/qdnaseq_tasks.tsv") - 1))
NC=$(($(wc -l < "$SIM1_MAN/cnvkit_tasks.tsv") - 1))
NI=$(($(wc -l < "$SIM1_MAN/ichor_tasks.tsv") - 1))
NS=$(($(wc -l < "$SIM1_MAN/spectre_tasks.tsv") - 1))

echo "QDNAseq tasks  = $NQ"
echo "CNVkit tasks   = $NC"
echo "ichorCNA tasks = $NI"
echo "Spectre tasks  = $NS"

[[ "$NQ" -eq 684 ]] || exit 30
[[ "$NC" -eq 684 ]] || exit 31
[[ "$NI" -eq 114 ]] || exit 32
[[ "$NS" -eq 114 ]] || exit 33

JREF=$(sbatch --parsable \
"$PROJECT/00_metadata/02_prepare_reference.sh")

JQ=$(sbatch --parsable \
--dependency=afterok:$JREF \
--array=0-$((NQ-1))%18 \
"$PROJECT/00_metadata/10_qdnaseq_array.sh")

JC=$(sbatch --parsable \
--dependency=afterok:$JREF \
--array=0-$((NC-1))%18 \
"$PROJECT/00_metadata/11_cnvkit_array.sh")

JI=$(sbatch --parsable \
--dependency=afterok:$JREF \
--array=0-$((NI-1))%12 \
"$PROJECT/00_metadata/12_ichorcna_array.sh")

JS=$(sbatch --parsable \
--dependency=afterok:$JREF \
--array=0-$((NS-1))%12 \
"$PROJECT/00_metadata/13_spectre_array.sh")

# BIC-seq2 installation/reference preflight can run independently
JB=$(sbatch --parsable \
"$PROJECT/00_metadata/40_bicseq2_preflight.sh")

# ------------------------------------------------------------
# 8. Automatically score after all core callers succeed
# ------------------------------------------------------------

echo
echo "[8/8] Scheduling automatic scoring..."

JSC=$(sbatch --parsable \
--dependency=afterok:$JQ:$JC:$JI:$JS \
"$PROJECT/00_metadata/31_score.sh")

STATUS="$SIM1_MAN/submitted_jobs.tsv"

cat > "$STATUS" <<TABLE
component	job_id
reference	$JREF
QDNAseq	$JQ
CNVkit	$JC
ichorCNA	$JI
Spectre	$JS
BICseq2_preflight	$JB
scoring	$JSC
TABLE

echo
echo "=========================================="
echo "SUBMISSION COMPLETE"
echo "=========================================="
cat "$STATUS"
echo
echo "You may log out now."
echo
echo "Check later with:"
echo "  squeue -u \$USER"
echo
echo "Results will appear under:"
echo "  $SIM1_HARM"
echo "  $SIM1_MET"
echo "=========================================="
