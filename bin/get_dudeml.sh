#!/usr/bin/env bash
# dudeML has no licence file, so cfCNV does not redistribute it. At run time this downloads dudeML.py
# (pinned commit) directly from the authors' GitHub repository into the task folder and applies the two
# compatibility edits used for the cfCNV paper (joblib import; pandas row-wise comparison).
# Prints the folder to use as DUDEML_DIR.
set -euo pipefail
C=${DUDEML_COMMIT:-a4c7c0bc58966f9eda024ecab1729d865fa39819}
D=$PWD/dudeML_src; mkdir -p "$D"
if [[ ! -s "$D/dudeML.py" ]]; then
  python - "$C" "$D/dudeML.py" <<'PY'
import sys, urllib.request
c, out = sys.argv[1:]
src = urllib.request.urlopen(f"https://raw.githubusercontent.com/tomh1lll/dudeML/{c}/dudeML.py").read().decode()
src = src.replace("\r\n", "\n").replace("from sklearn.externals import joblib", "import joblib")
src = src.replace("out_bs_1[out_bs_1 == bs_1]", "out_bs_1[out_bs_1.eq(bs_1, axis=0)]").replace("out_bs_2[out_bs_2 == bs_2]", "out_bs_2[out_bs_2.eq(bs_2, axis=0)]")
assert "import joblib" in src and "eq(bs_1, axis=0)" in src, "dudeML.py differs from the expected version"
open(out, "w").write(src)
PY
fi
echo "$D"
