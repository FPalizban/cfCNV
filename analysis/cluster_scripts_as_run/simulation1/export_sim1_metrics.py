from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from collections import Counter
import os
import re

root = Path.cwd()
out = root / "13_exports/cfCNV_Sim1_METRICS_ONLY.zip"

MAX_FILE = 20 * 1024**2
MAX_TOTAL = 250 * 1024**2

extensions = (
    ".csv", ".tsv", ".json", ".txt", ".xlsx",
    ".csv.gz", ".tsv.gz"
)

exclude_names = (
    "read_depth", "readdepth", "coverage_matrix",
    "exon_predictions", "exon_truth_overlaps",
    "per_base", "per_exon_predictions",
    "processed", "raw_reads", "samtools"
)

summary_pattern = re.compile(
    r"metric|summary|runtime|performance|evaluation|"
    r"benchmark|recall|precision|specificity|"
    r"truth_event|negative_control|f1|mcc|"
    r"by_depth|by_size|by_background",
    re.I
)

candidates = {}

def add(p, priority):
    if not p.is_file():
        return

    name = p.name.lower()

    if not name.endswith(extensions):
        return

    if any(x in name for x in exclude_names):
        return

    candidates[p] = min(priority, candidates.get(p, priority))

# 1. Main benchmarking metrics — highest priority.
for folder in (
    "09_metrics/simulation1",
    "12_tables/simulation1",
):
    d = root / folder
    if d.exists():
        for p in d.rglob("*"):
            add(p, 0)

# 2. Caller-specific final summaries and per-sample metrics.
d = root / "07_tool_outputs/simulation1"

if d.exists():
    for base, dirs, names in os.walk(d):
        dirs[:] = [
            x for x in dirs
            if x.lower() not in {
                "read_depths", "processed_native",
                "processed_hg38", "raw", "tmp", "temp",
                "cache", "__pycache__", "model", "calls",
                "reference_adapted"
            }
        ]

        for name in names:
            if summary_pattern.search(name):
                add(Path(base) / name, 1)

# 3. Sample manifests.
for folder in (
    "00_metadata/cfcnv_ml",
    "00_metadata/simulation1",
):
    d = root / folder
    if d.exists():
        for p in d.glob("*"):
            if p.is_file() and (
                "manifest" in p.name.lower()
                or "samples" in p.name.lower()
            ):
                add(p, 2)

selected = []
omitted = []
total = 0

for p, priority in sorted(
    candidates.items(),
    key=lambda x: (x[1], x[0].stat().st_size, str(x[0]))
):
    size = p.stat().st_size

    if size > MAX_FILE:
        omitted.append((p, size, "file_size_limit"))
        continue

    if total + size > MAX_TOTAL:
        omitted.append((p, size, "archive_size_limit"))
        continue

    selected.append((p, size))
    total += size

if not selected:
    raise SystemExit("No analysis-ready files found")

groups = Counter()

with ZipFile(
    out, "w", ZIP_DEFLATED,
    compresslevel=3, allowZip64=True
) as z:

    for p, size in selected:
        rel = str(p.relative_to(root))
        z.write(p, rel)
        groups[rel.split("/")[0]] += 1

    inventory = "relative_path\tbytes\n"
    inventory += "".join(
        f"{p.relative_to(root)}\t{size}\n"
        for p, size in selected
    )
    z.writestr("INCLUDED_FILES.tsv", inventory)

    skipped = "relative_path\tbytes\treason\n"
    skipped += "".join(
        f"{p.relative_to(root)}\t{size}\t{reason}\n"
        for p, size, reason in omitted
    )
    z.writestr("EXCLUDED_FILES.tsv", skipped)

print("===== METRICS-ONLY EXPORT COMPLETE =====")
print("Files included:", len(selected))
print("Files excluded by size:", len(omitted))
print("Uncompressed size (MB):", round(total / 1e6, 2))
print("ZIP size (MB):", round(out.stat().st_size / 1e6, 2))

for group, count in sorted(groups.items()):
    print(group, count, "files")

print("ZIP:", out)
