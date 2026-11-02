import csv, re, urllib.request
from pathlib import Path
from html.parser import HTMLParser

P = Path.cwd()
M = P / "00_metadata/simulation2"
R = P / "01_raw_data/NIST_GIAB_HG008/benchmark"

with (M / "HG008_mixtures.tsv").open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["sample_id", "depth_x", "tumor_read_fraction"])
    for depth in (0.5, 1, 5, 10):
        for frac in (0, .01, .05, .10, .20, .30):
            name = "HG008_f{:03d}_d{}".format(
                round(frac * 100), str(depth).replace(".", "p")
            )
            w.writerow([name, depth, frac])

BASE = ("https://ftp-trace.ncbi.nlm.nih.gov/"
        "ReferenceSamples/giab/data_somatic/HG008/"
        "Liss_lab/analysis/")

V06 = "NIST_HG008-T_somatic-stvar-CNV_DraftBenchmark_V0.6-20261006/"
V05 = "NIST_HG008-T_somatic-stvar-CNV_DraftBenchmark_V0.5-20260318/"

class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
    def handle_starttag(self, tag, attrs):
        if tag == "a":
            self.links += [v for k, v in attrs if k == "href" and v]

url = None
try:
    html = urllib.request.urlopen(BASE + V06, timeout=30).read().decode(
        "utf-8", errors="replace")
    parser = Links()
    parser.feed(html)
    candidates = [
        x.rsplit("/", 1)[-1] for x in parser.links
        if re.search(
            r"^GRCh38.*somatic-CNV.*PASS.*calls\.bed(?:\.gz)?$",
            x.rsplit("/", 1)[-1], re.I
        )
    ]
    if len(candidates) == 1:
        url = BASE + V06 + candidates[0]
except Exception as e:
    print("V0.6 directory lookup unavailable:", str(e)[:150])

if not url:
    # Explicit, documented backup release, never an invented filename.
    url = (BASE + V05 +
       "GRCh38_HG008-T-V0.5_somatic-CNV_PASS.draftbenchmark.calls.bed")

print("BENCHMARK_URL:", url)

dest = R / url.rsplit("/", 1)[-1]
if not dest.exists() or dest.stat().st_size == 0:
    temp = Path(str(dest) + ".part")
    urllib.request.urlretrieve(url, temp)
    temp.replace(dest)

import gzip
reader = gzip.open if dest.suffix == ".gz" else open
truth = R / "HG008_somatic_CNV_truth_GRCh38.bed"

n = 0
with reader(dest, "rt") as src, truth.open("w") as out:
    for line in src:
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 7:
            raise ValueError("Unexpected CNV benchmark format: " + line[:150])

        chrom = fields[0]
        if not re.fullmatch(r"(?:chr)?(?:[1-9]|1[0-9]|2[0-2])", chrom):
            continue

        start, end = int(fields[1]), int(fields[2])
        cn = float(fields[3])
        if start < 0 or end <= start:
            raise ValueError("Invalid benchmark coordinates")
        if cn == 2:
            continue

        direction = "gain" if cn > 2 else "loss"
        out.write(f"{chrom}\t{start}\t{end}\t{direction}\t{cn:g}\n")
        n += 1

assert n > 0, "No autosomal somatic CNA intervals found"
(M / "truth_source_url.txt").write_text(url + "\n")

print("MIXTURE_MANIFEST_READY: 24")
print("TRUTH_CNA_INTERVALS:", n)
print("TRUTH_FILE:", truth)
