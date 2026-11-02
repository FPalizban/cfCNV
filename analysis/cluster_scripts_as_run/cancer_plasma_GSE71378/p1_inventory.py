#!/usr/bin/env python3
"""
Inventory of GSE71378 / PRJNA291063 (Snyder et al., Cell 2016) from the ENA portal API.
Only the samples named in samples.txt (next to this script, one sample title per line, e.g. IC15)
are used. To add samples later: add their titles to samples.txt, run this script again
("$SYSPY" p1_inventory.py) and then bash 5_resume_stage1.sh.

Writes to $S2_META: p3_runs_all.tsv, p3_samples_wgs.tsv, tasks_download.tsv, tasks_align.tsv
(same formats as Phase 3, so the Phase 3 scripts run unchanged).
"""
import csv, io, json, os, sys, time, urllib.request
from pathlib import Path

META = Path(os.environ["S2_META"]); FQ = Path(os.environ["P3_FASTQ"]); HERE = Path(__file__).resolve().parent
TARGET = float(os.environ.get("TARGET_MAX_X", "12")); GENOME = 3.1e9
FIELDS = ("run_accession,sample_accession,sample_title,library_strategy,library_layout,instrument_model,"
          "read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes")
DEEP = {"IC15": "small-cell lung cancer", "IC17": "hepatocellular carcinoma", "IC20": "squamous-cell lung cancer",
        "IC35": "breast cancer (DCIS; reported as colorectal adenocarcinoma elsewhere in the paper)",
        "IC37": "colorectal adenocarcinoma"}


def get(url, tries=5):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode()
        except Exception as e:
            print("  retry", i + 1, e, file=sys.stderr); time.sleep(5 * (i + 1))
    raise SystemExit(f"cannot reach {url}")


def main():
    want = [l.split()[0] for l in open(HERE / "samples.txt") if l.strip() and not l.startswith("#")]
    runs = list(csv.DictReader(io.StringIO(get(
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA291063&result=read_run"
        f"&fields={FIELDS}&format=tsv&limit=0")), delimiter="\t"))
    if not runs:
        raise SystemExit("ENA returned no runs for PRJNA291063")
    with open(META / "p3_runs_all.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(runs[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(runs)
    by = {}
    for r in runs:
        if r["sample_title"] in want:
            by.setdefault(r["sample_title"], []).append(r)
    missing = [s for s in want if s not in by]
    if missing:
        print("NOT FOUND in PRJNA291063:", ", ".join(missing))
    samples, dl = [], []
    for s in want:
        if s not in by:
            continue
        rr = by[s]
        bases = sum(int(r["base_count"] or 0) for r in rr); reads = sum(int(r["read_count"] or 0) for r in rr)
        layout = rr[0]["library_layout"].upper()
        rlen = round(bases / reads / (2 if layout == "PAIRED" else 1)) if reads else 0
        depth = bases / GENOME; frac = min(1.0, TARGET / depth) if depth > 0 else 1.0
        files = []
        for r in rr:
            urls = [u for u in r["fastq_ftp"].split(";") if u]
            md5s = r["fastq_md5"].split(";") if r["fastq_md5"] else [""] * len(urls)
            if layout == "PAIRED" and len(urls) == 3:
                urls, md5s = urls[1:], md5s[1:]
            for u, m in zip(urls, md5s):
                local = FQ / s / Path(u).name; files.append(str(local))
                dl.append({"sample": s, "run": r["run_accession"], "url": "https://" + u, "md5": m, "local": str(local)})
            if not urls:
                dl.append({"sample": s, "run": r["run_accession"], "url": "SRA", "md5": "", "local": str(FQ / s / r["run_accession"])})
        group = "healthy" if s[:2] in ("IH", "BH") else ("cancer_deep" if s in DEEP else "cancer_shallow")
        samples.append({"sample_id": s, "study_group": group, "title": DEEP.get(s, ""), "alias": s,
                        "runs": ";".join(r["run_accession"] for r in rr), "layout": layout,
                        "instrument": rr[0]["instrument_model"], "read_len": rlen,
                        "raw_depth_x": round(depth, 3), "subsample_fraction": round(frac, 4),
                        "fastq_files": ";".join(files), "attributes": json.dumps({"library_strategy": rr[0]["library_strategy"]})})
    with open(META / "p3_samples_wgs.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(samples[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(samples)
    with open(META / "tasks_download.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=["sample", "run", "url", "md5", "local"], delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(dl)
    with open(META / "tasks_align.tsv", "w") as f:
        f.write("sample_id\tfastq_files\tlayout\tparam\n")
        for x in samples:
            f.write(f"{x['sample_id']}\t{x['fastq_files']}\t{x['layout']}\t{x['subsample_fraction']}\n")
    gb = sum(int(b or 0) for s in want for r in by.get(s, []) for b in (r["fastq_bytes"] or "0").split(";")) / 1e9
    print(f"samples: {len(samples)}  |  FASTQ files: {len(dl)}  |  download size: {gb:.0f} GB")
    print(f"{'sample':8s} {'group':15s} {'layout':7s} {'len':>4s} {'depth':>8s} {'keep':>6s}  files")
    for x in samples:
        print(f"{x['sample_id']:8s} {x['study_group']:15s} {x['layout']:7s} {x['read_len']:4d} {x['raw_depth_x']:7.2f}x "
              f"{x['subsample_fraction']:6.3f}  {len(x['fastq_files'].split(';'))}")


if __name__ == "__main__":
    main()
