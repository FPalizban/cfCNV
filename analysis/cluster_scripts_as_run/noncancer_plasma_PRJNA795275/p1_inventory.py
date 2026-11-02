#!/usr/bin/env python3
"""
Phase 3 inventory of PRJNA795275 from the ENA portal API (no login needed).

Writes to $S2_META:
  p3_runs_all.tsv      every run of the project with its metadata and sample attributes
  p3_samples_wgs.tsv   the WGS samples used here (exome and other strategies excluded),
                       one row per sample: runs, download URLs, md5, layout, read length,
                       estimated depth and the subsampling fraction to reach ~TARGET_MAX_X
  tasks_download.tsv   one row per FASTQ file to download
  tasks_align.tsv      one row per sample to align
"""
import csv, io, json, os, sys, time, urllib.request, xml.etree.ElementTree as ET
from pathlib import Path

META = Path(os.environ["S2_META"]); FQ = Path(os.environ["P3_FASTQ"])
TARGET = float(os.environ.get("TARGET_MAX_X", "12"))
GENOME = 3.1e9
FIELDS = ("run_accession,sample_accession,secondary_sample_accession,experiment_accession,sample_title,"
          "sample_alias,experiment_title,library_name,library_strategy,library_source,library_selection,"
          "library_layout,instrument_platform,instrument_model,read_count,base_count,fastq_ftp,fastq_md5,"
          "fastq_bytes,submitted_ftp")


def get(url, tries=5):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode()
        except Exception as e:
            print("  retry", i + 1, url, e, file=sys.stderr); time.sleep(5 * (i + 1))
    raise SystemExit(f"cannot reach {url}")


def attributes(sample):
    try:
        x = ET.fromstring(get(f"https://www.ebi.ac.uk/ena/browser/api/xml/{sample}"))
        return {a.findtext("TAG"): a.findtext("VALUE") for a in x.iter("SAMPLE_ATTRIBUTE")}
    except SystemExit:
        return {}


def main():
    url = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA795275&result=read_run"
           f"&fields={FIELDS}&format=tsv&limit=0")
    runs = list(csv.DictReader(io.StringIO(get(url)), delimiter="\t"))
    if not runs:
        raise SystemExit("ENA returned no runs for PRJNA795275")
    att = {}
    for s in sorted({r["sample_accession"] for r in runs}):
        att[s] = attributes(s); time.sleep(0.2)
    keys = sorted({k for a in att.values() for k in a})
    with open(META / "p3_runs_all.tsv", "w") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(list(runs[0]) + [f"attr:{k}" for k in keys])
        for r in runs:
            w.writerow(list(r.values()) + [att[r["sample_accession"]].get(k, "") for k in keys])

    wgs = [r for r in runs if r["library_strategy"].upper() == "WGS"]
    by = {}
    for r in wgs:
        by.setdefault(r["sample_accession"], []).append(r)
    samples, dl = [], []
    for s, rr in sorted(by.items()):
        bases = sum(int(r["base_count"] or 0) for r in rr); reads = sum(int(r["read_count"] or 0) for r in rr)
        layout = rr[0]["library_layout"].upper()
        rlen = round(bases / reads / (2 if layout == "PAIRED" else 1)) if reads else 0
        depth = bases / GENOME
        frac = min(1.0, TARGET / depth) if depth > 0 else 1.0
        files = []
        for r in rr:
            urls = [u for u in r["fastq_ftp"].split(";") if u]
            md5s = r["fastq_md5"].split(";") if r["fastq_md5"] else [""] * len(urls)
            if layout == "PAIRED" and len(urls) == 3:            # ENA lists unpaired reads first
                urls, md5s = urls[1:], md5s[1:]
            for u, m in zip(urls, md5s):
                local = FQ / s / Path(u).name
                files.append(str(local))
                dl.append({"sample": s, "run": r["run_accession"], "url": "https://" + u, "md5": m, "local": str(local)})
            if not urls:
                dl.append({"sample": s, "run": r["run_accession"], "url": "SRA", "md5": "", "local": str(FQ / s / r["run_accession"])})
        a = att[s]
        # study groups as described in Methods 2.3.2 (3 deep WGS, 3 additional WGS, 24 shallow WGS):
        # the additional WGS samples form a separate BioSample series (SAMN245920xx); among the rest,
        # deep = > 10x raw depth, shallow = the remainder
        group = "additional_wgs" if s.startswith("SAMN245920") else ("deep_wgs" if depth > 10 else "shallow_wgs")
        samples.append({"sample_id": s, "study_group": group, "title": rr[0]["sample_title"], "alias": rr[0]["sample_alias"],
                        "runs": ";".join(r["run_accession"] for r in rr), "layout": layout,
                        "instrument": rr[0]["instrument_model"], "read_len": rlen,
                        "raw_depth_x": round(depth, 3), "subsample_fraction": round(frac, 4),
                        "fastq_files": ";".join(files),
                        "attributes": json.dumps(a, ensure_ascii=False)})
    with open(META / "p3_samples_wgs.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(samples[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(samples)
    with open(META / "tasks_download.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=["sample", "run", "url", "md5", "local"], delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(dl)
    with open(META / "tasks_align.tsv", "w") as f:
        f.write("sample_id\tfastq_files\tlayout\tparam\n")
        for s in samples:
            f.write(f"{s['sample_id']}\t{s['fastq_files']}\t{s['layout']}\t{s['subsample_fraction']}\n")

    strat = {}
    for r in runs:
        strat[r["library_strategy"]] = strat.get(r["library_strategy"], 0) + 1
    gb = sum(int(x or 0) for r in wgs for x in (r["fastq_bytes"] or "0").split(";")) / 1e9
    print("runs by library strategy:", strat)
    from collections import Counter
    print("WGS study groups:", dict(Counter(x["study_group"] for x in samples)), "(Methods: 3 deep, 3 additional, 24 shallow)")
    print(f"WGS samples: {len(samples)}  |  FASTQ files to download: {len(dl)}  |  download size: {gb:.0f} GB")
    print(f"{'sample':14s} {'group':15s} {'layout':7s} {'len':>4s} {'depth':>8s} {'keep':>6s}  title")
    for s in samples:
        print(f"{s['sample_id']:14s} {s['study_group']:15s} {s['layout']:7s} {s['read_len']:4d} {s['raw_depth_x']:7.2f}x {s['subsample_fraction']:6.3f}  {s['title'][:60]}")


if __name__ == "__main__":
    main()
