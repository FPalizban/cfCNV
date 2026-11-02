#!/usr/bin/env python3
"""Helpers to run the pretrained LYCEUM model (hg38 exon targets) on cfDNA BAMs.

  targets  write the exon BED to measure depth on, in the BAM's build and chromosome naming.
           GRCh38: the released hg38 targets (chr prefix adapted to the BAM).
           GRCh37: exons lifted one-to-one with the UCSC hg38ToHg19 chain (whole exon inside one '+'
           strand block, same length), as in the paper's Simulation 1 adaptation; a map is written.
  depth    convert a 'sambamba depth base' file to the hg38 coordinates and 'chr' names LYCEUM expects
           (positions are shifted exon by exon, so the model sees its native positional input).
  calls    turn LYCEUM's exon-level output into gain/loss segments (same-state exons <= 100 kb apart are
           merged; a NO-CALL exon ends a segment, as in the paper) in the BAM's build.

usage:
  lyceum_tools.py targets --hg38-bed T.bed --bam-chroms chroms.txt --out-bed targets.bam.bed --out-map map.tsv [--chain C.gz]
  lyceum_tools.py depth   --in depth.txt --map map.tsv --out read_depths/LYC.bam.txt --out-bed targets.hg38.bed
  lyceum_tools.py calls   --in LYC_exonlevel.csv --map map.tsv --out segments.tsv [--merge-gap 100000]
"""
import argparse, bisect, gzip, sys
from collections import defaultdict


def norm(c): return c[3:] if c.startswith("chr") else c


def load_chain(path):
    blocks = defaultdict(list)   # hg38 chrom -> (t0, t1, q0, qchrom, qstrand, score)
    with gzip.open(path, "rt") as f:
        for l in f:
            x = l.split()
            if not x: continue
            if x[0] == "chain":
                sc, tN, tSt, tB = float(x[1]), x[2], x[4], int(x[5]); qN, qSt, qB = x[7], x[9], int(x[10])
                t, q = tB, qB; hdr = (sc, tN, qN, qSt, tSt); continue
            sc, tN, qN, qSt, tSt = hdr; n = int(x[0])
            if tSt == "+": blocks[tN].append((t, t + n, q, qN, qSt, sc))
            if len(x) == 3: t += n + int(x[1]); q += n + int(x[2])
    for k in blocks: blocks[k].sort()
    return blocks


def cmd_targets(a):
    have = {l.strip() for l in open(a.bam_chroms) if l.strip()}
    pref = "chr" if any(c.startswith("chr") for c in have) else ""
    rows = [l.split()[:3] for l in open(a.hg38_bed) if l.strip()]
    out, mp, n_drop = [], [], 0
    if not a.chain:      # same build
        for c, s, e in rows:
            bc = pref + norm(c)
            if bc in have:
                out.append((bc, int(s), int(e))); mp.append((c, s, e, bc, s, e))
    else:
        B = load_chain(a.chain); starts = {k: [b[0] for b in v] for k, v in B.items()}
        for c, s, e in rows:
            s, e = int(s), int(e); bl = B.get(c, []); i = bisect.bisect_right(starts.get(c, []), s) - 1
            best, k = None, 0
            while 0 <= i < len(bl) and k < 200:          # blocks starting at or before the exon
                t0, t1, q0, qN, qSt, sc = bl[i]
                if t1 >= e and qSt == "+" and (best is None or sc > best[5]): best = bl[i]
                i -= 1; k += 1
            if best is None: n_drop += 1; continue
            t0, t1, q0, qN, qSt, sc = best; qs = q0 + (s - t0); bc = pref + norm(qN)
            if bc not in have or norm(qN) != norm(c): n_drop += 1; continue
            out.append((bc, qs, qs + (e - s))); mp.append((c, s, e, bc, qs, qs + (e - s)))
        print(f"lifted {len(out)} of {len(rows)} exons one-to-one; {n_drop} dropped", file=sys.stderr)
    with open(a.out_bed, "w") as f:
        for r in sorted(out, key=lambda r: (r[0], r[1])): f.write("%s\t%d\t%d\n" % r)
    with open(a.out_map, "w") as f:
        f.write("chrom38\tstart38\tend38\tchrom_bam\tstart_bam\tend_bam\n")
        for r in mp: f.write("\t".join(map(str, r)) + "\n")


def read_map(path):
    m = [l.rstrip("\n").split("\t") for l in open(path)][1:]
    return [(r[0], int(r[1]), int(r[2]), r[3], int(r[4]), int(r[5])) for r in m]


def cmd_depth(a):
    M = read_map(a.map); by = defaultdict(list)
    for r in M: by[r[3]].append(r)
    st = {k: [r[4] for r in sorted(v, key=lambda r: r[4])] for k, v in by.items()}
    by = {k: sorted(v, key=lambda r: r[4]) for k, v in by.items()}
    with open(a.out_bed, "w") as f:
        for r in sorted(M, key=lambda r: (r[0], r[1])): f.write(f"{r[0]}\t{r[1]}\t{r[2]}\n")
    n = 0
    with open(a.input) as fin, open(a.out, "w") as fo:
        hdr = fin.readline().rstrip("\n").split("\t")
        hdr[0] = "REF"; fo.write("\t".join(hdr) + "\n")
        for l in fin:
            p = l.rstrip("\n").split("\t"); c, pos = p[0], int(p[1])
            L = st.get(c)
            if not L: continue
            i = bisect.bisect_right(L, pos) - 1
            # exons do not overlap in the released target set; check the closest one(s)
            for j in (i, i - 1):
                if 0 <= j < len(L):
                    r = by[c][j]
                    if r[4] <= pos < r[5]:
                        p[0] = r[0]; p[1] = str(r[1] + pos - r[4]); fo.write("\t".join(p) + "\n"); n += 1; break
    print(f"{n} per-base depth rows written in hg38 coordinates", file=sys.stderr)


def cmd_calls(a):
    M = {(r[0], r[1], r[2]): r for r in read_map(a.map)} if a.map else {}
    rows = []
    for l in open(a.input):
        p = l.rstrip("\n").split("\t")
        if len(p) < 5: continue
        try: s, e = int(float(p[2])), int(float(p[3]))
        except ValueError: continue
        lab = p[4].upper(); st = "loss" if "DEL" in lab else "gain" if "DUP" in lab else None
        r = M.get((p[1], s, e))
        if a.map and r is None: continue
        c, s2, e2 = (norm(r[3]), r[4], r[5]) if r else (norm(p[1]), s, e)
        rows.append((c, s2, e2, st))
    rows.sort(key=lambda r: (int(r[0]) if r[0].isdigit() else 99, r[1]))
    segs, cur = [], None
    for c, s, e, st in rows:
        if st is None:
            if cur: segs.append(cur); cur = None
            continue
        if cur and cur[0] == c and cur[3] == st and s - cur[2] <= a.merge_gap: cur[2] = max(cur[2], e); cur[4] += 1
        else:
            if cur: segs.append(cur)
            cur = [c, s, e, st, 1]
    if cur: segs.append(cur)
    with open(a.out, "w") as f:
        f.write("chrom\tstart\tend\ttype\tn_exons\n")
        for x in segs: f.write("\t".join(map(str, x)) + "\n")
    print(f"{len(rows)} exons -> {len(segs)} LYCEUM segments", file=sys.stderr)


ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest="cmd", required=True)
t = sub.add_parser("targets"); t.add_argument("--hg38-bed", required=True); t.add_argument("--bam-chroms", required=True)
t.add_argument("--out-bed", required=True); t.add_argument("--out-map", required=True); t.add_argument("--chain")
d = sub.add_parser("depth"); d.add_argument("--in", dest="input", required=True); d.add_argument("--map", required=True)
d.add_argument("--out", required=True); d.add_argument("--out-bed", required=True)
c = sub.add_parser("calls"); c.add_argument("--in", dest="input", required=True); c.add_argument("--map")
c.add_argument("--out", required=True); c.add_argument("--merge-gap", type=int, default=100_000)
a = ap.parse_args()
{"targets": cmd_targets, "depth": cmd_depth, "calls": cmd_calls}[a.cmd](a)
