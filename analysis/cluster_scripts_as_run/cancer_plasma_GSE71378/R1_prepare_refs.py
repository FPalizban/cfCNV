#!/usr/bin/env python3
"""Condense the raw GRCh37 references into small tables (standard library only).
usage: R1_prepare_refs.py <raw_dir> <out_dir>"""
import bisect, csv, gzip, json, os, re, sys
from collections import defaultdict
raw, out = sys.argv[1:]
J = lambda f: os.path.join(raw, f); O = lambda f: os.path.join(out, f)
COHORTS = ["LIHC", "LUSC", "LUAD", "COAD", "READ", "BRCA"]
AUTO = [str(i) for i in range(1, 23)]
op = lambda f: gzip.open(f, "rt") if f.endswith(".gz") else open(f)
ok = lambda f: os.path.exists(f) and os.path.getsize(f) > 0

# chromosome sizes (GRCh37) from cytobands; arm boundaries
size, arms = {}, {}
if ok(J("hg19_cytoband.txt.gz")):
    for l in op(J("hg19_cytoband.txt.gz")):
        c, s, e, band, stain = l.rstrip("\n").split("\t")[:5]; c = c.replace("chr", "")
        if c not in AUTO: continue
        size[c] = max(size.get(c, 0), int(e))
        if stain == "acen" and band.startswith("p"): arms[c] = int(e)
    with open(O("arms_hg19.tsv"), "w") as f:
        f.write("chrom\tarm\tstart\tend\n")
        for c in AUTO:
            f.write(f"{c}\t{c}p\t0\t{arms[c]}\n{c}\t{c}q\t{arms[c]}\t{size[c]}\n")
    print("arms_hg19.tsv")

# TCGA arm calls -> per cohort gain / loss frequency per arm
type_of = {}
if ok(J("tcga_arm_calls.txt")):
    rows = list(csv.DictReader(open(J("tcga_arm_calls.txt")), delimiter="\t"))
    armcols = [k for k in rows[0] if re.fullmatch(r"\d+[pq]", k)]
    cnt = defaultdict(lambda: defaultdict(lambda: [0, 0, 0]))
    for r in rows:
        type_of[r["Sample"][:15]] = r["Type"]
        if r["Type"] not in COHORTS: continue
        for a in armcols:
            v = r[a]
            if v in ("", "NA"): continue
            x = cnt[r["Type"]][a]; x[2] += 1
            if float(v) > 0: x[0] += 1
            elif float(v) < 0: x[1] += 1
    with open(O("tcga_arm_freq.tsv"), "w") as f:
        f.write("cohort\tarm\tn\tgain_freq\tloss_freq\n")
        for t in COHORTS:
            for a in armcols:
                g, l, n = cnt[t][a]
                if n: f.write(f"{t}\t{a}\t{n}\t{g/n:.4f}\t{l/n:.4f}\n")
    print("tcga_arm_freq.tsv", {t: max(x[2] for x in cnt[t].values()) for t in COHORTS if cnt[t]})

# TCGA SNP6 segments -> per cohort 1-Mb window gain / loss frequency (|log2 ratio| >= 0.2, primary tumours)
if ok(J("tcga_snp6_whitelisted.seg")) and type_of and size:
    segs = defaultdict(lambda: defaultdict(list))
    with open(J("tcga_snp6_whitelisted.seg")) as f:
        next(f)
        for l in f:
            s, c, a, b, n, m = l.rstrip("\n").split("\t")[:6]
            t = type_of.get(s[:15])
            if t not in COHORTS or s[13:15] != "01" or c not in AUTO: continue
            segs[s][c].append((int(float(a)), int(float(b)), float(m)))
    W = 1_000_000; acc = defaultdict(lambda: [0, 0, 0])
    for s, bychr in segs.items():
        t = type_of[s[:15]]
        for c, L in bychr.items():
            L.sort(); starts = [x[0] for x in L]
            for w in range(0, size[c], W):
                mid = w + W // 2; i = bisect.bisect_right(starts, mid) - 1
                if i < 0 or L[i][1] < mid: continue
                x = acc[(t, c, w)]; x[2] += 1
                if L[i][2] >= 0.2: x[0] += 1
                elif L[i][2] <= -0.2: x[1] += 1
    with open(O("tcga_window_freq_1Mb.tsv"), "w") as f:
        f.write("cohort\tchrom\tstart\tend\tn\tgain_freq\tloss_freq\n")
        for (t, c, w), (g, l, n) in sorted(acc.items(), key=lambda k: (k[0][0], int(k[0][1]), k[0][2])):
            f.write(f"{t}\t{c}\t{w}\t{w+W}\t{n}\t{g/n:.4f}\t{l/n:.4f}\n")
    print("tcga_window_freq_1Mb.tsv", len(segs), "tumours")

# genes
if ok(J("gencode.v19.gtf.gz")):
    with op(J("gencode.v19.gtf.gz")) as f, open(O("genes_hg19.tsv"), "w") as g:
        g.write("chrom\tstart\tend\tgene\tgene_type\n")
        for l in f:
            if l.startswith("#"): continue
            x = l.split("\t")
            if x[2] != "gene": continue
            c = x[0].replace("chr", "")
            if c not in AUTO + ["X"]: continue
            name = re.search(r'gene_name "([^"]+)"', x[8]).group(1); typ = re.search(r'gene_type "([^"]+)"', x[8]).group(1)
            if typ == "protein_coding": g.write(f"{c}\t{x[3]}\t{x[4]}\t{name}\t{typ}\n")
    print("genes_hg19.tsv")

# OncoKB
if ok(J("oncokb_cancer_genes.json")):
    d = json.load(open(J("oncokb_cancer_genes.json")))
    with open(O("cancer_genes.tsv"), "w") as g:
        g.write("gene\toncogene\ttsg\n")
        for r in d:
            gt = str(r.get("geneType") or "").upper()      # ONCOGENE / TSG / ONCOGENE_AND_TSG / ... (newer API)
            onc = bool(r.get("oncogene")) or "ONCOGENE" in gt; tsg = bool(r.get("tsg")) or "TSG" in gt
            g.write(f"{r.get('hugoSymbol')}\t{int(onc)}\t{int(tsg)}\n")
    print("cancer_genes.tsv", len(d))

# gnomAD-SV common copy-number variants
if ok(J("gnomad_v2.1_sv.bed.gz")):
    with op(J("gnomad_v2.1_sv.bed.gz")) as f, open(O("gnomad_sv_common_cnv.bed"), "w") as g:
        h = f.readline().lstrip("#").rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
        n = 0
        for l in f:
            x = l.rstrip("\n").split("\t")
            if x[ix["svtype"]] not in ("DEL", "DUP", "CNV"): continue
            try: af = float(x[ix["AF"]])
            except ValueError: continue
            if af >= 0.01: g.write(f"{x[0]}\t{x[1]}\t{x[2]}\t{x[ix['svtype']]}\t{af}\n"); n += 1
    print("gnomad_sv_common_cnv.bed", n)

# artefact-prone regions
def bed(src, dst, cols):
    with op(J(src)) as f, open(O(dst), "w") as g:
        for l in f:
            x = l.rstrip("\n").split("\t"); c = x[cols[0]].replace("chr", "")
            if c in AUTO: g.write(f"{c}\t{x[cols[1]]}\t{x[cols[2]]}\n")
    print(dst)
if ok(J("hg19_segdups.txt.gz")): bed("hg19_segdups.txt.gz", "segdups_hg19.bed", (1, 2, 3))
if ok(J("hg19_blacklist.v2.bed.gz")): bed("hg19_blacklist.v2.bed.gz", "blacklist_hg19.bed", (0, 1, 2))
if ok(J("hallmark.gmt")):
    open(O("hallmark.gmt"), "w").write(open(J("hallmark.gmt")).read()); print("hallmark.gmt")
