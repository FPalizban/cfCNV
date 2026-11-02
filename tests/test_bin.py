"""Unit tests for the helper scripts in bin/ (run: pytest -q tests)."""
import subprocess, sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BIN, ASSETS = ROOT / "bin", ROOT / "assets"


def run(*args):
    subprocess.run([sys.executable, *map(str, args)], check=True)


def test_truth_sets_are_reproduced_exactly(tmp_path):
    """make_truth.py with GRCh37 sizes regenerates the published Simulation 1 truth sets byte for byte."""
    for rep in range(1, 6):
        out = tmp_path / f"rep{rep}.tsv"
        run(BIN / "make_truth.py", "--chrom-sizes", ASSETS / "GRCh37.autosomes.sizes", "--replicate", rep, "--output", out)
        assert out.read_text() == (ASSETS / "truth" / f"rep{rep}.truth.tsv").read_text()


def test_scoring_rule(tmp_path):
    truth = ASSETS / "truth" / "rep1.truth.tsv"
    t = pd.read_csv(truth, sep="\t")
    t.columns = [c.lstrip("#") for c in t.columns]
    calls = t[["chrom", "start", "end", "cnv_type"]].head(10).copy()
    calls.columns = ["chrom", "start", "end", "type"]
    calls.loc[0, "end"] = calls.loc[0, "start"] + (calls.loc[0, "end"] - calls.loc[0, "start"]) // 3   # < 50% overlap
    calls.loc[1, "type"] = "gain" if calls.loc[1, "type"] == "loss" else "loss"                    # wrong direction
    bed = tmp_path / "calls.bed"; calls.to_csv(bed, sep="\t", index=False, header=False)
    (tmp_path / "tm.tsv").write_text(f"user\t{truth}\n")
    run(BIN / "score_calls.py", "--calls", bed, "--truth-map", tmp_path / "tm.tsv",
        "--fai", ASSETS / "GRCh37.autosomes.sizes", "--sample", "user", "--outdir", tmp_path)
    s = pd.read_csv(tmp_path / "event_scores.tsv", sep="\t").iloc[0]
    assert (s.TP, s.FP, s.FN) == (8, 2, 22)
    assert abs(s.F1 - 2 * 0.8 * (8 / 30) / (0.8 + 8 / 30)) < 1e-9


def test_collect_and_consensus(tmp_path):
    # three callers agree on a 3-Mb loss on chr1; one caller adds a gain on chr2
    files = {
        "a.calls.vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                       "1\t10000001\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=13000000\n",
        "b.call.cns": "chromosome\tstart\tend\tgene\tlog2\tcn\nchr1\t10000000\t13000000\t-\t-0.9\t1\n"
                      "chr2\t5000000\t9000000\t-\t0.5\t3\n",
        "c.seg.txt": "ID\tchrom\tstart\tend\tCorrected_Copy_Number\nS\t1\t10000000\t13000000\t1\n",
        "d.params.txt": "Tumor Fraction:\t0.031\nPloidy:\t2.1\n",
    }
    for n, txt in files.items(): (tmp_path / n).write_text(txt)
    man = ["sample_id\tbase\tlevel\tgroup\texpected_tf\ttool\tconfig\tkind\tfile",
           f"S\tS\tfull\tcancer\t\tQDNAseq\t100kb\tcalls\t{tmp_path/'a.calls.vcf'}",
           f"S\tS\tfull\tcancer\t\tCNVkit\t30kb\tcalls\t{tmp_path/'b.call.cns'}",
           f"S\tS\tfull\tcancer\t\tichorCNA\t1Mb low-TF\tcalls\t{tmp_path/'c.seg.txt'}",
           f"S\tS\tfull\tcancer\t\tichorCNA\t1Mb low-TF\ttf\t{tmp_path/'d.params.txt'}"]
    (tmp_path / "m.tsv").write_text("\n".join(man) + "\n")
    run(BIN / "collect_calls.py", "--manifest", tmp_path / "m.tsv", "--fai", ASSETS / "GRCh37.autosomes.sizes", "--outdir", tmp_path)
    calls = pd.read_csv(tmp_path / "all_calls.tsv", sep="\t", dtype={"chrom": str})
    assert len(calls) == 4 and set(calls.chrom) == {"1", "2"}
    tf = pd.read_csv(tmp_path / "tumour_fraction.tsv", sep="\t")
    assert abs(tf.TF.iloc[0] - 0.031) < 1e-9
    run(BIN / "consensus.py", "--calls", tmp_path / "all_calls.tsv", "--fai", ASSETS / "GRCh37.autosomes.sizes",
        "--primary", "QDNAseq:100kb,CNVkit:30kb,ichorCNA:1Mb low-TF", "--outdir", tmp_path)
    cons = pd.read_csv(tmp_path / "consensus_calls.tsv", sep="\t", dtype={"chrom": str})
    assert list(cons[["chrom", "start", "end", "type"]].itertuples(index=False, name=None)) == [("1", 10_000_000, 13_000_000, "loss")]
    ag = pd.read_csv(tmp_path / "caller_agreement.tsv", sep="\t").set_index("tool")
    assert ag.loc["QDNAseq", "F1"] == 1.0 and ag.loc["CNVkit", "precision"] < 1.0


def test_cnvmeann_train_and_predict(tmp_path):
    """CNV-MEANN network learns a strong simulated gain and calls it in an unlabelled sample."""
    import numpy as np
    rng = np.random.default_rng(1)

    def make(d, gain_chr, label=True):
        d.mkdir()
        for c in range(1, 23):
            n = 200; rd = rng.normal(1, 0.1, n); lab = np.zeros(n, int)
            if c == gain_chr: rd[50:90] *= 1.8; lab[50:90] = 1
            if c == gain_chr + 7: rd[100:150] *= 0.3; lab[100:150] = 2
            x = np.column_stack([np.arange(1, n + 1), rd, rng.uniform(0.35, 0.55, n), np.abs(rng.normal(0, 0.1, n)),
                                 np.ones(n), lab if label else np.zeros(n, int)])
            np.savetxt(d / f"{d.name}.chr{c}.bin500000.txt", x, delimiter="\t", fmt=["%d", "%.5f", "%.5f", "%.5f", "%.5f", "%d"])
    dirs = []
    for i in range(4):
        make(tmp_path / f"tr{i}", 3 + i); dirs.append(tmp_path / f"tr{i}")
    make(tmp_path / "test", 5, label=False)
    run(BIN / "cnvmeann_train.py", "--bin", 500000, "--out", tmp_path / "m.joblib", *dirs)
    run(BIN / "cnvmeann_predict.py", "--model", tmp_path / "m.joblib", "--features", tmp_path / "test", "--out", tmp_path / "ev.tsv")
    ev = pd.read_csv(tmp_path / "ev.tsv", sep="\t", dtype={"chrom": str})
    assert ((ev.chrom == "5") & (ev.type == "gain")).any()


def test_lyceum_exon_calls_to_segments(tmp_path):
    rows = ["Sample Name\tChromosome\tCNV Start Index\tCNV End Index\tLYCEUM Prediction",
            "LYC\tchr1\t1000\t1200\t<DEL>", "LYC\tchr1\t5000\t5200\t<DEL>", "LYC\tchr1\t9000\t9100\t<NO-CALL>",
            "LYC\tchr1\t9500\t9600\t<DEL>", "LYC\tchr2\t100\t300\t<DUP>"]
    (tmp_path / "c.csv").write_text("\n".join(rows) + "\n")
    run(BIN / "lyceum_tools.py", "calls", "--in", tmp_path / "c.csv", "--out", tmp_path / "s.tsv")
    s = pd.read_csv(tmp_path / "s.tsv", sep="\t", dtype={"chrom": str})
    assert list(s[["chrom", "start", "end", "type"]].itertuples(index=False, name=None)) == \
        [("1", 1000, 5200, "loss"), ("1", 9500, 9600, "loss"), ("2", 100, 300, "gain")]
