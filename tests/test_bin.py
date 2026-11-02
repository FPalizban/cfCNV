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
