"""Offline tests for 01_primers.py. Each guard is shown to fire.

QIIME 2 is replaced by a fake runner, except in the parser test, which reads a real
q2-cutadapt report excerpt with interleaved QIIME output (fixtures/README.md). The
timeout test starts real processes.
"""

import csv
import importlib.util
import pathlib
import re
import sys
import time
import zipfile

import pytest

HERE = pathlib.Path(__file__).resolve().parent
FIX = HERE / "fixtures"
spec = importlib.util.spec_from_file_location("primers", HERE.parent / "01_primers.py")
pr = importlib.util.module_from_spec(spec)
sys.modules["primers"] = pr
spec.loader.exec_module(pr)


def make_qza(path, samples, drop_reverse=()):
    """A zip laid out like a QIIME 2 paired-end demux artifact, MANIFEST only."""
    lines = ["sample-id,filename,direction"]
    for i, s in enumerate(samples):
        lines.append(f"{s},{s}_{i}_L001_R1_001.fastq.gz,forward")
        if s not in drop_reverse:
            lines.append(f"{s},{s}_{i}_L001_R2_001.fastq.gz,reverse")
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("uuid/data/MANIFEST", "\n".join(lines) + "\n")
        zf.writestr("uuid/metadata.yaml", "uuid: x\n")
    return str(path)


def block(sample, idx, pairs_in, r1=0, r2=0, pairs_out=None, prefix=""):
    """One report in the layout cutadapt 5.1 prints inside q2-cutadapt."""
    pairs_out = pairs_in if pairs_out is None else pairs_out
    return (f"{prefix}This is cutadapt 5.1 with Python 3.10.14\n"
            f"Command line parameters: -o /tmp/out/{sample}_{idx}_L001_R1_001.fastq.gz "
            f"-p /tmp/out/{sample}_{idx}_L001_R2_001.fastq.gz --front ^AAA -G ^CCC "
            f"/tmp/in/{sample}_{idx}_L001_R1_001.fastq.gz /tmp/in/{sample}_{idx}_L001_R2_001.fastq.gz\n"
            "Processing paired-end reads on 4 cores ...\n\n=== Summary ===\n\n"
            f"Total read pairs processed:          {pairs_in:,}\n"
            f"  Read 1 with adapter:                 {r1:,} (0.0%)\n"
            f"  Read 2 with adapter:                 {r2:,} (0.0%)\n\n"
            "== Read fate breakdown ==\n"
            f"Pairs written (passing filters):     {pairs_out:,} (100.0%)\n\n")


def mapping(samples):
    return {f"{s}_{i}_L001_R1_001.fastq.gz": s for i, s in enumerate(samples)}


class FakeRunner:
    """Stands in for pr.run_cmd. Writes trimmed.qza and returns cutadapt reports."""

    def __init__(self, samples, reports=None, rc=0, write_output=True, out_samples=None):
        self.samples = samples
        self.reports = reports
        self.rc = rc
        self.write_output = write_output
        self.out_samples = out_samples or samples
        self.calls = []

    def __call__(self, cmd, timeout):
        self.calls.append(cmd)
        assert timeout > 0
        if self.write_output:
            make_qza(cmd[cmd.index("--o-trimmed-sequences") + 1], self.out_samples)
        reports = self.reports
        if reports is None:
            reports = "".join(block(s, i, 1000) for i, s in enumerate(self.samples))
        return self.rc, reports, "Running external command line application(s).\n"


def run_main(monkeypatch, tmp_path, runner, samples, extra=()):
    qza = make_qza(tmp_path / "demux.qza", samples)
    monkeypatch.setattr(pr, "run_cmd", runner)
    return pr.main(["-i", qza, "-o", str(tmp_path / "out"), *extra])


# Input guards

def test_bad_primer_is_fatal():
    with pytest.raises(pr.PrimerError, match="IUPAC"):
        pr.validate_primer("GTGCCAGC-XX")


def test_short_primer_is_fatal():
    with pytest.raises(pr.PrimerError, match="at least 10"):
        pr.validate_primer("ACGT")


def test_missing_input_is_fatal(tmp_path):
    with pytest.raises(pr.PrimerError, match="not found"):
        pr.read_manifest(str(tmp_path / "nope.qza"))


def test_non_artifact_is_fatal(tmp_path):
    (tmp_path / "x.qza").write_bytes(b"not a zip")
    with pytest.raises(pr.PrimerError, match="not a QIIME 2 artifact"):
        pr.read_manifest(str(tmp_path / "x.qza"))


def test_artifact_without_manifest_is_fatal(tmp_path):
    with zipfile.ZipFile(tmp_path / "x.qza", "w") as zf:
        zf.writestr("uuid/data/table.biom", "x")
    with pytest.raises(pr.PrimerError, match="MANIFEST"):
        pr.read_manifest(str(tmp_path / "x.qza"))


def test_sample_without_reverse_is_fatal(tmp_path):
    qza = make_qza(tmp_path / "d.qza", ["S1", "S2"], drop_reverse={"S2"})
    with pytest.raises(pr.PrimerError, match="lack a forward or reverse"):
        pr.read_manifest(qza)


# Command

def test_primers_are_anchored_and_untrimmed_kept(monkeypatch, tmp_path):
    runner = FakeRunner(["S1"])
    assert run_main(monkeypatch, tmp_path, runner, ["S1"]) == 0
    cmd = runner.calls[0]
    assert cmd[cmd.index("--p-front-f") + 1] == "^" + pr.PRIMER_F
    assert cmd[cmd.index("--p-front-r") + 1] == "^" + pr.PRIMER_R
    assert "--p-no-discard-untrimmed" in cmd and "--p-discard-untrimmed" not in cmd


def test_timeout_kills_the_whole_process_group():
    start = time.monotonic()
    with pytest.raises(pr.PrimerError, match="timed out"):
        pr.run_cmd(["bash", "-c", "sleep 30 & sleep 30; wait"], timeout=1)
    assert time.monotonic() - start < 10


# Report guards

def test_real_report_parses_including_interleaved_output():
    text = (FIX / "cutadapt_excerpt.log").read_text()
    # Build the MANIFEST mapping from the R1 inputs named in the excerpt.
    r1s = re.findall(r"(\S+_L001_R1_001\.fastq\.gz) \S+_L001_R2_001\.fastq\.gz\n", text)
    r1_to_sample = {pathlib.Path(p).name: pathlib.Path(p).name.split("_")[0] for p in r1s}
    assert "Command: This is cutadapt" in text  # the interleaved case is in the fixture
    rows, version = pr.parse_reports(text, r1_to_sample)
    assert version == "5.1"
    assert len(rows) == text.count("Total read pairs processed") == 14
    by = {r["sample"]: r for r in rows}
    assert by["SRR2143538"] == {"sample": "SRR2143538", "pairs_in": 42883,
                                "r1_with_primer": 0, "r2_with_primer": 0, "pairs_out": 42883}
    assert pr.interpret(rows).startswith("primers absent")


def test_missing_report_is_fatal():
    text = block("S1", 0, 100)
    with pytest.raises(pr.PrimerError, match="1 of 2 samples have no cutadapt report"):
        pr.parse_reports(text, mapping(["S1", "S2"]))


def test_duplicate_report_is_fatal():
    text = block("S1", 0, 100) * 2
    with pytest.raises(pr.PrimerError, match="two cutadapt reports"):
        pr.parse_reports(text, mapping(["S1"]))


def test_report_for_unknown_file_is_fatal():
    with pytest.raises(pr.PrimerError, match="not in the MANIFEST"):
        pr.parse_reports(block("S9", 0, 100), mapping(["S1"]))


def test_report_missing_a_count_is_fatal():
    text = block("S1", 0, 100).replace("Pairs written", "Pairs kept")
    with pytest.raises(pr.PrimerError, match="pairs_out"):
        pr.parse_reports(text, mapping(["S1"]))


def test_no_reports_is_fatal():
    with pytest.raises(pr.PrimerError, match="no cutadapt reports"):
        pr.parse_reports("Running external command line application(s).\n", mapping(["S1"]))


def test_lost_reads_without_discard_is_fatal(monkeypatch, tmp_path, caplog):
    reports = block("S1", 0, 1000) + block("S2", 1, 1000, pairs_out=990)
    runner = FakeRunner(["S1", "S2"], reports=reports)
    assert run_main(monkeypatch, tmp_path, runner, ["S1", "S2"]) == 1
    assert "S2: 1000 in, 990 out" in caplog.text


def test_loss_is_allowed_with_discard(monkeypatch, tmp_path):
    reports = block("S1", 0, 1000, r1=980, r2=975, pairs_out=970)
    runner = FakeRunner(["S1"], reports=reports)
    assert run_main(monkeypatch, tmp_path, runner, ["S1"], ["--discard-untrimmed"]) == 0
    assert "--p-discard-untrimmed" in runner.calls[0]


def test_cutadapt_failure_is_fatal(monkeypatch, tmp_path, caplog):
    assert run_main(monkeypatch, tmp_path, FakeRunner(["S1"], rc=1), ["S1"]) == 1
    assert "exited with code 1" in caplog.text


def test_missing_output_is_fatal(monkeypatch, tmp_path, caplog):
    runner = FakeRunner(["S1"], write_output=False)
    assert run_main(monkeypatch, tmp_path, runner, ["S1"]) == 1
    assert "wrote no" in caplog.text


def test_output_with_different_samples_is_fatal(monkeypatch, tmp_path, caplog):
    runner = FakeRunner(["S1", "S2"], out_samples=["S1"])
    assert run_main(monkeypatch, tmp_path, runner, ["S1", "S2"]) == 1
    assert "different sample set" in caplog.text


# Results

def test_interpret_mixed_is_flagged():
    rows = [{"sample": "S1", "pairs_in": 100, "r1_with_primer": 0, "r2_with_primer": 0},
            {"sample": "S2", "pairs_in": 100, "r1_with_primer": 95, "r2_with_primer": 96}]
    assert pr.interpret(rows).startswith("MIXED")


def test_interpret_present():
    rows = [{"sample": "S1", "pairs_in": 100, "r1_with_primer": 97, "r2_with_primer": 96}]
    assert pr.interpret(rows).startswith("primers present")


def test_happy_path_writes_table(monkeypatch, tmp_path, caplog):
    reports = block("S2", 1, 500, prefix="Command: ") + block("S1", 0, 1000, r1=3)
    runner = FakeRunner(["S1", "S2"], reports=reports)
    assert run_main(monkeypatch, tmp_path, runner, ["S1", "S2"]) == 0
    out = tmp_path / "out"
    rows = list(csv.DictReader(open(out / "primer_summary.tsv"), delimiter="\t"))
    assert [(r["sample"], r["pairs_in"], r["r1_with_primer"], r["pct_pairs_kept"])
            for r in rows] == [("S1", "1000", "3", "100.00"), ("S2", "500", "0", "100.00")]
    assert (out / "cutadapt_report.log").is_file()
    assert (out / "primers_log.txt").is_file()
    assert "primers absent" in caplog.text
