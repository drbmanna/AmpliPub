"""Offline tests for 00_qc_raw.py. Each guard is shown to fire.

FastQC and MultiQC are replaced by a fake runner that writes report zips, except in
the parser test, which reads a real FastQC 0.12.1 report (fixtures/README.md). The
timeout tests start real processes.
"""

import csv
import importlib.util
import pathlib
import sys
import time
import zipfile

import pytest

HERE = pathlib.Path(__file__).resolve().parent
FIX = HERE / "fixtures"
spec = importlib.util.spec_from_file_location("qc_raw", HERE.parent / "00_qc_raw.py")
qc = importlib.util.module_from_spec(spec)
sys.modules["qc_raw"] = qc
spec.loader.exec_module(qc)

MODULES = [
    "Basic Statistics", "Per base sequence quality", "Per sequence quality scores",
    "Per base sequence content", "Per sequence GC content", "Per base N content",
    "Sequence Length Distribution", "Sequence Duplication Levels",
    "Overrepresented sequences", "Adapter Content",
]


def make_zip(outdir, fastq_name, total=1000, statuses=None, filename=None):
    """Write a report zip laid out like FastQC 0.12.1 output."""
    statuses = statuses or {}
    s = qc.stem(fastq_name)
    listed = filename or fastq_name
    summary = "".join(f"{statuses.get(m, 'PASS')}\t{m}\t{listed}\n" for m in MODULES)
    data = (f"##FastQC\t0.12.1\n>>Basic Statistics\tpass\n#Measure\tValue\n"
            f"Filename\t{listed}\nTotal Sequences\t{total}\n>>END_MODULE\n")
    with zipfile.ZipFile(pathlib.Path(outdir) / f"{s}_fastqc.zip", "w") as zf:
        zf.writestr(f"{s}_fastqc/summary.txt", summary)
        zf.writestr(f"{s}_fastqc/fastqc_data.txt", data)


class FakeRunner:
    """Stands in for qc.run_cmd. Writes the files FastQC and MultiQC would write."""

    def __init__(self, totals=None, statuses=None, skip=(), fastqc_rc=0, multiqc_writes=True):
        self.totals = totals or {}
        self.statuses = statuses or {}
        self.skip = set(skip)
        self.fastqc_rc = fastqc_rc
        self.multiqc_writes = multiqc_writes
        self.calls = []

    def __call__(self, cmd, timeout):
        self.calls.append(cmd)
        assert timeout > 0
        tool = cmd[4]
        if "--version" in cmd:
            return 0, f"{tool} version x\n", ""
        outdir = cmd[cmd.index("-o") + 1]
        if tool == "fastqc":
            if self.fastqc_rc:
                return self.fastqc_rc, "", "java.lang.OutOfMemoryError\n"
            for f in cmd[cmd.index("-o") + 2:]:
                name = pathlib.Path(f).name
                if name not in self.skip:
                    make_zip(outdir, name, self.totals.get(name, 1000),
                             self.statuses.get(name))
            return 0, "", ""
        if tool == "multiqc":
            if self.multiqc_writes:
                # Skip conda's own "-n <env>" at index 2.
                (pathlib.Path(outdir) / cmd[cmd.index("-n", 5) + 1]).write_text("<html></html>")
            return 0, "", ""
        raise AssertionError(f"unexpected command {cmd}")


def fastq_dir(tmp_path, names=("S1_1.fastq.gz", "S1_2.fastq.gz", "S2_1.fastq.gz",
                               "S2_2.fastq.gz")):
    d = tmp_path / "fastq"
    d.mkdir()
    for n in names:
        (d / n).write_bytes(b"@r\nACGT\n+\nFFFF\n")
    return d


def run_main(monkeypatch, tmp_path, indir, runner):
    monkeypatch.setattr(qc, "run_cmd", runner)
    return qc.main(["-i", str(indir), "-o", str(tmp_path / "out")])


# Input guards

def test_missing_input_dir_is_fatal(tmp_path):
    with pytest.raises(qc.QCError, match="not found"):
        qc.find_fastqs(str(tmp_path / "nope"))


def test_zero_fastqs_is_fatal(tmp_path):
    (tmp_path / "notes.txt").write_text("x")
    with pytest.raises(qc.QCError, match="no FASTQ files"):
        qc.find_fastqs(str(tmp_path))


def test_zero_fastqs_never_calls_fastqc(monkeypatch, tmp_path):
    runner = FakeRunner()
    empty = tmp_path / "empty"
    empty.mkdir()
    assert run_main(monkeypatch, tmp_path, empty, runner) == 1
    assert runner.calls == []


def test_empty_file_is_fatal(tmp_path):
    d = fastq_dir(tmp_path)
    (d / "S3_1.fastq.gz").write_bytes(b"")
    with pytest.raises(qc.QCError, match="empty FASTQ"):
        qc.find_fastqs(str(d))


def test_colliding_report_names_are_fatal(tmp_path):
    d = fastq_dir(tmp_path, names=("S1_1.fastq", "S1_1.fastq.gz"))
    with pytest.raises(qc.QCError, match="same FastQC report"):
        qc.find_fastqs(str(d))


# Process guards, with real processes

def test_timeout_kills_the_whole_process_group():
    # The background sleep keeps the output pipe open. If only the parent were killed,
    # reading the pipe would block for 30 s.
    start = time.monotonic()
    with pytest.raises(qc.QCError, match="timed out"):
        qc.run_cmd(["bash", "-c", "sleep 30 & sleep 30; wait"], timeout=1)
    assert time.monotonic() - start < 10


def test_stdin_is_closed():
    # cat reading an inherited stdin could wait for ever; with stdin closed it ends at once.
    rc, out, _ = qc.run_cmd(["cat"], timeout=10)
    assert (rc, out) == (0, "")


def test_missing_executable_is_fatal():
    with pytest.raises(qc.QCError, match="cannot start"):
        qc.run_cmd(["no-such-tool-amplipub"], timeout=10)


# Output guards

def test_fastqc_failure_is_fatal(monkeypatch, tmp_path, caplog):
    assert run_main(monkeypatch, tmp_path, fastq_dir(tmp_path), FakeRunner(fastqc_rc=1)) == 1
    assert "FastQC exited with code 1" in caplog.text


def test_missing_report_is_fatal(monkeypatch, tmp_path, caplog):
    runner = FakeRunner(skip={"S2_2.fastq.gz"})
    assert run_main(monkeypatch, tmp_path, fastq_dir(tmp_path), runner) == 1
    assert "1 of 4 reports are missing" in caplog.text


def test_report_for_wrong_file_is_fatal(tmp_path):
    make_zip(tmp_path, "S1_1.fastq.gz", filename="other.fastq.gz")
    with pytest.raises(qc.QCError, match="describes other.fastq.gz"):
        qc.collect_reports([str(tmp_path / "S1_1.fastq.gz")], str(tmp_path))


def test_corrupt_zip_is_fatal(tmp_path):
    (tmp_path / "S1_1_fastqc.zip").write_bytes(b"not a zip")
    with pytest.raises(qc.QCError, match="not a readable zip"):
        qc.read_report(str(tmp_path / "S1_1_fastqc.zip"))


def test_pair_count_mismatch_is_fatal(monkeypatch, tmp_path, caplog):
    runner = FakeRunner(totals={"S2_2.fastq.gz": 999})
    assert run_main(monkeypatch, tmp_path, fastq_dir(tmp_path), runner) == 1
    assert "S2: 1000 vs 999" in caplog.text


def test_lone_mate_is_fatal():
    reports = [{"file": "S1_1.fastq.gz", "total_sequences": 5}]
    with pytest.raises(qc.QCError, match="only one mate"):
        qc.check_pairs(reports)


def test_qiime_export_names_are_paired():
    # Real names from `qiime tools export` of the dev demux.qza: the two mates of a
    # sample carry different file numbers.
    reports = [{"file": "SRR2143538_1_L001_R1_001.fastq.gz", "total_sequences": 5},
               {"file": "SRR2143538_46_L001_R2_001.fastq.gz", "total_sequences": 5},
               {"file": "SRR2143541_38_L001_R1_001.fastq.gz", "total_sequences": 7},
               {"file": "SRR2143541_83_L001_R2_001.fastq.gz", "total_sequences": 6}]
    with pytest.raises(qc.QCError, match="SRR2143541: 7 vs 6"):
        qc.check_pairs(reports)
    assert qc.check_pairs(reports[:2]) == 1


def test_single_end_files_are_not_paired():
    reports = [{"file": "S1.fastq.gz", "total_sequences": 5},
               {"file": "S2.fastq.gz", "total_sequences": 7}]
    assert qc.check_pairs(reports) == 0


def test_missing_multiqc_report_is_fatal(monkeypatch, tmp_path, caplog):
    runner = FakeRunner(multiqc_writes=False)
    assert run_main(monkeypatch, tmp_path, fastq_dir(tmp_path), runner) == 1
    assert "wrote no report" in caplog.text


# Results

def test_real_fastqc_report_parses():
    r = qc.read_report(str(FIX / "SRR2144126_1_fastqc.zip"))
    assert r["filename"] == "SRR2144126_1.fastq.gz"
    assert r["total_sequences"] == 41261
    assert len(r["modules"]) == 10
    assert r["modules"]["Per base sequence quality"] == "PASS"
    # The four failures on this real library are exactly the amplicon-expected ones.
    assert sorted(qc.flagged(r, expected=True)) == sorted(qc.EXPECTED_FOR_AMPLICONS)
    assert qc.flagged(r, expected=False) == []


def test_happy_path_writes_table_and_report(monkeypatch, tmp_path, caplog):
    runner = FakeRunner(statuses={
        "S1_1.fastq.gz": {"Per base sequence content": "FAIL"},
        "S2_2.fastq.gz": {"Adapter Content": "WARN", "Sequence Duplication Levels": "FAIL"},
    })
    assert run_main(monkeypatch, tmp_path, fastq_dir(tmp_path), runner) == 0
    out = tmp_path / "out"
    rows = list(csv.DictReader(open(out / "qc_summary.tsv"), delimiter="\t"))
    assert [r["file"] for r in rows] == ["S1_1.fastq.gz", "S1_2.fastq.gz",
                                         "S2_1.fastq.gz", "S2_2.fastq.gz"]
    by = {r["file"]: r for r in rows}
    assert by["S1_1.fastq.gz"]["to_check"] == ""
    assert by["S1_1.fastq.gz"]["expected_for_amplicons"] == "Per base sequence content"
    assert by["S2_2.fastq.gz"]["to_check"] == "Adapter Content"
    assert (out / "multiqc_report.html").is_file()
    assert (out / "qc_log.txt").is_file()
    assert "1 of 4 files have a module worth checking" in caplog.text
    fastqc_cmd = next(c for c in runner.calls if c[4] == "fastqc" and "--version" not in c)
    assert fastqc_cmd[:4] == ["conda", "run", "-n", "amplipub-qc"]
    assert "--noextract" in fastqc_cmd
