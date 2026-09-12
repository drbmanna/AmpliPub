"""Offline tests for 00_demux.py. Each guard is shown to fire.

QIIME 2 is replaced by a fake runner that writes the artifacts emp-paired and
tabulate-read-counts would write. The orientation check is the important one: it is the
guard that stands between a wrong flag and a near-empty feature table.
"""

import csv
import gzip
import importlib.util
import pathlib
import sys
import time
import zipfile

import pytest

HERE = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("demux", HERE.parent / "00_demux.py")
d = importlib.util.module_from_spec(spec)
sys.modules["demux"] = d
spec.loader.exec_module(d)

BARCODES = {"s1": "AAAACCCCGGGG", "s2": "TTTTCCCCGGGG", "s3": "GGGGCCCCAAAA"}


def emp_dir(tmp_path, reads=100, files=d.EMP_FILES):
    p = tmp_path / "emp"
    p.mkdir(exist_ok=True)
    for name in files:
        with gzip.open(p / name, "wt") as fh:
            for i in range(reads):
                fh.write(f"@r{i}\nACGT\n+\nIIII\n")
    return str(p)


def metadata(tmp_path, barcodes=None, column="barcode-sequence", crlf=False,
             name="md.tsv"):
    barcodes = BARCODES if barcodes is None else barcodes
    lines = [f"sample-id\t{column}\tdx"]
    for s, b in barcodes.items():
        lines.append(f"{s}\t{b}\tcase")
    text = ("\r\n" if crlf else "\n").join(lines) + ("\r\n" if crlf else "\n")
    p = tmp_path / name
    p.write_bytes(text.encode())
    return str(p)


def counts_qzv(path, counts):
    lines = ["id\tforward sequence count\treverse sequence count",
             "#q2:types\tnumeric\tnumeric"]
    for s, n in counts.items():
        lines.append(f"{s}\t{n}\t{n}")
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("uuid/data/metadata.tsv", "\n".join(lines) + "\n")
    return str(path)


def details_qza(path, records=3, corrected=1, orphan=1):
    lines = ["id\tsample\tbarcode\terrors", "#q2:types\tcategorical\tcategorical\tnumeric"]
    for i in range(records):
        if i < orphan:
            lines.append(f"r{i}\t\tAAAACCCCGGGG\t0")
        elif i < orphan + corrected:
            lines.append(f"r{i}\ts1\tAAAACCCCGGGG\t1")
        else:
            lines.append(f"r{i}\ts1\tAAAACCCCGGGG\t0")
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("uuid/data/metadata.tsv", "\n".join(lines) + "\n")
    return str(path)


# ---- the barcode column --------------------------------------------------

def test_barcodes_are_read(tmp_path):
    assert d.read_barcodes(metadata(tmp_path), "barcode-sequence") == BARCODES


def test_cr_line_endings_are_handled(tmp_path):
    assert d.read_barcodes(metadata(tmp_path, crlf=True), "barcode-sequence") == BARCODES


def test_a_missing_column_names_the_columns(tmp_path):
    with pytest.raises(d.DemuxError, match="no column 'barcode'"):
        d.read_barcodes(metadata(tmp_path), "barcode")


def test_two_samples_sharing_a_barcode_fires(tmp_path):
    md = metadata(tmp_path, {"s1": "AAAACCCCGGGG", "s2": "AAAACCCCGGGG"})
    with pytest.raises(d.DemuxError, match="share the barcode"):
        d.read_barcodes(md, "barcode-sequence")


def test_mixed_barcode_lengths_fire(tmp_path):
    md = metadata(tmp_path, {"s1": "AAAACCCCGGGG", "s2": "TTTTCCCC"})
    with pytest.raises(d.DemuxError, match="mixed lengths"):
        d.read_barcodes(md, "barcode-sequence")


def test_a_non_dna_column_fires(tmp_path):
    md = metadata(tmp_path, {"s1": "normalcase12", "s2": "adenomacas12"})
    with pytest.raises(d.DemuxError, match="not a DNA sequence"):
        d.read_barcodes(md, "barcode-sequence")


def test_an_empty_barcode_fires(tmp_path):
    md = metadata(tmp_path, {"s1": "AAAACCCCGGGG", "s2": ""})
    with pytest.raises(d.DemuxError, match="has no barcode"):
        d.read_barcodes(md, "barcode-sequence")


# ---- the input directory -------------------------------------------------

def test_a_missing_emp_file_says_which_and_suggests_the_manifest(tmp_path):
    partial = emp_dir(tmp_path, files=("forward.fastq.gz", "reverse.fastq.gz"))
    with pytest.raises(d.DemuxError, match="barcodes.fastq.gz"):
        d.check_input_dir(partial)
    with pytest.raises(d.DemuxError, match="import with a manifest instead"):
        d.check_input_dir(partial)


def test_an_empty_emp_file_fires(tmp_path):
    p = emp_dir(tmp_path)
    (pathlib.Path(p) / "barcodes.fastq.gz").write_bytes(b"")
    with pytest.raises(d.DemuxError, match="is empty"):
        d.check_input_dir(p)


def test_reads_are_counted(tmp_path):
    assert d.count_reads(str(pathlib.Path(emp_dir(tmp_path, reads=42)) / "barcodes.fastq.gz")) == 42


def test_a_truncated_fastq_fires(tmp_path):
    p = tmp_path / "bad.fastq.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("@r1\nACGT\n+\n")
    with pytest.raises(d.DemuxError, match="not a whole number of FASTQ records"):
        d.count_reads(str(p))


# ---- parsing -------------------------------------------------------------

def test_per_sample_counts_skip_the_types_row(tmp_path):
    path = counts_qzv(tmp_path / "c.qzv", {"s1": 10, "s2": 20})
    got = d.parse_per_sample_counts(d.read_tsv_from_artifact(path, "/data/metadata.tsv"),
                                    "c")
    assert got == {"s1": 10, "s2": 20}


def test_counts_without_a_count_column_fire():
    with pytest.raises(d.DemuxError, match="no read count column"):
        d.parse_per_sample_counts("id\tsomething\ns1\t5\n", "x")


def test_error_corrections_are_summarised(tmp_path):
    path = details_qza(tmp_path / "e.qza", records=5, corrected=2, orphan=1)
    got = d.parse_error_corrections(d.read_tsv_from_artifact(path, "/data/metadata.tsv"))
    assert got == {"records": 5, "corrected": 2, "uncorrectable": 1}


# ---- end to end ----------------------------------------------------------

class FakeRunner:
    """Stands in for d.run_cmd. `assign` maps (rev_barcodes, rev_mapping) to counts."""

    def __init__(self, assign=None, rc=0, write=True):
        self.assign = assign if assign is not None else {
            (False, False): {"s1": 30, "s2": 30, "s3": 30}}
        self.rc, self.write = rc, write
        self.calls = []

    def __call__(self, cmd, timeout):
        self.calls.append(cmd)
        assert timeout > 0
        if not self.write:
            return self.rc, "", "Plugin error from demux\n"
        if "import" in cmd:
            pathlib.Path(cmd[cmd.index("--output-path") + 1]).write_bytes(b"emp")
        elif "emp-paired" in cmd:
            rb = "--p-rev-comp-barcodes" in cmd
            rm = "--p-rev-comp-mapping-barcodes" in cmd
            seqs = pathlib.Path(cmd[cmd.index("--o-per-sample-sequences") + 1])
            seqs.write_bytes(b"seqs")
            self._last = (rb, rm)
            details_qza(cmd[cmd.index("--o-error-correction-details") + 1])
        else:  # tabulate-read-counts
            seqs = pathlib.Path(cmd[cmd.index("--i-sequences") + 1])
            tag = seqs.stem.replace("per_sample_sequences", "")
            key = (("_rb1" in tag), ("_rm1" in tag)) if tag else (False, False)
            counts_qzv(cmd[cmd.index("--o-visualization") + 1],
                       self.assign.get(key, {}))
        return self.rc, "", ""


def run_main(tmp_path, runner, extra=(), reads=100, monkeypatch=None):
    monkeypatch.setattr(d, "run_cmd", runner)
    return d.main(["-i", emp_dir(tmp_path, reads=reads), "-m", metadata(tmp_path),
                   "-o", str(tmp_path / "out"), *extra])


def test_end_to_end_accounts_for_every_read(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(), monkeypatch=monkeypatch) == 0
    out = tmp_path / "out"
    rows = {r["sample-id"]: r for r in
            csv.DictReader((out / "demux_counts.tsv").open(), delimiter="\t")}
    assert rows["s1"]["reads"] == "30" and rows["s1"]["barcode"] == BARCODES["s1"]
    log = (out / "demux_log.txt").read_text()
    assert "assigned 90 of 100 read pairs (90.00%); 10 unassigned (10.00%)" in log
    assert "barcode error correction" in log
    assert "unassigned fraction is a measurement, not waste" in log


def test_flags_are_passed_explicitly_including_the_defaults(tmp_path, monkeypatch):
    runner = FakeRunner()
    assert run_main(tmp_path, runner, monkeypatch=monkeypatch) == 0
    cmd = [c for c in runner.calls if "emp-paired" in c][0]
    assert "--p-no-rev-comp-barcodes" in cmd
    assert "--p-no-rev-comp-mapping-barcodes" in cmd
    assert "--p-golay-error-correction" in cmd


def test_no_golay_is_honoured(tmp_path, monkeypatch):
    runner = FakeRunner()
    assert run_main(tmp_path, runner, extra=["--no-golay"], monkeypatch=monkeypatch) == 0
    cmd = [c for c in runner.calls if "emp-paired" in c][0]
    assert "--p-no-golay-error-correction" in cmd


def test_a_wrong_orientation_is_found_and_named(tmp_path, monkeypatch):
    """The guard this stage exists for: the given flags assign almost nothing, another
    orientation assigns nearly everything, and nothing downstream is written."""
    runner = FakeRunner(assign={(False, False): {"s1": 1, "s2": 1},
                                (False, True): {"s1": 45, "s2": 45}})
    assert run_main(tmp_path, runner, monkeypatch=monkeypatch) == 1
    log = (tmp_path / "out" / "demux_log.txt").read_text()
    assert "the barcode orientation looks wrong" in log
    assert "rev-comp-mapping-barcodes=True assigns 90 reads" in log
    assert "near-empty table is the failure this check exists to prevent" in log


def test_low_assignment_with_no_better_orientation_says_so(tmp_path, monkeypatch):
    runner = FakeRunner(assign={(False, False): {"s1": 2, "s2": 2}})
    assert run_main(tmp_path, runner, monkeypatch=monkeypatch) == 1
    log = (tmp_path / "out" / "demux_log.txt").read_text()
    assert "no barcode orientation does better" in log


def test_the_orientation_check_can_be_switched_off(tmp_path, monkeypatch):
    runner = FakeRunner(assign={(False, False): {"s1": 2, "s2": 2}})
    assert run_main(tmp_path, runner, extra=["--no-orientation-check"],
                    monkeypatch=monkeypatch) == 0
    assert len([c for c in runner.calls if "emp-paired" in c]) == 1


def test_more_assigned_than_input_is_fatal(tmp_path, monkeypatch):
    runner = FakeRunner(assign={(False, False): {"s1": 500}})
    assert run_main(tmp_path, runner, reads=100, monkeypatch=monkeypatch) == 1
    assert "accounting cannot be trusted" in (tmp_path / "out" / "demux_log.txt").read_text()


def test_a_sample_with_no_reads_is_named(tmp_path, monkeypatch):
    runner = FakeRunner(assign={(False, False): {"s1": 50, "s2": 40}})
    assert run_main(tmp_path, runner, monkeypatch=monkeypatch) == 0
    assert "got zero reads: s3" in (tmp_path / "out" / "demux_log.txt").read_text()


def test_non_golay_length_warns(tmp_path, monkeypatch):
    md = metadata(tmp_path, {"s1": "AAAACCCC", "s2": "TTTTCCCC"}, name="short.tsv")
    monkeypatch.setattr(d, "run_cmd", FakeRunner(assign={(False, False): {"s1": 50,
                                                                          "s2": 40}}))
    assert d.main(["-i", emp_dir(tmp_path), "-m", md, "-o", str(tmp_path / "out")]) == 0
    log = (tmp_path / "out" / "demux_log.txt").read_text()
    assert "Golay error correction is on but the barcodes are 8 nt" in log


def test_a_qiime_failure_is_reported(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(rc=1, write=False), monkeypatch=monkeypatch) == 1
    assert "Plugin error from demux" in (tmp_path / "out" / "demux_log.txt").read_text()


@pytest.mark.parametrize("bad", [["--min-assigned", "0"], ["--min-assigned", "1.5"]])
def test_out_of_range_options_fire(tmp_path, monkeypatch, bad):
    assert run_main(tmp_path, FakeRunner(), extra=bad, monkeypatch=monkeypatch) == 1


# ---- the no-hang rule ----------------------------------------------------

def test_a_stalled_command_is_killed():
    start = time.monotonic()
    with pytest.raises(d.DemuxError, match="timed out after 1 s"):
        d.run_cmd(["sh", "-c", "sleep 30"], timeout=1)
    assert time.monotonic() - start < 10
