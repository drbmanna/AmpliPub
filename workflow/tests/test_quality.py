"""Offline tests for 02_quality.py. Each guard is shown to fire.

The real-profile tests read the forward and reverse quality tables from
qiime demux summarize on the 45-run dev subset (fixtures/README.md). QIIME 2 is replaced
by a fake runner. The timeout test starts real processes.
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
spec = importlib.util.spec_from_file_location("quality", HERE.parent / "02_quality.py")
q = importlib.util.module_from_spec(spec)
sys.modules["quality"] = q
spec.loader.exec_module(q)

FWD = (FIX / "dev_forward-seven-number-summaries.tsv").read_text()
REV = (FIX / "dev_reverse-seven-number-summaries.tsv").read_text()


def table(medians, count=100):
    """A seven-number table with the given medians, other rows copied from them."""
    n = len(medians)
    head = "\t" + "\t".join(str(i) for i in range(1, n + 1))
    rows = [head, "count\t" + "\t".join(str(float(count)) for _ in medians)]
    for label in ("2%", "9%", "25%", "50%", "75%", "91%", "98%"):
        rows.append(label + "\t" + "\t".join(str(float(m)) for m in medians))
    return "\n".join(rows) + "\n"


def make_qzv(path, fwd=FWD, rev=REV):
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("uuid/data/per-sample-fastq-counts.tsv", "sample\tcount\n")
        if fwd is not None:
            zf.writestr("uuid/data/forward-seven-number-summaries.tsv", fwd)
        if rev is not None:
            zf.writestr("uuid/data/reverse-seven-number-summaries.tsv", rev)
    return str(path)


class FakeRunner:
    """Stands in for q.run_cmd. Writes quality.qzv from the given tables."""

    def __init__(self, fwd=FWD, rev=REV, rc=0, write=True):
        self.fwd, self.rev, self.rc, self.write = fwd, rev, rc, write
        self.calls = []

    def __call__(self, cmd, timeout):
        self.calls.append(cmd)
        assert timeout > 0
        if self.write:
            make_qzv(cmd[cmd.index("--o-visualization") + 1], self.fwd, self.rev)
        return self.rc, "", "Plugin error from demux\n" if self.rc else ""


def run_main(monkeypatch, tmp_path, runner, extra=()):
    qza = tmp_path / "trimmed.qza"
    qza.write_bytes(b"x")
    monkeypatch.setattr(q, "run_cmd", runner)
    return q.main(["-i", str(qza), "-o", str(tmp_path / "out"), *extra])


# The rule on real data

def test_real_profiles_give_the_pre_registered_lengths():
    fwd = q.parse_summary(FWD, "fwd")
    rev = q.parse_summary(REV, "rev")
    assert len(fwd["50%"]) == len(rev["50%"]) == 251
    # Checked by eye in the table: no R1 median falls below 30; the first R2 median
    # below 30 is at position 219.
    assert q.trunc_len(fwd["50%"], 30, "R1") == 251
    assert rev["50%"][218] < 30 and all(m >= 30 for m in rev["50%"][:218])
    assert q.trunc_len(rev["50%"], 30, "R2") == 218


def test_rule_stops_at_the_first_dip_even_if_quality_recovers():
    assert q.trunc_len([35, 35, 29, 35, 35], 30, "R2") == 2


# Guards

def test_overlap_floor_fires_below_285():
    with pytest.raises(q.QualityError, match="would not merge"):
        q.check_overlap(150, 134, 253, 12, 20)


def test_overlap_floor_boundary_passes():
    assert q.check_overlap(150, 135, 253, 12, 20) == 32


def test_bad_first_position_is_fatal():
    with pytest.raises(q.QualityError, match="position 1"):
        q.trunc_len([20, 35, 35], 30, "R1")


def test_positions_out_of_order_are_fatal():
    bad = table([35, 35, 35]).replace("\t1\t2\t3", "\t1\t3\t2", 1)
    with pytest.raises(q.QualityError, match="not 1..N"):
        q.parse_summary(bad, "x")


def test_missing_median_row_is_fatal():
    bad = "\n".join(ln for ln in table([35, 35]).splitlines() if not ln.startswith("50%"))
    with pytest.raises(q.QualityError, match="'50%'"):
        q.parse_summary(bad, "x")


def test_ragged_row_is_fatal():
    bad = table([35, 35, 35]).replace("50%\t35.0\t35.0\t35.0", "50%\t35.0\t35.0")
    with pytest.raises(q.QualityError, match="has 2 values for 3 positions"):
        q.parse_summary(bad, "x")


def test_single_end_input_is_fatal(tmp_path):
    with pytest.raises(q.QualityError, match="paired-end"):
        q.read_profiles(make_qzv(tmp_path / "v.qzv", rev=None))


def test_missing_input_is_fatal(tmp_path):
    assert q.main(["-i", str(tmp_path / "nope.qza"), "-o", str(tmp_path / "out")]) == 1


def test_summarize_failure_is_fatal(monkeypatch, tmp_path, caplog):
    assert run_main(monkeypatch, tmp_path, FakeRunner(rc=1)) == 1
    assert "exited with code 1" in caplog.text


def test_reads_that_cannot_merge_stop_the_run(monkeypatch, tmp_path, caplog):
    early_drop = table([35] * 120 + [25] * 131)
    assert run_main(monkeypatch, tmp_path, FakeRunner(fwd=early_drop, rev=early_drop)) == 1
    assert "240 bp, below the 285 bp needed" in caplog.text
    assert not (tmp_path / "out" / "trunc_len.tsv").exists()


def test_timeout_kills_the_whole_process_group():
    start = time.monotonic()
    with pytest.raises(q.QualityError, match="timed out"):
        q.run_cmd(["bash", "-c", "sleep 30 & sleep 30; wait"], timeout=1)
    assert time.monotonic() - start < 10


# Results

def test_reach():
    assert q.reach([100.0, 100.0, 80.0], 3) == 80.0


def test_happy_path_on_real_profiles(monkeypatch, tmp_path, caplog):
    runner = FakeRunner()
    assert run_main(monkeypatch, tmp_path, runner) == 0
    cmd = runner.calls[0]
    assert cmd[4:7] == ["qiime", "demux", "summarize"]
    assert cmd[cmd.index("--p-n") + 1] == "10000"
    out = tmp_path / "out"
    values = dict(csv.reader(open(out / "trunc_len.tsv"), delimiter="\t"))
    assert (values["trunc_len_f"], values["trunc_len_r"], values["expected_overlap"]) == \
        ("251", "218", "216")
    profile = list(csv.DictReader(open(out / "quality_profile.tsv"), delimiter="\t"))
    assert len(profile) == 251
    assert "--p-trunc-len-f 251 --p-trunc-len-r 218" in caplog.text
