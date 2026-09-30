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


# The retention floor. table() above holds the count row constant, which is what untrimmed
# fixed-length Illumina reads look like, so these build their own declining count rows.

def table_with_counts(medians, counts):
    """A seven-number table whose count row declines, as variable-length reads do."""
    head = "\t" + "\t".join(str(i) for i in range(1, len(medians) + 1))
    rows = [head, "count\t" + "\t".join(str(float(c)) for c in counts)]
    for label in ("2%", "9%", "25%", "50%", "75%", "91%", "98%"):
        rows.append(label + "\t" + "\t".join(str(float(m)) for m in medians))
    return "\n".join(rows) + "\n"


def short_tail_table():
    """251 positions, all medians well above Q30, but only 40% of reads reach the end.

    Medians stay high so truncation lands at 251 and the 285 bp overlap floor is cleared,
    which leaves retention as the only thing that can fail.
    """
    return table_with_counts([35] * 251, [10000.0] * 200 + [4000.0] * 51)


def test_check_reach_reports_the_surviving_fraction():
    assert q.check_reach([("R1", [100.0, 100.0, 90.0], 3)], 0.75) == {"R1": 0.9}


def test_retention_floor_fires_when_truncation_discards_too_much():
    with pytest.raises(q.QualityError, match="truncation discards too much"):
        q.check_reach([("R2", [100.0, 100.0, 50.0], 3)], 0.75)


def test_retention_floor_boundary_is_inclusive():
    assert q.check_reach([("R1", [100.0, 75.0], 2)], 0.75) == {"R1": 0.75}
    with pytest.raises(q.QualityError):
        q.check_reach([("R1", [100.0, 74.0], 2)], 0.75)


def test_retention_floor_names_every_failing_read_not_just_the_first():
    with pytest.raises(q.QualityError) as exc:
        q.check_reach([("R1", [100.0, 10.0], 2), ("R2", [100.0, 20.0], 2)], 0.75)
    assert "R1 keeps 10.0% at position 2" in str(exc.value)
    assert "R2 keeps 20.0% at position 2" in str(exc.value)


def test_the_cap_saves_a_run_the_floor_would_otherwise_refuse(monkeypatch, tmp_path, caplog):
    """Capping is the mechanism; check_reach is only the backstop behind it.

    With counts non-increasing the cap always lands somewhere the floor is satisfied, so
    check_reach cannot fire through this path. It stays for malformed profiles, and its
    own refusal is proven directly in the unit tests above.
    """
    t = short_tail_table()
    assert run_main(monkeypatch, tmp_path, FakeRunner(fwd=t, rev=t)) == 0
    assert "truncation capped at 200" in caplog.text
    assert "40.0% of reads reach 251" in caplog.text
    values = dict(csv.reader(open(tmp_path / "out" / "trunc_len.tsv"), delimiter="\t"))
    assert values["trunc_len_f"] == "200"
    assert values["trunc_len_f_from_quality"] == "251"


def test_a_lower_floor_lets_the_same_run_through(monkeypatch, tmp_path):
    t = short_tail_table()
    assert run_main(monkeypatch, tmp_path, FakeRunner(fwd=t, rev=t),
                    extra=("--min-reach", "0.3")) == 0


def test_trunc_len_records_the_floor_and_what_was_retained(monkeypatch, tmp_path):
    assert run_main(monkeypatch, tmp_path, FakeRunner()) == 0
    values = dict(csv.reader(open(tmp_path / "out" / "trunc_len.tsv"), delimiter="\t"))
    assert values["min_reach"] == "0.75"
    assert (values["reach_f"], values["reach_r"]) == ("1.0", "1.0")


def test_retention_cap_finds_the_last_position_the_reads_reach():
    assert q.retention_cap([100.0, 100.0, 80.0, 10.0], 0.75) == 3
    assert q.retention_cap([100.0, 100.0, 100.0], 0.75) == 3


def test_retention_cap_is_zero_when_nothing_qualifies():
    assert q.retention_cap([0.0, 0.0], 0.75) == 0
    assert q.retention_cap([], 0.75) == 0


def test_the_cap_only_binds_when_the_quality_rule_overshoots():
    counts = [100.0, 100.0, 80.0, 10.0]
    assert q.apply_retention_cap("R1", 4, 3, counts, 0.75) == 3   # binds
    assert q.apply_retention_cap("R1", 2, 3, counts, 0.75) == 2   # does not


def test_the_v4o_case_end_to_end(monkeypatch, tmp_path, caplog):
    """The real shape that exposed this: a long tail of untrimmed reads.

    Numbers from PRJNA643648 V4O. Median quality stays above Q30 to position 281 because
    past 273 it is the median of a small minority, while 90% of reads stop around 274.
    """
    medians = [38.0] * 281 + [22.0] * 20          # first dip at 282 -> quality picks 281
    counts = [10000.0] * 273 + [1000.0] * 28      # 10% past 273
    t = table_with_counts(medians, counts)
    assert run_main(monkeypatch, tmp_path, FakeRunner(fwd=t, rev=t)) == 0
    assert "truncation capped at 273, down from the 281" in caplog.text
    values = dict(csv.reader(open(tmp_path / "out" / "trunc_len.tsv"), delimiter="\t"))
    assert values["trunc_len_f"] == "273"                 # what dada2 will read
    assert values["trunc_len_f_from_quality"] == "281"    # what quality alone proposed
    assert values["retention_cap_f"] == "273"


def test_the_cap_is_applied_before_the_overlap_check(monkeypatch, tmp_path, caplog):
    """Overlap must be judged on the lengths that will really be used.

    Capped to 150 + 150 = 300, which clears the 285 bp floor. Uncapped it would have been
    260 + 260 and the overlap check would have passed on lengths the reads do not have.
    """
    medians = [38.0] * 260 + [10.0] * 40
    counts = [10000.0] * 150 + [1000.0] * 150
    t = table_with_counts(medians, counts)
    assert run_main(monkeypatch, tmp_path, FakeRunner(fwd=t, rev=t)) == 0
    assert "capped at 150" in caplog.text
    values = dict(csv.reader(open(tmp_path / "out" / "trunc_len.tsv"), delimiter="\t"))
    assert (values["trunc_len_f"], values["trunc_len_r"]) == ("150", "150")
    assert values["expected_overlap"] == str(150 + 150 - 253)


def test_a_capped_run_that_still_cannot_merge_is_refused(monkeypatch, tmp_path, caplog):
    """Capping must not rescue a run that genuinely has too little overlap."""
    medians = [38.0] * 300
    counts = [10000.0] * 100 + [1000.0] * 200      # cap 100 each, 200 < 285
    t = table_with_counts(medians, counts)
    assert run_main(monkeypatch, tmp_path, FakeRunner(fwd=t, rev=t)) == 1
    assert "below the 285 bp needed" in caplog.text


def test_the_dev_fixture_loses_no_reads_to_truncation():
    """Why 0.75 is safe as a default: the validated profiles keep every read.

    These are fixed-length Illumina reads, so both cutoffs retain 100%. The floor is a net
    for variable-length input, not a threshold the usual case sits near.
    """
    assert q.reach(q.parse_summary(FWD, "fwd")["count"], 251) == 100.0
    assert q.reach(q.parse_summary(REV, "rev")["count"], 218) == 100.0
