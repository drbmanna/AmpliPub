"""Offline tests for 07_collapse.py. Each guard is shown to fire.

The read-conservation check is the one that matters: pooling moves reads between columns
and must not create or destroy any. A fake runner stands in for QIIME 2 and can be told
to pool wrongly, so the guard is shown firing rather than assumed to work.
"""

import csv
import importlib.util
import pathlib
import sys
import time

import pytest

HERE = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("collapse", HERE.parent / "07_collapse.py")
c = importlib.util.module_from_spec(spec)
sys.modules["collapse"] = c
spec.loader.exec_module(c)

# r1 and r2 are two runs of sample A; r3 is sample B; r4 is sample C
RUNS = {"r1": "A", "r2": "A", "r3": "B", "r4": "C"}
COUNTS = {"r1": {"f1": 10, "f2": 5}, "r2": {"f1": 20, "f2": 0},
          "r3": {"f1": 0, "f2": 30}, "r4": {"f1": 7, "f2": 8}}


def biom_tsv(counts):
    """A biom convert --to-tsv table from {column: {feature: count}}."""
    cols = list(counts)
    features = sorted({f for v in counts.values() for f in v})
    lines = ["# Constructed from biom file", "#OTU ID\t" + "\t".join(cols)]
    for f in features:
        lines.append(f + "\t" + "\t".join(str(float(counts[col].get(f, 0)))
                                          for col in cols))
    return "\n".join(lines) + "\n"


def run_map_file(tmp_path, runs=None, column="sample_title", name="map.tsv"):
    runs = RUNS if runs is None else runs
    lines = [f"sample-id\tsample_accession\t{column}"]
    for run, sample in runs.items():
        lines.append(f"{run}\tSAMN{run}\t{sample}")
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\n")
    return str(p)


def metadata_file(tmp_path, ids=("A", "B"), name="md.tsv", crlf=True):
    lines = ["sample\tdx\tage"] + [f"{i}\tcase\t50" for i in ids]
    text = ("\r\n" if crlf else "\n").join(lines) + "\n"
    p = tmp_path / name
    p.write_bytes(text.encode())
    return str(p)


# ---- the run map ---------------------------------------------------------

def test_run_map_is_read(tmp_path):
    assert c.read_run_map(run_map_file(tmp_path), "sample_title") == RUNS


def test_a_missing_group_column_names_the_columns(tmp_path):
    with pytest.raises(c.CollapseError, match="no column 'nope'"):
        c.read_run_map(run_map_file(tmp_path), "nope")


def test_a_run_with_no_sample_fires(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("sample-id\tsample_title\nr1\tA\nr2\t\n")
    with pytest.raises(c.CollapseError, match="cannot be assigned to a sample"):
        c.read_run_map(str(p), "sample_title")


def test_a_duplicated_run_fires(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("sample-id\tsample_title\nr1\tA\nr1\tB\n")
    with pytest.raises(c.CollapseError, match="appears more than once"):
        c.read_run_map(str(p), "sample_title")


def test_cr_line_endings_are_handled(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_bytes(b"sample-id\tsample_title\rr1\tA\rr2\tB\r")
    assert c.read_run_map(str(p), "sample_title") == {"r1": "A", "r2": "B"}


# ---- id safety -----------------------------------------------------------

@pytest.mark.parametrize("bad", ["Mock community sample 2", "#A", " A", "A\tB", "",
                                 "sample-id", "id"])
def test_unusable_ids_are_found(bad):
    assert c.unsafe_ids(["good", bad]) == [bad]


def test_good_ids_pass():
    assert c.unsafe_ids(["A", "B_1", "2017660", "mock1"]) == []


def test_sanitize_collapses_whitespace_runs():
    assert c.sanitize("Mock community  sample 2") == "Mock_community_sample_2"


# ---- the exported table --------------------------------------------------

def test_table_totals_are_summed_per_column(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text(biom_tsv(COUNTS))
    assert c.read_table_totals(str(p)) == {"r1": 15.0, "r2": 20.0, "r3": 30.0, "r4": 15.0}


def test_a_wrong_first_column_fires(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("#OTU\ts1\nf1\t1.0\n")
    with pytest.raises(c.CollapseError, match="expected the feature id column"):
        c.read_table_totals(str(p))


def test_a_ragged_table_row_fires(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("#OTU ID\ts1\ts2\nf1\t1.0\n")
    with pytest.raises(c.CollapseError, match="values for 2 columns"):
        c.read_table_totals(str(p))


# ---- end to end ----------------------------------------------------------

class FakeRunner:
    """Stands in for c.run_cmd. `pool` decides how the grouped table comes back."""

    def __init__(self, pool=None, rc=0, write=True):
        self.pool = pool          # {sample: {feature: count}}, or None to pool correctly
        self.rc, self.write = rc, write
        self.calls = []
        self.grouping = None

    def _pooled(self):
        if self.pool is not None:
            return self.pool
        out: dict[str, dict[str, float]] = {}
        for run, sample in (self.grouping or RUNS).items():
            for feature, n in COUNTS.get(run, {}).items():
                out.setdefault(sample, {})
                out[sample][feature] = out[sample].get(feature, 0) + n
        return out

    def __call__(self, cmd, timeout):
        self.calls.append(cmd)
        assert timeout > 0
        if not self.write:
            return self.rc, "", "Plugin error from feature-table\n"
        if "group" in cmd:
            md = pathlib.Path(cmd[cmd.index("--m-metadata-file") + 1])
            rows = list(csv.DictReader(md.open(), delimiter="\t"))
            col = [k for k in rows[0] if k != "sample-id"][0]
            self.grouping = {r["sample-id"]: r[col] for r in rows}
            pathlib.Path(cmd[cmd.index("--o-grouped-table") + 1]).write_bytes(b"x")
        elif "export" in cmd:
            src = pathlib.Path(cmd[cmd.index("--input-path") + 1])
            dest = pathlib.Path(cmd[cmd.index("--output-path") + 1])
            dest.mkdir(parents=True, exist_ok=True)
            (dest / "feature-table.biom").write_bytes(b"biom")
            self._is_grouped = "by_sample" in src.name
        else:  # biom convert
            out = pathlib.Path(cmd[cmd.index("-o") + 1])
            counts = self._pooled() if self._is_grouped else COUNTS
            out.write_text(biom_tsv(counts))
        return self.rc, "", ""


def run_main(tmp_path, runner, extra=(), monkeypatch=None, runs=None):
    table = tmp_path / "table.qza"
    table.write_bytes(b"x")
    monkeypatch.setattr(c, "run_cmd", runner)
    return c.main(["-b", str(table), "-r", run_map_file(tmp_path, runs),
                   "-o", str(tmp_path / "out"), *extra])


def test_end_to_end_pools_and_proves_the_total(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(), monkeypatch=monkeypatch,
                    extra=["--expect-samples", "3"]) == 0
    out = tmp_path / "out"
    rows = {r["sample-id"]: r for r in
            csv.DictReader((out / "runs_per_sample.tsv").open(), delimiter="\t")}
    assert rows["A"]["n_runs"] == "2" and rows["A"]["runs"] == "r1;r2"
    assert rows["A"]["reads"] == "35"
    log = (out / "collapse_log.txt").read_text()
    assert "4 runs pooled into 3 samples; 1 sample(s) had more than one run" in log
    assert "80 reads in and 80 out, unchanged" in log


def test_read_loss_during_pooling_is_fatal(tmp_path, monkeypatch):
    """The guard the stage exists for: a mode that drops reads must not pass."""
    runner = FakeRunner(pool={"A": {"f1": 20, "f2": 5}, "B": {"f2": 30},
                              "C": {"f1": 7, "f2": 8}})
    assert run_main(tmp_path, runner, monkeypatch=monkeypatch) == 1
    log = (tmp_path / "out" / "collapse_log.txt").read_text()
    assert "pooling changed the read total: 80 before, 70 after" in log
    assert "must conserve reads" in log


def test_a_sample_pooled_to_the_wrong_total_is_fatal(tmp_path, monkeypatch):
    """Totals can match overall while an individual sample is wrong."""
    runner = FakeRunner(pool={"A": {"f1": 30, "f2": 10}, "B": {"f2": 25},
                              "C": {"f1": 7, "f2": 8}})
    assert run_main(tmp_path, runner, monkeypatch=monkeypatch) == 1
    assert "pooled to" in (tmp_path / "out" / "collapse_log.txt").read_text()


def test_a_run_missing_from_the_map_is_fatal(tmp_path, monkeypatch):
    runs = {"r1": "A", "r2": "A", "r3": "B"}          # r4 is in the table, not the map
    assert run_main(tmp_path, FakeRunner(), monkeypatch=monkeypatch, runs=runs) == 1
    log = (tmp_path / "out" / "collapse_log.txt").read_text()
    assert "not in the map" in log and "dropped without a word" in log


def test_a_run_missing_from_the_table_is_warned_and_ignored(tmp_path, monkeypatch):
    runs = dict(RUNS, r9="D")                          # r9 is in the map, not the table
    assert run_main(tmp_path, FakeRunner(), monkeypatch=monkeypatch, runs=runs) == 0
    assert "not in the table and are ignored: r9" in \
        (tmp_path / "out" / "collapse_log.txt").read_text()


def test_unsafe_ids_stop_the_run_by_default(tmp_path, monkeypatch):
    runs = dict(RUNS, r4="Mock community sample 2")
    assert run_main(tmp_path, FakeRunner(), monkeypatch=monkeypatch, runs=runs) == 1
    log = (tmp_path / "out" / "collapse_log.txt").read_text()
    assert "not usable as QIIME 2 sample ids" in log
    assert "--sanitize-ids" in log


def test_sanitizing_is_opt_in_and_recorded(tmp_path, monkeypatch):
    runs = dict(RUNS, r4="Mock community sample 2")
    assert run_main(tmp_path, FakeRunner(), monkeypatch=monkeypatch, runs=runs,
                    extra=["--sanitize-ids"]) == 0
    mapped = list(csv.DictReader((tmp_path / "out" / "sanitized_ids.tsv").open(),
                                 delimiter="\t"))
    assert mapped == [{"original_id": "Mock community sample 2",
                       "sanitized_id": "Mock_community_sample_2"}]


def test_sanitizing_into_a_collision_is_fatal(tmp_path, monkeypatch):
    runs = {"r1": "A_B", "r2": "A B", "r3": "C", "r4": "D"}
    assert run_main(tmp_path, FakeRunner(), monkeypatch=monkeypatch, runs=runs,
                    extra=["--sanitize-ids"]) == 1
    assert "would make" in (tmp_path / "out" / "collapse_log.txt").read_text()


def test_expect_samples_mismatch_is_fatal(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(), monkeypatch=monkeypatch,
                    extra=["--expect-samples", "4"]) == 1
    assert "expected 4 samples, pooling gave 3" in \
        (tmp_path / "out" / "collapse_log.txt").read_text()


# ---- the metadata join ---------------------------------------------------

def test_metadata_is_joined_and_gaps_reported(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(), monkeypatch=monkeypatch,
                    extra=["--metadata", metadata_file(tmp_path, ids=("A", "B", "Z"))]) == 0
    out = tmp_path / "out"
    rows = {r["sample-id"]: r for r in
            csv.DictReader((out / "sample_metadata.tsv").open(), delimiter="\t")}
    assert rows["A"]["in_study_metadata"] == "yes" and rows["A"]["dx"] == "case"
    assert rows["C"]["in_study_metadata"] == "no" and rows["C"]["dx"] == ""
    log = (out / "collapse_log.txt").read_text()
    assert "2 of 3 pooled samples matched" in log
    assert "no metadata row" in log and "match no sequenced sample: Z" in log
    assert "controls and mocks normally land in the unmatched list" in log


def test_metadata_that_matches_nothing_is_fatal(tmp_path, monkeypatch):
    md = metadata_file(tmp_path, ids=("X", "Y"))
    assert run_main(tmp_path, FakeRunner(), monkeypatch=monkeypatch,
                    extra=["--metadata", md]) == 1
    assert "not one pooled sample matches" in \
        (tmp_path / "out" / "collapse_log.txt").read_text()


def test_a_duplicated_metadata_id_is_fatal(tmp_path, monkeypatch):
    md = tmp_path / "dup.tsv"
    md.write_text("sample\tdx\nA\tcase\nA\tcontrol\n")
    assert run_main(tmp_path, FakeRunner(), monkeypatch=monkeypatch,
                    extra=["--metadata", str(md)]) == 1
    assert "appears more than once" in (tmp_path / "out" / "collapse_log.txt").read_text()


def test_a_qiime_failure_is_reported(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(rc=1, write=False), monkeypatch=monkeypatch) == 1
    assert "Plugin error from feature-table" in \
        (tmp_path / "out" / "collapse_log.txt").read_text()


# ---- the no-hang rule ----------------------------------------------------

def test_a_stalled_command_is_killed():
    start = time.monotonic()
    with pytest.raises(c.CollapseError, match="timed out after 1 s"):
        c.run_cmd(["sh", "-c", "sleep 30"], timeout=1)
    assert time.monotonic() - start < 10
