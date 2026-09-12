"""Offline tests for 08_filter.py. Each guard is shown to fire.

The important one: qiime feature-table filter-features drops empty samples by default,
so a feature filter can silently change the sample count. The fake runner can reproduce
that, so the guard is shown firing rather than assumed.
"""

import csv
import importlib.util
import pathlib
import sys
import time

import pytest

HERE = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("filterstage", HERE.parent / "08_filter.py")
f = importlib.util.module_from_spec(spec)
sys.modules["filterstage"] = f
spec.loader.exec_module(f)

# s1 and s2 are study samples; mock1 is a control
TABLES = {
    "00_input": {"s1": {"a": 100, "b": 50, "rare": 1}, "s2": {"a": 200, "b": 0, "rare": 0},
                 "mock1": {"a": 10, "b": 10, "rare": 0}},
    "01_taxa": {"s1": {"a": 100, "b": 50, "rare": 1}, "s2": {"a": 200, "rare": 0},
                "mock1": {"a": 10, "b": 10, "rare": 0}},
    "02_samples": {"s1": {"a": 100, "b": 50, "rare": 1}, "s2": {"a": 200, "rare": 0}},
    "03_prevalence": {"s1": {"a": 100, "b": 50}, "s2": {"a": 200}},
    "04_depth": {"s1": {"a": 100, "b": 50}, "s2": {"a": 200}},
}


def biom_tsv(counts):
    cols = list(counts)
    features = sorted({x for v in counts.values() for x in v})
    lines = ["# Constructed from biom file", "#OTU ID\t" + "\t".join(cols)]
    for feat in features:
        lines.append(feat + "\t" + "\t".join(str(float(counts[col].get(feat, 0)))
                                             for col in cols))
    return "\n".join(lines) + "\n"


def inputs(tmp_path):
    for name in ("table.qza", "rep_seqs.qza", "taxonomy.qza"):
        (tmp_path / name).write_bytes(b"x")
    md = tmp_path / "md.tsv"
    md.write_text("sample-id\tin_study_metadata\ns1\tyes\ns2\tyes\nmock1\tno\n")
    return str(md)


class FakeRunner:
    """Stands in for f.run_cmd, serving a table per export tag."""

    def __init__(self, tables=None, rc=0, write=True, fail_on=None):
        self.tables = dict(TABLES) if tables is None else tables
        self.rc, self.write, self.fail_on = rc, write, fail_on
        self.calls = []
        self._tag = "00_input"

    def __call__(self, cmd, timeout):
        self.calls.append(cmd)
        assert timeout > 0
        joined = " ".join(cmd)
        if self.fail_on and self.fail_on in joined:
            return 1, "", f"Plugin error from {self.fail_on}\n"
        if not self.write:
            return self.rc, "", "Plugin error\n"
        if "export" in cmd:
            dest = pathlib.Path(cmd[cmd.index("--output-path") + 1])
            dest.mkdir(parents=True, exist_ok=True)
            (dest / "feature-table.biom").write_bytes(b"biom")
            self._tag = dest.name.replace("export_", "")
        elif cmd[0] == "cp":
            pathlib.Path(cmd[2]).write_bytes(pathlib.Path(cmd[1]).read_bytes())
        elif "convert" in joined:
            table = self.tables.get(self._tag, {})
            pathlib.Path(cmd[cmd.index("-o") + 1]).write_text(biom_tsv(table))
        else:
            for flag in ("--o-filtered-table", "--o-filtered-data"):
                if flag in cmd:
                    pathlib.Path(cmd[cmd.index(flag) + 1]).write_bytes(b"x")
        return self.rc, "", ""


def run_main(tmp_path, runner, extra=(), monkeypatch=None):
    md = inputs(tmp_path)
    monkeypatch.setattr(f, "run_cmd", runner)
    return f.main(["-b", str(tmp_path / "table.qza"), "-r", str(tmp_path / "rep_seqs.qza"),
                   "-t", str(tmp_path / "taxonomy.qza"), "-m", md,
                   "-o", str(tmp_path / "out"), *extra])


ALL_FILTERS = ["--drop-where", "in_study_metadata='no'",
               "--min-samples-fraction", "0.5", "--min-sample-reads", "100"]


# ---- the full chain ------------------------------------------------------

def test_every_step_is_accounted_for(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(), extra=ALL_FILTERS,
                    monkeypatch=monkeypatch) == 0
    out = tmp_path / "out"
    # the file opens with comment lines recording the settings, so skip them first
    body = [ln for ln in (out / "filter_summary.tsv").read_text().splitlines()
            if not ln.startswith("#")]
    rows = {r["step"]: r for r in csv.DictReader(body, delimiter="\t")}
    assert rows["input"]["samples"] == "3" and rows["input"]["reads"] == "371"
    assert rows["samples"]["samples"] == "2" and rows["samples"]["samples_lost"] == "1"
    assert rows["prevalence"]["features"] == "2"
    assert rows["prevalence"]["features_lost"] == "1"
    assert (out / "table_filtered.qza").exists()
    assert (out / "rep_seqs_filtered.qza").exists()


def test_the_settings_land_next_to_the_result(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(), extra=ALL_FILTERS,
                    monkeypatch=monkeypatch) == 0
    text = (tmp_path / "out" / "filter_summary.tsv").read_text()
    assert "# min_samples_fraction: 0.5" in text
    assert "# drop_where: in_study_metadata='no'" in text
    assert "choices recorded with the result" in text


def test_per_sample_before_and_after(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(), extra=ALL_FILTERS,
                    monkeypatch=monkeypatch) == 0
    rows = {r["sample-id"]: r for r in
            csv.DictReader((tmp_path / "out" / "filter_per_sample.tsv").open(),
                           delimiter="\t")}
    assert rows["mock1"]["dropped"] == "yes" and rows["mock1"]["reads_after"] == "0"
    assert rows["s1"]["reads_before"] == "151" and rows["s1"]["reads_after"] == "150"


def test_taxonomy_filter_flags_are_explicit(tmp_path, monkeypatch):
    runner = FakeRunner()
    assert run_main(tmp_path, runner, monkeypatch=monkeypatch) == 0
    cmd = [c for c in runner.calls if "filter-table" in c][0]
    assert cmd[cmd.index("--p-include") + 1] == "Bacteria"
    assert cmd[cmd.index("--p-exclude") + 1] == "mitochondria,chloroplast"
    assert cmd[cmd.index("--p-mode") + 1] == "contains"


def test_sequences_are_filtered_against_the_final_table(tmp_path, monkeypatch):
    runner = FakeRunner()
    assert run_main(tmp_path, runner, extra=ALL_FILTERS, monkeypatch=monkeypatch) == 0
    cmd = [c for c in runner.calls if "filter-seqs" in c][0]
    assert cmd[cmd.index("--i-table") + 1].endswith("table_filtered.qza")


# ---- the guards ----------------------------------------------------------

def test_samples_lost_to_a_feature_filter_are_called_out(tmp_path, monkeypatch):
    """filter-features drops empty samples by default, which would otherwise change the
    sample count without a word."""
    tables = dict(TABLES)
    tables["03_prevalence"] = {"s1": {"a": 100, "b": 50}}      # s2 vanishes
    tables["04_depth"] = tables["03_prevalence"]
    assert run_main(tmp_path, FakeRunner(tables), extra=ALL_FILTERS,
                    monkeypatch=monkeypatch) == 0
    log = (tmp_path / "out" / "filter_log.txt").read_text()
    assert "removed by the *feature* filter" in log and "s2" in log
    assert "--p-filter-empty-samples on by default" in log


def test_an_emptied_table_is_fatal(tmp_path, monkeypatch):
    tables = dict(TABLES)
    tables["01_taxa"] = {}
    assert run_main(tmp_path, FakeRunner(tables), monkeypatch=monkeypatch) == 1
    log = (tmp_path / "out" / "filter_log.txt").read_text()
    assert "filtering emptied the table" in log
    assert "a domain label that does not appear removes everything" in log


def test_a_drop_where_that_matches_nothing_warns(tmp_path, monkeypatch):
    tables = dict(TABLES)
    tables["02_samples"] = tables["01_taxa"]                    # nothing removed
    assert run_main(tmp_path, FakeRunner(tables),
                    extra=["--drop-where", "in_study_metadata='nope'"],
                    monkeypatch=monkeypatch) == 0
    assert "matched no samples" in (tmp_path / "out" / "filter_log.txt").read_text()


def test_drop_where_without_metadata_fires(tmp_path, monkeypatch):
    monkeypatch.setattr(f, "run_cmd", FakeRunner())
    for name in ("table.qza", "rep_seqs.qza", "taxonomy.qza"):
        (tmp_path / name).write_bytes(b"x")
    assert f.main(["-b", str(tmp_path / "table.qza"), "-r", str(tmp_path / "rep_seqs.qza"),
                   "-t", str(tmp_path / "taxonomy.qza"), "-o", str(tmp_path / "out"),
                   "--drop-where", "x='y'"]) == 1


def test_prevalence_is_off_by_default(tmp_path, monkeypatch):
    runner = FakeRunner()
    assert run_main(tmp_path, runner, monkeypatch=monkeypatch) == 0
    assert not any("filter-features" in c for c in runner.calls)


def test_a_ragged_exported_table_fires(tmp_path, monkeypatch):
    class Ragged(FakeRunner):
        def __call__(self, cmd, timeout):
            if "convert" in " ".join(cmd):
                self.calls.append(cmd)
                pathlib.Path(cmd[cmd.index("-o") + 1]).write_text(
                    "#OTU ID\ts1\ts2\na\t1.0\n")
                return 0, "", ""
            return super().__call__(cmd, timeout)
    assert run_main(tmp_path, Ragged(), monkeypatch=monkeypatch) == 1
    assert "values for 2 columns" in (tmp_path / "out" / "filter_log.txt").read_text()


@pytest.mark.parametrize("step", ["taxa filter-table", "filter-seqs"])
def test_a_failing_qiime_step_is_reported(tmp_path, monkeypatch, step):
    assert run_main(tmp_path, FakeRunner(fail_on=step), monkeypatch=monkeypatch) == 1
    assert "Plugin error from" in (tmp_path / "out" / "filter_log.txt").read_text()


@pytest.mark.parametrize("bad", [["--min-samples-fraction", "1"],
                                 ["--min-samples-fraction", "-0.1"],
                                 ["--min-sample-reads", "-5"]])
def test_out_of_range_options_fire(tmp_path, monkeypatch, bad):
    assert run_main(tmp_path, FakeRunner(), extra=bad, monkeypatch=monkeypatch) == 1


def test_missing_input_fires(tmp_path, monkeypatch):
    monkeypatch.setattr(f, "run_cmd", FakeRunner())
    assert f.main(["-b", str(tmp_path / "nope.qza"), "-r", str(tmp_path / "nope.qza"),
                   "-t", str(tmp_path / "nope.qza"), "-o", str(tmp_path / "out")]) == 1


# ---- the no-hang rule ----------------------------------------------------

def test_a_stalled_command_is_killed():
    start = time.monotonic()
    with pytest.raises(f.FilterError, match="timed out after 1 s"):
        f.run_cmd(["sh", "-c", "sleep 30"], timeout=1)
    assert time.monotonic() - start < 10
