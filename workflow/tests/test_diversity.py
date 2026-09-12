"""Offline tests for 10_diversity.py. Each guard is shown to fire.

The interesting ones: the stage must refuse to invent a depth, it must find the natural
break Schloss 2024 describes, and it must say out loud that core-metrics subsamples once.
"""

import csv
import importlib.util
import pathlib
import sys
import time

import pytest

HERE = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("diversity", HERE.parent / "10_diversity.py")
dv = importlib.util.module_from_spec(spec)
sys.modules["diversity"] = dv
spec.loader.exec_module(dv)

# s4 and s5 are deep; s1..s3 sit below a clear gap, as in the real study
TABLE = {"s1": [900.0, 100.0], "s2": [950.0, 50.0], "s3": [980.0, 20.0],
         "s4": [3000.0, 2000.0], "s5": [4000.0, 2000.0], "s6": [5000.0, 1000.0]}
# the pre-filter table keeps singletons, which Good's coverage needs
COVERAGE_TABLE = {name: values + [1.0, 1.0] for name, values in TABLE.items()}

def wide(n_shallow, n_deep):
    """A realistic depth distribution: a short low tail below a clear gap."""
    out = {}
    for i in range(n_shallow):
        out[f"low{i:02d}"] = [900.0 + i * 10, 100.0]
    for i in range(n_deep):
        out[f"deep{i:02d}"] = [4000.0 + i * 100, 1000.0]
    return out


# 1 of 13 below the break: found, and cheap enough to use
WIDE_CHEAP = wide(1, 12)
# 4 of 20 below the break: found, but 20% loss is past the ceiling
WIDE_COSTLY = wide(4, 16)

CORE_FILES = ["faith_pd_vector.qza", "shannon_vector.qza",
              "observed_features_vector.qza", "unweighted_unifrac_distance_matrix.qza",
              "weighted_unifrac_distance_matrix.qza",
              "bray_curtis_distance_matrix.qza"]


def biom_tsv(counts):
    cols = list(counts)
    n = max(len(v) for v in counts.values())
    lines = ["# Constructed from biom file", "#OTU ID\t" + "\t".join(cols)]
    for i in range(n):
        row = [f"f{i}"]
        for col in cols:
            values = counts[col]
            row.append(str(values[i] if i < len(values) else 0.0))
        lines.append("\t".join(row))
    return "\n".join(lines) + "\n"


# ---- the depth distribution ----------------------------------------------

def test_depths_are_summed_per_sample():
    assert dv.depths(TABLE)["s1"] == 1000 and dv.depths(TABLE)["s5"] == 6000


def test_the_natural_break_is_found():
    """Schloss 2024's stated method, made countable."""
    values = [1000, 1000, 1000, 5000, 5100, 5200, 5300, 6000]
    got = dv.natural_break(values, window=0.5)
    assert got["below"] == 3 and got["low"] == 1000 and got["high"] == 5000
    assert got["ratio"] == pytest.approx(5.0)


def test_a_smooth_distribution_has_no_break():
    values = [1000 + 100 * i for i in range(40)]
    got = dv.natural_break(values)
    # every gap is tiny and roughly equal; the largest is still barely above 1
    assert got is None or got["ratio"] < 1.2


def test_a_break_needs_enough_samples():
    assert dv.natural_break([10, 20, 30]) is None


def test_only_the_low_tail_is_searched():
    """A gap between two deep samples says nothing about where to cut."""
    values = [1000, 1100, 1200, 1300, 1400, 1500, 1600, 90000]
    got = dv.natural_break(values, window=0.25)
    assert got is None or got["high"] != 90000


# ---- Good's coverage ------------------------------------------------------

def test_goods_coverage_counts_singletons():
    got = dv.goods_coverage({"s1": [98.0, 1.0, 1.0]})
    assert got["s1"] == pytest.approx(1 - 2 / 100)


def test_coverage_of_an_empty_sample_is_zero():
    assert dv.goods_coverage({"s1": [0.0, 0.0]})["s1"] == 0.0


def test_a_prevalence_filtered_table_gives_a_perfect_score():
    """Why the stage insists on the pre-filter table: no singletons, coverage 1.0."""
    assert dv.goods_coverage({"s1": [500.0, 500.0]})["s1"] == 1.0


# ---- end to end -----------------------------------------------------------

class FakeRunner:
    def __init__(self, rc=0, write=True, fail_on=None, core_files=None, table=None):
        self.rc, self.write, self.fail_on = rc, write, fail_on
        self.table = TABLE if table is None else table
        self.core_files = CORE_FILES if core_files is None else core_files
        self.calls = []
        self._tag = "table"

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
        elif "convert" in joined:
            counts = (COVERAGE_TABLE if self._tag == "coverage"
                      else self.table)
            pathlib.Path(cmd[cmd.index("-o") + 1]).write_text(biom_tsv(counts))
        elif "core-metrics-phylogenetic" in cmd:
            core = pathlib.Path(cmd[cmd.index("--output-dir") + 1])
            core.mkdir(parents=True, exist_ok=True)
            for name in self.core_files:
                (core / name).write_bytes(b"x")
        elif "--o-visualization" in cmd:
            pathlib.Path(cmd[cmd.index("--o-visualization") + 1]).write_bytes(b"x")
        return self.rc, "", ""


def run_main(tmp_path, runner, extra=(), monkeypatch=None):
    for name in ("table.qza", "tree.qza"):
        (tmp_path / name).write_bytes(b"x")
    md = tmp_path / "md.tsv"
    md.write_text("sample-id\tdx\ns1\tcase\ns2\tcase\ns3\tcontrol\ns4\tcontrol\n"
                  "s5\tcase\ns6\tcontrol\n")
    monkeypatch.setattr(dv, "run_cmd", runner)
    return dv.main(["-b", str(tmp_path / "table.qza"), "-p", str(tmp_path / "tree.qza"),
                    "-m", str(md), "-o", str(tmp_path / "out"), *extra])


def test_it_refuses_to_invent_a_depth(tmp_path, monkeypatch, caplog):
    """The argument check fires before any output directory is made, so the message
    goes to the logger's fallback rather than to a file."""
    with caplog.at_level("ERROR", logger="diversity"):
        assert run_main(tmp_path, FakeRunner(), monkeypatch=monkeypatch) == 1
    assert "will not invent one" in caplog.text


def test_auto_depth_refuses_a_break_that_costs_too_many_samples(tmp_path, monkeypatch):
    """The break and the sample-loss ceiling are independent checks, and the ceiling
    still applies to a depth the stage chose itself."""
    assert run_main(tmp_path, FakeRunner(table=WIDE_COSTLY), extra=["--auto-depth"],
                    monkeypatch=monkeypatch) == 1
    log = (tmp_path / "out" / "diversity_log.txt").read_text()
    assert "natural break in the low tail: 1030 -> 5000" in log
    assert "past the 10% you allowed" in log


def test_no_break_with_auto_depth_is_refused_not_guessed(tmp_path, monkeypatch):
    """A distribution with no low tail gives no break, and the stage says the choice is
    arbitrary instead of picking something."""
    assert run_main(tmp_path, FakeRunner(), extra=["--auto-depth"],
                    monkeypatch=monkeypatch) == 1
    log = (tmp_path / "out" / "diversity_log.txt").read_text()
    assert "no natural break stands out" in log
    assert "the data does not point anywhere" in log


def test_auto_depth_uses_the_break_and_says_so(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(table=WIDE_CHEAP), extra=["--auto-depth"],
                    monkeypatch=monkeypatch) == 0
    log = (tmp_path / "out" / "diversity_log.txt").read_text()
    assert "natural break in the low tail: 1000 -> 5000" in log
    assert "depth 5000, the natural break" in log
    assert "keeps 12 of 13 samples" in log
    rows = {r["setting"]: r for r in
            csv.DictReader((tmp_path / "out" / "diversity_settings.tsv").open(),
                           delimiter="\t")}
    assert rows["sampling_depth"]["value"] == "5000"
    assert "natural break" in rows["sampling_depth"]["source"]


def test_an_explicit_depth_is_honoured(tmp_path, monkeypatch):
    runner = FakeRunner()
    assert run_main(tmp_path, runner, extra=["--depth", "1000"],
                    monkeypatch=monkeypatch) == 0
    cmd = [c for c in runner.calls if "core-metrics-phylogenetic" in c][0]
    assert cmd[cmd.index("--p-sampling-depth") + 1] == "1000"


def test_the_single_subsample_caveat_is_always_logged(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(), extra=["--depth", "1000"],
                    monkeypatch=monkeypatch) == 0
    log = (tmp_path / "out" / "diversity_log.txt").read_text()
    assert "subsamples ONCE" in log and "Schloss 2024" in log
    rows = {r["setting"]: r["value"] for r in
            csv.DictReader((tmp_path / "out" / "diversity_settings.tsv").open(),
                           delimiter="\t")}
    assert rows["core_metrics_subsampling"] == "single"


def test_curve_iterations_are_passed(tmp_path, monkeypatch):
    runner = FakeRunner()
    assert run_main(tmp_path, runner, extra=["--depth", "1000", "--iterations", "250"],
                    monkeypatch=monkeypatch) == 0
    cmd = [c for c in runner.calls if "alpha-rarefaction" in c][0]
    assert cmd[cmd.index("--p-iterations") + 1] == "250"


def test_too_much_sample_loss_is_fatal(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(), extra=["--depth", "5500"],
                    monkeypatch=monkeypatch) == 1
    log = (tmp_path / "out" / "diversity_log.txt").read_text()
    assert "past the 10% you allowed" in log
    assert "record why" in log


def test_the_loss_ceiling_can_be_raised_deliberately(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(),
                    extra=["--depth", "5500", "--max-sample-loss", "0.9"],
                    monkeypatch=monkeypatch) == 0


def test_a_depth_above_every_sample_is_fatal(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(),
                    extra=["--depth", "99999", "--max-sample-loss", "0.99"],
                    monkeypatch=monkeypatch) == 1
    assert "above every sample's read count" in \
        (tmp_path / "out" / "diversity_log.txt").read_text()


def test_the_candidate_table_is_written(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(), extra=["--depth", "1000"],
                    monkeypatch=monkeypatch) == 0
    rows = {int(r["depth"]): r for r in
            csv.DictReader((tmp_path / "out" / "depth_candidates.tsv").open(),
                           delimiter="\t")}
    assert rows[1000]["samples_kept"] == "6" and rows[1000]["samples_lost"] == "0"
    assert rows[5000]["samples_kept"] == "3"


# ---- Good's coverage plumbing --------------------------------------------

def test_coverage_is_reported_from_the_pre_filter_table(tmp_path, monkeypatch):
    (tmp_path / "pre.qza").write_bytes(b"x")
    assert run_main(tmp_path, FakeRunner(),
                    extra=["--depth", "1000", "--coverage-table",
                           str(tmp_path / "pre.qza")], monkeypatch=monkeypatch) == 0
    log = (tmp_path / "out" / "diversity_log.txt").read_text()
    assert "Good's coverage over 6 shared sample(s)" in log
    assert "only place it can mean anything" in log
    rows = {r["sample-id"]: r for r in
            csv.DictReader((tmp_path / "out" / "sample_depths.tsv").open(),
                           delimiter="\t")}
    assert float(rows["s1"]["goods_coverage"]) == pytest.approx(1 - 2 / 1002)


def test_a_table_with_no_singletons_is_called_degenerate(tmp_path, monkeypatch):
    """DADA2 never emits singleton ASVs, so Good's coverage is 1.0 by construction and
    must not be presented as evidence of depth."""
    no_singletons = {name: [500.0, 500.0] for name in TABLE}
    class NoSingletons(FakeRunner):
        def __call__(self, cmd, timeout):
            if "convert" in " ".join(cmd) and self._tag == "coverage":
                self.calls.append(cmd)
                pathlib.Path(cmd[cmd.index("-o") + 1]).write_text(biom_tsv(no_singletons))
                return 0, "", ""
            return super().__call__(cmd, timeout)
    (tmp_path / "pre.qza").write_bytes(b"x")
    assert run_main(tmp_path, NoSingletons(),
                    extra=["--depth", "1000", "--coverage-table",
                           str(tmp_path / "pre.qza")], monkeypatch=monkeypatch) == 0
    log = (tmp_path / "out" / "diversity_log.txt").read_text()
    assert "not one singleton in this table" in log
    assert "NOT evidence that anything was sequenced deeply enough" in log
    assert "issue 1491" in log
    assert "Do not report it here" in log
    rows = {r["setting"]: r["value"] for r in
            csv.DictReader((tmp_path / "out" / "diversity_settings.tsv").open(),
                           delimiter="	")}
    assert rows["goods_coverage_status"] == "degenerate: no singletons in the table"


def test_singletons_present_reports_the_count(tmp_path, monkeypatch):
    (tmp_path / "pre.qza").write_bytes(b"x")
    assert run_main(tmp_path, FakeRunner(),
                    extra=["--depth", "1000", "--coverage-table",
                           str(tmp_path / "pre.qza")], monkeypatch=monkeypatch) == 0
    rows = {r["setting"]: r["value"] for r in
            csv.DictReader((tmp_path / "out" / "diversity_settings.tsv").open(),
                           delimiter="	")}
    assert rows["goods_coverage_status"] == "12 singleton(s) in the table"
    assert "not one singleton" not in (tmp_path / "out" / "diversity_log.txt").read_text()


def test_missing_coverage_table_is_explained_not_faked(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(), extra=["--depth", "1000"],
                    monkeypatch=monkeypatch) == 0
    log = (tmp_path / "out" / "diversity_log.txt").read_text()
    assert "no --coverage-table" in log
    assert "printed as if it were real" in log
    rows = {r["sample-id"]: r for r in
            csv.DictReader((tmp_path / "out" / "sample_depths.tsv").open(),
                           delimiter="\t")}
    assert rows["s1"]["goods_coverage"] == ""


# ---- group tests ---------------------------------------------------------

def test_group_tests_run_for_each_column(tmp_path, monkeypatch):
    runner = FakeRunner()
    assert run_main(tmp_path, runner, extra=["--depth", "1000", "--group-column", "dx"],
                    monkeypatch=monkeypatch) == 0
    beta = [c for c in runner.calls if "beta-group-significance" in c]
    assert len(beta) == 3
    assert beta[0][beta[0].index("--p-method") + 1] == "permanova"
    assert beta[0][beta[0].index("--p-permutations") + 1] == "999"
    assert len([c for c in runner.calls if "alpha-group-significance" in c]) == 3
    assert "3 group test(s)" not in (tmp_path / "out" / "diversity_log.txt").read_text()


def test_no_group_column_says_nothing_was_compared(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(), extra=["--depth", "1000"],
                    monkeypatch=monkeypatch) == 0
    assert "nothing was compared" in \
        (tmp_path / "out" / "diversity_log.txt").read_text()


def test_a_missing_core_metric_is_skipped_with_a_warning(tmp_path, monkeypatch):
    runner = FakeRunner(core_files=["shannon_vector.qza",
                                    "bray_curtis_distance_matrix.qza"])
    assert run_main(tmp_path, runner, extra=["--depth", "1000", "--group-column", "dx"],
                    monkeypatch=monkeypatch) == 0
    log = (tmp_path / "out" / "diversity_log.txt").read_text()
    assert "faith_pd_vector not found" in log


@pytest.mark.parametrize("step", ["alpha-rarefaction", "core-metrics-phylogenetic",
                                  "beta-group-significance"])
def test_a_failing_step_is_reported(tmp_path, monkeypatch, step):
    assert run_main(tmp_path, FakeRunner(fail_on=step),
                    extra=["--depth", "1000", "--group-column", "dx"],
                    monkeypatch=monkeypatch) == 1
    assert "Plugin error from" in (tmp_path / "out" / "diversity_log.txt").read_text()


@pytest.mark.parametrize("bad", [["--depth", "0"], ["--iterations", "0"],
                                 ["--max-sample-loss", "1"]])
def test_out_of_range_options_fire(tmp_path, monkeypatch, bad):
    assert run_main(tmp_path, FakeRunner(), extra=bad, monkeypatch=monkeypatch) == 1


def test_missing_input_fires(tmp_path, monkeypatch):
    monkeypatch.setattr(dv, "run_cmd", FakeRunner())
    assert dv.main(["-b", str(tmp_path / "nope.qza"), "-p", str(tmp_path / "nope.qza"),
                    "-m", str(tmp_path / "nope.tsv"), "-o", str(tmp_path / "out"),
                    "--depth", "1000"]) == 1


# ---- the no-hang rule ----------------------------------------------------

def test_a_stalled_command_is_killed():
    start = time.monotonic()
    with pytest.raises(dv.DiversityError, match="timed out after 1 s"):
        dv.run_cmd(["sh", "-c", "sleep 30"], timeout=1)
    assert time.monotonic() - start < 10
