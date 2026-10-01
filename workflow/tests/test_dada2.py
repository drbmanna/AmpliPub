"""Offline tests for 03_dada2.py. Each guard is shown to fire.

QIIME 2 is replaced by a fake runner that writes the three artifacts denoise-paired
would write. The stats table has the columns read out of q2_dada2/_denoise.py in the
2025.7 build. The timeout test starts real processes.
"""

import csv
import importlib.util
import logging
import pathlib
import sys
import time
import zipfile

import pytest

HERE = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("dada2", HERE.parent / "03_dada2.py")
d = importlib.util.module_from_spec(spec)
sys.modules["dada2"] = d
spec.loader.exec_module(d)

HEADER = ["sample-id", "input", "filtered", "percentage of input passed filter",
          "denoised", "merged", "percentage of input merged", "non-chimeric",
          "percentage of input non-chimeric"]
TYPES = ["#q2:types", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
         "numeric", "numeric"]


def stats_tsv(samples, header=HEADER, types=True):
    """samples: (id, input, filtered, denoised, merged, non-chimeric) tuples."""
    lines = ["\t".join(header)]
    if types:
        lines.append("\t".join(types if isinstance(types, list) else TYPES))
    for sid, inp, filt, den, mer, nochim in samples:
        pf = 100 * filt / inp if inp else 0
        pm = 100 * mer / inp if inp else 0
        pn = 100 * nochim / inp if inp else 0
        lines.append("\t".join(str(x) for x in
                               [sid, inp, filt, round(pf, 2), den, mer, round(pm, 2),
                                nochim, round(pn, 2)]))
    return "\n".join(lines) + "\n"


HEALTHY = [("s1", 10000, 9000, 8800, 8500, 8200),
           ("s2", 20000, 18000, 17600, 17000, 16500)]


def make_stats_qza(path, text):
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("uuid/metadata.yaml", "type: SampleData[DADA2Stats]\n")
        zf.writestr("uuid/data/stats.tsv", text)
    return str(path)


def trunc_file(tmp_path, f=251, r=218, amplicon=253, overlap=12, margin=20, extra=None):
    p = tmp_path / "trunc_len.tsv"
    rows = [("trunc_len_f", f), ("trunc_len_r", r), ("min_q", 30),
            ("amplicon_len", amplicon), ("min_overlap", overlap), ("margin", margin),
            ("expected_overlap", f + r - amplicon), ("n_sampled", 10000)]
    if extra is not None:
        rows = extra
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerows(rows)
    return str(p)


def read_tsv(path):
    """Rows of a TSV, ignoring the comment lines the stage writes above the header."""
    lines = [l for l in pathlib.Path(path).read_text().splitlines()
             if l and not l.startswith("#")]
    return list(csv.DictReader(lines, delimiter="	"))


def read_tsv(path):
    """Rows of a TSV, ignoring the comment lines a stage writes above the header."""
    lines = [l for l in pathlib.Path(path).read_text().splitlines()
             if l and not l.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def fasta(lengths):
    """A FASTA of ASVs with the given lengths, named asv1, asv2, ... in order."""
    return "".join(f">asv{i}\n{'A' * L}\n" for i, L in enumerate(lengths, 1))


# A V4 catalogue that sits on the configured 253 bp, within the 20 bp of natural
# variation 02_quality.py's MARGIN allows for.
IN_RANGE = [253, 253, 252, 254, 248, 258]

# The shape found on PRJNA643648's V4 arm: a median that is exactly right, and a tail of
# identical-length ASVs at 1.7x the amplicon that no check in the pipeline looked at.
OVERSIZE_TAIL = IN_RANGE + [441, 441, 446, 456]


class FakeRunner:
    """Stands in for d.run_cmd. Writes the artifacts denoise-paired would write.

    It also answers `qiime tools export`, which the ASV length check uses to get at the
    sequences, by writing the FASTA the export would have produced.
    """

    def __init__(self, stats=None, rc=0, write=True, skip=(), asvs=None):
        self.stats = stats if stats is not None else stats_tsv(HEALTHY)
        self.asvs = IN_RANGE if asvs is None else asvs
        self.rc, self.write, self.skip = rc, write, skip
        self.calls = []

    def __call__(self, cmd, timeout):
        self.calls.append(cmd)
        assert timeout > 0
        if "export" in cmd:
            if self.write and "export" not in self.skip:
                dest = pathlib.Path(cmd[cmd.index("--output-path") + 1])
                dest.mkdir(parents=True, exist_ok=True)
                (dest / "dna-sequences.fasta").write_text(fasta(self.asvs))
            return self.rc, "", "Plugin error from tools\n" if self.rc else ""
        if self.write:
            for flag in ("--o-table", "--o-representative-sequences"):
                if flag not in self.skip:
                    pathlib.Path(cmd[cmd.index(flag) + 1]).write_bytes(b"artifact")
            if "--o-denoising-stats" not in self.skip:
                make_stats_qza(cmd[cmd.index("--o-denoising-stats") + 1], self.stats)
        return self.rc, "", "Plugin error from dada2\n" if self.rc else ""


def run_main(monkeypatch, tmp_path, runner, extra=(), trunc=None):
    qza = tmp_path / "trimmed.qza"
    qza.write_bytes(b"x")
    monkeypatch.setattr(d, "run_cmd", runner)
    return d.main(["-i", str(qza), "-t", trunc or trunc_file(tmp_path),
                   "-o", str(tmp_path / "out"), *extra])


# ---- trunc_len.tsv -------------------------------------------------------

def test_trunc_len_is_read(tmp_path):
    v = d.read_trunc_len(trunc_file(tmp_path))
    assert v["trunc_len_f"] == 251 and v["trunc_len_r"] == 218 and v["amplicon_len"] == 253


def test_missing_key_fires(tmp_path):
    path = trunc_file(tmp_path, extra=[("trunc_len_f", 251), ("trunc_len_r", 218)])
    with pytest.raises(d.Dada2Error, match="amplicon_len"):
        d.read_trunc_len(path)


def test_non_numeric_value_fires(tmp_path):
    path = trunc_file(tmp_path, extra=[("trunc_len_f", "long")])
    with pytest.raises(d.Dada2Error, match="not a number"):
        d.read_trunc_len(path)


def test_wrong_column_count_fires(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("trunc_len_f\t251\textra\n")
    with pytest.raises(d.Dada2Error, match="expected two columns"):
        d.read_trunc_len(str(p))


def test_zero_trunc_len_fires(tmp_path):
    path = trunc_file(tmp_path, f=0)
    with pytest.raises(d.Dada2Error, match="trunc_len_f is 0"):
        d.read_trunc_len(path)


# ---- the overlap guard ---------------------------------------------------

def test_overlap_guard_fires_below_the_floor():
    v = {"trunc_len_f": 150, "trunc_len_r": 130, "amplicon_len": 253,
         "min_overlap": 12, "margin": 20}
    with pytest.raises(d.Dada2Error, match="would not merge"):
        d.check_overlap(v)


def test_overlap_guard_passes_at_the_floor():
    v = {"trunc_len_f": 200, "trunc_len_r": 85, "amplicon_len": 253,
         "min_overlap": 12, "margin": 20}
    assert d.check_overlap(v) == 32  # 285 - 253


def test_dev_lengths_pass_the_guard(tmp_path):
    assert d.check_overlap(d.read_trunc_len(trunc_file(tmp_path))) == 216


# ---- the stats table -----------------------------------------------------

def test_types_row_is_not_a_sample():
    rows = d.parse_stats(stats_tsv(HEALTHY), "x")
    assert [r["sample-id"] for r in rows] == ["s1", "s2"]
    assert rows[0]["non-chimeric"] == 8200


def test_single_end_table_is_refused():
    header = [c for c in HEADER if c not in ("merged", "percentage of input merged")]
    text = "\n".join(["\t".join(header), "s1\t10\t9\t90.0\t9\t8\t80.0"]) + "\n"
    with pytest.raises(d.Dada2Error, match="no merged column"):
        d.parse_stats(text, "x")


def test_unexpected_first_column_fires():
    header = ["id"] + HEADER[1:]
    with pytest.raises(d.Dada2Error, match="expected 'sample-id'"):
        d.parse_stats(stats_tsv(HEALTHY, header=header), "x")


def test_non_numeric_count_fires():
    text = stats_tsv(HEALTHY).replace("\t8200\t", "\tmany\t")
    with pytest.raises(d.Dada2Error, match="non-numeric"):
        d.parse_stats(text, "x")


def test_table_with_only_a_types_row_fires():
    with pytest.raises(d.Dada2Error, match="holds no samples"):
        d.parse_stats(stats_tsv([]), "x")


def test_bad_artifact_fires(tmp_path):
    p = tmp_path / "stats.qza"
    p.write_bytes(b"not a zip")
    with pytest.raises(d.Dada2Error, match="not a readable artifact"):
        d.read_stats(str(p))


def test_artifact_without_stats_fires(tmp_path):
    p = tmp_path / "stats.qza"
    with zipfile.ZipFile(p, "w") as zf:
        zf.writestr("uuid/metadata.yaml", "type: x\n")
    with pytest.raises(d.Dada2Error, match="expected one stats table"):
        d.read_stats(str(p))


# ---- the flag rules ------------------------------------------------------

def test_no_flags_on_a_healthy_run():
    flags, totals = d.judge(d.parse_stats(stats_tsv(HEALTHY), "x"), 0.5)
    assert flags == []
    assert totals["input"] == 30000 and totals["non-chimeric"] == 24700


def test_majority_lost_at_one_step_is_flagged():
    rows = d.parse_stats(stats_tsv([("s1", 10000, 9000, 8800, 4000, 3900)]), "x")
    flags, _ = d.judge(rows, 0.5)
    assert [f["step"] for f in flags] == ["denoised->merged"]
    assert flags[0]["lost_fraction"] == pytest.approx(0.5455, abs=1e-4)


def test_chimera_loss_is_flagged():
    rows = d.parse_stats(stats_tsv([("s1", 10000, 9000, 8800, 8500, 2000)]), "x")
    flags, _ = d.judge(rows, 0.5)
    assert [f["step"] for f in flags] == ["merged->non-chimeric"]


def test_filtering_loss_is_never_flagged():
    """The tutorial rule is explicitly "outside of filtering"."""
    rows = d.parse_stats(stats_tsv([("s1", 10000, 1000, 980, 950, 940)]), "x")
    flags, _ = d.judge(rows, 0.5)
    assert flags == []


def test_lost_fraction_of_nothing_is_zero():
    assert d.lost_fraction(0, 0) == 0.0


# ---- end to end ----------------------------------------------------------

def test_healthy_run_writes_everything(tmp_path, monkeypatch):
    runner = FakeRunner()
    assert run_main(monkeypatch, tmp_path, runner) == 0
    out = tmp_path / "out"
    rows = list(csv.DictReader((out / "dada2_stats.tsv").open(), delimiter="\t"))
    assert [r["sample-id"] for r in rows] == ["s1", "s2"]
    assert rows[0]["retained_of_input"] == "0.82"
    assert (out / "dada2_flags.tsv").read_text().strip().endswith("note")  # header only
    criteria = (out / "criteria.tsv").read_text()
    assert "majority of reads are lost" in criteria and "hard fail" in criteria
    assert (out / "table.qza").exists() and (out / "rep_seqs.qza").exists()


def test_parameters_are_passed_explicitly(tmp_path, monkeypatch):
    runner = FakeRunner()
    assert run_main(monkeypatch, tmp_path, runner, extra=["--threads", "8"]) == 0
    cmd = runner.calls[0]
    assert cmd[cmd.index("--p-n-threads") + 1] == "8"
    assert cmd[cmd.index("--p-trunc-len-f") + 1] == "251"
    assert cmd[cmd.index("--p-trunc-len-r") + 1] == "218"
    assert cmd[cmd.index("--p-min-overlap") + 1] == "12"


def test_zero_read_library_is_a_hard_fail(tmp_path, monkeypatch):
    runner = FakeRunner(stats=stats_tsv(HEALTHY + [("s3", 900, 100, 90, 0, 0)]))
    assert run_main(monkeypatch, tmp_path, runner) == 1
    assert "s3" in (tmp_path / "out" / "dada2_log.txt").read_text()


def test_zero_read_library_can_be_overridden(tmp_path, monkeypatch):
    runner = FakeRunner(stats=stats_tsv(HEALTHY + [("s3", 900, 100, 90, 0, 0)]))
    assert run_main(monkeypatch, tmp_path, runner,
                    extra=["--allow-zero-read-samples"]) == 0
    log = (tmp_path / "out" / "dada2_log.txt").read_text()
    assert "hard check was overridden" in log
    assert "override passed" in (tmp_path / "out" / "criteria.tsv").read_text()


def test_empty_table_is_a_hard_fail(tmp_path, monkeypatch):
    runner = FakeRunner(stats=stats_tsv([("s1", 10000, 9000, 8800, 0, 0)]))
    assert run_main(monkeypatch, tmp_path, runner,
                    extra=["--allow-zero-read-samples"]) == 1
    assert "feature table is empty" in (tmp_path / "out" / "dada2_log.txt").read_text()


def test_flags_are_written_and_logged(tmp_path, monkeypatch):
    runner = FakeRunner(stats=stats_tsv([("s1", 10000, 9000, 8800, 8500, 2000),
                                         ("s2", 20000, 18000, 17600, 17000, 16500)]))
    assert run_main(monkeypatch, tmp_path, runner) == 0
    flags = list(csv.DictReader((tmp_path / "out" / "dada2_flags.tsv").open(), delimiter="\t"))
    assert len(flags) == 1 and flags[0]["sample-id"] == "s1"
    log = (tmp_path / "out" / "dada2_log.txt").read_text()
    assert "removed as chimeric" in log  # the quote travels with the flag


def test_qiime_failure_is_reported(tmp_path, monkeypatch):
    assert run_main(monkeypatch, tmp_path, FakeRunner(rc=1, write=False)) == 1
    assert "Plugin error from dada2" in (tmp_path / "out" / "dada2_log.txt").read_text()


def test_missing_output_artifact_fires(tmp_path, monkeypatch):
    runner = FakeRunner(skip=("--o-table",))
    assert run_main(monkeypatch, tmp_path, runner) == 1
    assert "wrote no" in (tmp_path / "out" / "dada2_log.txt").read_text()


def test_overlap_guard_stops_before_any_compute(tmp_path, monkeypatch):
    runner = FakeRunner()
    trunc = trunc_file(tmp_path, f=150, r=130)
    assert run_main(monkeypatch, tmp_path, runner, trunc=trunc) == 1
    assert runner.calls == []  # nothing was run


def test_missing_input_fires(tmp_path, monkeypatch):
    monkeypatch.setattr(d, "run_cmd", FakeRunner())
    assert d.main(["-i", str(tmp_path / "nope.qza"), "-t", trunc_file(tmp_path),
                   "-o", str(tmp_path / "out")]) == 1


def test_empty_input_fires(tmp_path, monkeypatch):
    qza = tmp_path / "trimmed.qza"
    qza.write_bytes(b"")
    monkeypatch.setattr(d, "run_cmd", FakeRunner())
    assert d.main(["-i", str(qza), "-t", trunc_file(tmp_path),
                   "-o", str(tmp_path / "out")]) == 1


@pytest.mark.parametrize("bad", [["--majority", "0"], ["--majority", "1"],
                                 ["--threads", "0"]])
def test_out_of_range_options_fire(tmp_path, monkeypatch, bad):
    assert run_main(monkeypatch, tmp_path, FakeRunner(), extra=bad) == 1


# ---- the ASV length check ------------------------------------------------
#
# Written after 8.9% of a V4 catalogue was found sitting at 294-456 bp on a 253 bp
# amplicon, having passed every check in the pipeline. The median was 253 bp exactly, so
# a median test alone does not fire: each test below says which of the two it exercises.

def test_asv_lengths_are_measured_and_written(tmp_path, monkeypatch):
    assert run_main(monkeypatch, tmp_path, FakeRunner()) == 0
    out = tmp_path / "out"
    rows = list(csv.DictReader((out / "asv_lengths.tsv").open(), delimiter="\t"))
    assert len(rows) == len(IN_RANGE)
    assert {r["class"] for r in rows} == {"in_range"}
    assert rows[0]["length"] == "258" and rows[0]["difference_from_amplicon_len"] == "5"
    # Nothing was off-length, so the flags file is a header and nothing else.
    body = [l for l in (out / "asv_length_flags.tsv").read_text().splitlines()
            if not l.startswith("#")]
    assert body == ["feature-id\tlength\tdifference_from_amplicon_len\tclass"]


def test_the_oversize_guard_fires_on_a_tail_the_median_cannot_see(tmp_path, monkeypatch,
                                                                 caplog):
    """The guard the finding asked for. Median 253 bp, four ASVs at up to 1.8x."""
    runner = FakeRunner(asvs=OVERSIZE_TAIL)
    with caplog.at_level(logging.WARNING, logger="dada2"):
        rc = run_main(monkeypatch, tmp_path, runner)
    assert rc == 0  # flagged and reported, not failed: the table is usable
    out = tmp_path / "out"
    flagged = read_tsv(out / "asv_length_flags.tsv")
    assert [r["length"] for r in flagged] == ["456", "446", "441", "441"]
    assert {r["class"] for r in flagged} == {"oversize"}
    assert "4 of 10 ASVs (40.0%) are oversize, above 303 bp: 441-456 bp" in caplog.text
    assert "2 at 441 bp" in caplog.text          # the cluster, not scattered noise
    # the threshold is stated as a choice, not implied to be a standard
    assert "no published source sets a length threshold" in caplog.text
    assert "flagged, not removed" in caplog.text
    lengths = {r["feature-id"]: r["class"] for r in
               csv.DictReader((out / "asv_lengths.tsv").open(), delimiter="\t")}
    assert sum(v == "oversize" for v in lengths.values()) == 4
    assert sum(v == "in_range" for v in lengths.values()) == len(IN_RANGE)


def test_the_oversize_guard_stays_quiet_inside_the_tolerance(tmp_path, monkeypatch):
    """At the cutoff exactly, not flagged: the guard must not fire on its own boundary.

    303 bp is 253 + 50. The 268-293 bp primer readthrough band 01_primers.py already
    handles sits below it on purpose and is not reported twice.
    """
    runner = FakeRunner(asvs=IN_RANGE + [268, 293, 303])
    assert run_main(monkeypatch, tmp_path, runner) == 0
    rows = list(csv.DictReader((tmp_path / "out" / "asv_lengths.tsv").open(),
                               delimiter="\t"))
    assert {r["class"] for r in rows} == {"in_range"}


def test_the_undersize_guard_fires_as_well(tmp_path, monkeypatch, caplog):
    runner = FakeRunner(asvs=IN_RANGE + [120, 202])
    with caplog.at_level(logging.WARNING, logger="dada2"):
        assert run_main(monkeypatch, tmp_path, runner) == 0
    assert "2 of 8 ASVs (25.0%) are undersize, below 203 bp: 120-202 bp" in caplog.text
    flagged = read_tsv(tmp_path / "out" / "asv_length_flags.tsv")
    assert [r["length"] for r in flagged] == ["120", "202"]
    assert {r["class"] for r in flagged} == {"undersize"}


def test_the_median_guard_hard_fails_when_the_region_is_wrong(tmp_path, monkeypatch):
    """V3-V4 ASVs denoised with amplicon_len left at V4's 253: the whole run is wrong."""
    runner = FakeRunner(asvs=[427, 427, 428, 426])
    assert run_main(monkeypatch, tmp_path, runner) == 1
    text = (tmp_path / "out" / "dada2_log.txt").read_text()
    assert "do not describe the configured region" in text
    assert "the denoised ASVs have a median of 427 bp" in text
    assert "what was actually sequenced" in text


def test_the_median_guard_fires_one_bp_past_the_tolerance():
    assert d.check_asv_lengths({"a": 303}, 253, 50)["median"] == 303
    with pytest.raises(d.Dada2Error, match="do not describe the configured region"):
        d.check_asv_lengths({"a": 304}, 253, 50)
    with pytest.raises(d.Dada2Error, match="do not describe the configured region"):
        d.check_asv_lengths({"a": 202}, 253, 50)


def test_an_empty_catalogue_fires(tmp_path, monkeypatch):
    assert run_main(monkeypatch, tmp_path, FakeRunner(asvs=[])) == 1
    assert "cannot read the exported ASVs" in (tmp_path / "out" / "dada2_log.txt").read_text()


def test_a_failed_export_fires(tmp_path, monkeypatch):
    assert run_main(monkeypatch, tmp_path, FakeRunner(skip=("export",))) == 1
    assert "wrote no" in (tmp_path / "out" / "dada2_log.txt").read_text()


def test_the_tolerance_is_a_flag_and_is_recorded(tmp_path, monkeypatch):
    runner = FakeRunner(asvs=OVERSIZE_TAIL)
    assert run_main(monkeypatch, tmp_path, runner,
                    extra=["--asv-len-tolerance", "250"]) == 0
    rows = list(csv.DictReader((tmp_path / "out" / "asv_lengths.tsv").open(),
                               delimiter="\t"))
    assert {r["class"] for r in rows} == {"in_range"}  # 3-503 bp, so nothing is off
    criteria = (tmp_path / "out" / "criteria.tsv").read_text()
    assert "asv_len_tolerance\tparameter\t250" in criteria
    assert "none of the six sources sets a length threshold" in criteria


def test_a_negative_tolerance_fires(tmp_path, monkeypatch):
    assert run_main(monkeypatch, tmp_path, FakeRunner(),
                    extra=["--asv-len-tolerance", "-1"]) == 1


# ---- the no-hang rule ----------------------------------------------------

def test_a_stalled_command_is_killed():
    start = time.monotonic()
    with pytest.raises(d.Dada2Error, match="timed out after 1 s"):
        d.run_cmd(["sh", "-c", "sleep 30"], timeout=1)
    assert time.monotonic() - start < 10
