"""Offline tests for 05_taxonomy.py. Each guard is shown to fire.

QIIME 2 is replaced by a fake runner that writes the artifacts and exports
classify-sklearn would produce. No classifier and no reference data is needed.
"""

import csv
import importlib.util
import pathlib
import sys
import time
import zipfile

import pytest

HERE = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("taxonomy", HERE.parent / "05_taxonomy.py")
t = importlib.util.module_from_spec(spec)
sys.modules["taxonomy"] = t
spec.loader.exec_module(t)

GG = "d__Bacteria; p__Bacillota; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__Streptococcus mutans"
SILVA = "d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__uncultured"


def artifact(path, kind, uuid="0000-1111"):
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr(f"{uuid}/metadata.yaml", f"uuid: {uuid}\ntype: {kind}\nformat: X\n")
    return str(path)


def taxonomy_tsv(rows, header=("Feature ID", "Taxon", "Confidence"), types=True):
    lines = ["\t".join(header)]
    if types:
        lines.append("#q2:types\tcategorical\tnumeric")
    for row in rows:
        lines.append("\t".join(str(c) for c in row))
    return "\n".join(lines) + "\n"


# ---- taxon strings -------------------------------------------------------

def test_prefixed_taxon_splits_by_rank():
    r = t.split_taxon(GG)
    assert r["domain"] == "Bacteria" and r["genus"] == "Streptococcus"
    assert r["species"] == "Streptococcus mutans"


def test_k_prefix_is_treated_as_domain():
    r = t.split_taxon("k__Bacteria; p__Firmicutes")
    assert r["domain"] == "Bacteria" and r["phylum"] == "Firmicutes"


def test_empty_and_placeholder_labels_do_not_count():
    r = t.split_taxon("d__Bacteria; p__; c__uncultured; o__Unassigned")
    assert r == {"domain": "Bacteria"}


def test_unprefixed_taxon_falls_back_to_position():
    r = t.split_taxon("Bacteria; Firmicutes; Bacilli")
    assert r == {"domain": "Bacteria", "phylum": "Firmicutes", "class": "Bacilli"}


def test_a_truncated_taxon_stops_where_it_stops():
    r = t.split_taxon("d__Bacteria; p__Bacillota; c__Bacilli")
    assert "genus" not in r and r["class"] == "Bacilli"


# ---- taxonomy.tsv --------------------------------------------------------

def test_types_row_is_not_a_feature():
    got = t.parse_taxonomy(taxonomy_tsv([("f1", GG, 0.99)]), "x")
    assert list(got) == ["f1"] and got["f1"]["confidence"] == 0.99


def test_missing_taxon_column_fires():
    with pytest.raises(t.TaxonomyError, match="no Taxon column"):
        t.parse_taxonomy("Feature ID\tConfidence\nf1\t0.9\n", "x")


def test_wrong_first_column_fires():
    with pytest.raises(t.TaxonomyError, match="expected 'Feature ID'"):
        t.parse_taxonomy("id\tTaxon\nf1\tx\n", "x")


def test_duplicate_feature_fires():
    with pytest.raises(t.TaxonomyError, match="more than once"):
        t.parse_taxonomy(taxonomy_tsv([("f1", GG, 0.9), ("f1", GG, 0.9)]), "x")


def test_missing_confidence_is_not_fatal():
    got = t.parse_taxonomy("Feature ID\tTaxon\nf1\t" + GG + "\n", "x")
    assert got["f1"]["confidence"] is None


# ---- artifacts -----------------------------------------------------------

def test_a_classifier_of_the_wrong_type_fires(tmp_path):
    path = artifact(tmp_path / "wrong.qza", "FeatureData[Sequence]")
    with pytest.raises(t.TaxonomyError, match="not a TaxonomicClassifier"):
        t.check_classifier("gg2", path)


def test_a_real_classifier_passes(tmp_path):
    path = artifact(tmp_path / "ok.qza", "TaxonomicClassifier")
    assert t.check_classifier("gg2", path)["uuid"] == "0000-1111"


def test_a_missing_classifier_fires(tmp_path):
    with pytest.raises(t.TaxonomyError, match="not found or empty"):
        t.check_classifier("gg2", str(tmp_path / "nope.qza"))


def test_a_non_artifact_fires(tmp_path):
    p = tmp_path / "junk.qza"
    p.write_bytes(b"not a zip")
    with pytest.raises(t.TaxonomyError, match="not a readable artifact"):
        t.check_classifier("gg2", str(p))


@pytest.mark.parametrize("bad", ["gg2", "=path", "gg2=", ""])
def test_bad_classifier_spec_fires(bad):
    with pytest.raises(t.TaxonomyError, match="NAME=PATH"):
        t.parse_classifiers([bad])


def test_repeated_classifier_name_fires():
    with pytest.raises(t.TaxonomyError, match="more than once"):
        t.parse_classifiers(["gg2=/a.qza", "gg2=/b.qza"])


# ---- coverage and disagreement ------------------------------------------

def test_coverage_counts_asvs_and_reads():
    tax = t.parse_taxonomy(taxonomy_tsv([("f1", GG, 0.99),
                                         ("f2", "d__Bacteria; p__Bacillota", 0.8)]), "x")
    cov = t.coverage(tax, {"f1": 100.0, "f2": 900.0})
    assert cov["domain"]["asvs"] == 2 and cov["domain"]["read_fraction"] == 1.0
    assert cov["genus"]["asvs"] == 1
    assert cov["genus"]["read_fraction"] == pytest.approx(0.1)


def test_coverage_without_reads_is_still_counted():
    tax = t.parse_taxonomy(taxonomy_tsv([("f1", GG, 0.99)]), "x")
    cov = t.coverage(tax, None)
    assert cov["species"]["asv_fraction"] == 1.0 and cov["species"]["read_fraction"] == 0.0


def test_disagreement_is_counted_only_where_both_named_it():
    a = t.parse_taxonomy(taxonomy_tsv([("f1", GG, 0.9), ("f2", GG, 0.9)]), "a")
    b = t.parse_taxonomy(taxonomy_tsv([("f1", SILVA, 0.9),
                                       ("f2", "d__Bacteria", 0.9)]), "b")
    d = t.disagreements(a, b)
    # same genus, different phylum naming, and only a names f2 below domain
    assert d["genus"] == {"compared": 1, "differ": 0, "differ_fraction": 0.0,
                          "only_first": 1, "only_second": 0, "examples": []}
    assert d["phylum"]["differ"] == 1 and d["phylum"]["compared"] == 1
    assert d["species"]["differ"] == 0  # SILVA says uncultured, which is not a name


# ---- end to end ----------------------------------------------------------

class FakeRunner:
    """Stands in for t.run_cmd. Writes the artifacts and exports QIIME 2 would."""

    def __init__(self, per_classifier=None, rc=0, write=True, biom=None):
        self.per_classifier = per_classifier or {
            "gg2": taxonomy_tsv([("f1", GG, 0.99), ("f2", GG, 0.9)]),
            "silva": taxonomy_tsv([("f1", SILVA, 0.95), ("f2", SILVA, 0.8)]),
        }
        self.rc, self.write = rc, write
        self.biom = biom if biom is not None else "#OTU ID\ts1\nf1\t10.0\nf2\t90.0\n"
        self.calls = []

    def __call__(self, cmd, timeout):
        self.calls.append(cmd)
        assert timeout > 0
        if not self.write:
            return self.rc, "", "Plugin error from feature-classifier\n"
        if "classify-sklearn" in cmd:
            pathlib.Path(cmd[cmd.index("--o-classification") + 1]).write_bytes(b"x")
        elif "export" in cmd:
            src = pathlib.Path(cmd[cmd.index("--input-path") + 1])
            dest = pathlib.Path(cmd[cmd.index("--output-path") + 1])
            dest.mkdir(parents=True, exist_ok=True)
            if src.name.startswith("taxonomy_"):
                name = src.stem.replace("taxonomy_", "")
                (dest / "taxonomy.tsv").write_text(self.per_classifier[name])
            elif src.name.startswith("table"):
                (dest / "feature-table.biom").write_bytes(b"biom")
            else:
                (dest / "dna-sequences.fasta").write_text(">f1\nACGT\n>f2\nTTTT\n")
        else:  # biom convert
            pathlib.Path(cmd[cmd.index("-o") + 1]).write_text(self.biom)
        return self.rc, "", ""


def setup(tmp_path, names=("gg2", "silva")):
    rep = tmp_path / "rep_seqs.qza"
    rep.write_bytes(b"x")
    specs = []
    for n in names:
        specs += ["-c", f"{n}={artifact(tmp_path / (n + '.qza'), 'TaxonomicClassifier')}"]
    return str(rep), specs


def run_main(monkeypatch, tmp_path, runner, extra=(), names=("gg2", "silva")):
    rep, specs = setup(tmp_path, names)
    monkeypatch.setattr(t, "run_cmd", runner)
    return t.main(["-r", rep, "-o", str(tmp_path / "out"), *specs, *extra])


def test_two_classifiers_end_to_end(tmp_path, monkeypatch):
    runner = FakeRunner()
    assert run_main(monkeypatch, tmp_path, runner) == 0
    out = tmp_path / "out"
    cov = list(csv.DictReader((out / "taxonomy_coverage.tsv").open(), delimiter="\t"))
    assert {r["classifier"] for r in cov} == {"gg2", "silva"}
    dis = {r["rank"]: r for r in
           csv.DictReader((out / "taxonomy_disagreements.tsv").open(), delimiter="\t")}
    assert dis["phylum"]["differ"] == "2" and dis["genus"]["differ"] == "0"
    calls = list(csv.DictReader((out / "taxonomy_calls.tsv").open(), delimiter="\t"))
    assert calls[0]["gg2_taxon"] == GG and calls[0]["silva_confidence"] == "0.95"


def test_parameters_are_passed_explicitly(tmp_path, monkeypatch):
    runner = FakeRunner()
    assert run_main(monkeypatch, tmp_path, runner,
                    extra=["--n-jobs", "4", "--confidence", "0.8"]) == 0
    cmd = [c for c in runner.calls if "classify-sklearn" in c][0]
    assert cmd[cmd.index("--p-n-jobs") + 1] == "4"
    assert cmd[cmd.index("--p-confidence") + 1] == "0.8"


def test_reads_weight_coverage_when_a_table_is_given(tmp_path, monkeypatch):
    runner = FakeRunner()
    table = tmp_path / "table.qza"
    table.write_bytes(b"x")
    assert run_main(monkeypatch, tmp_path, runner, extra=["-b", str(table)]) == 0
    cov = {(r["classifier"], r["rank"]): r for r in
           csv.DictReader((tmp_path / "out" / "taxonomy_coverage.tsv").open(), delimiter="\t")}
    assert cov[("gg2", "genus")]["read_fraction"] == "1.0"


def test_a_table_missing_an_asv_fires(tmp_path, monkeypatch):
    runner = FakeRunner(biom="#OTU ID\ts1\nf1\t10.0\n")   # f2 absent
    table = tmp_path / "table.qza"
    table.write_bytes(b"x")
    assert run_main(monkeypatch, tmp_path, runner, extra=["-b", str(table)]) == 1
    assert "not the table" in (tmp_path / "out" / "taxonomy_log.txt").read_text()


def test_an_asv_with_no_taxonomy_row_fires(tmp_path, monkeypatch):
    runner = FakeRunner(per_classifier={"gg2": taxonomy_tsv([("f1", GG, 0.9)])})
    assert run_main(monkeypatch, tmp_path, runner, names=("gg2",)) == 1
    assert "got no row at all" in (tmp_path / "out" / "taxonomy_log.txt").read_text()


def test_nothing_placed_at_domain_fires(tmp_path, monkeypatch):
    runner = FakeRunner(per_classifier={
        "gg2": taxonomy_tsv([("f1", "Unassigned", 0.1), ("f2", "Unassigned", 0.1)])})
    assert run_main(monkeypatch, tmp_path, runner, names=("gg2",)) == 1
    log = (tmp_path / "out" / "taxonomy_log.txt").read_text()
    assert "not one ASV was placed even at domain level" in log


def test_one_classifier_says_so_instead_of_faking_a_comparison(tmp_path, monkeypatch):
    runner = FakeRunner(per_classifier={"gg2": taxonomy_tsv([("f1", GG, 0.9),
                                                             ("f2", GG, 0.9)])})
    assert run_main(monkeypatch, tmp_path, runner, names=("gg2",)) == 0
    out = tmp_path / "out"
    assert not (out / "taxonomy_disagreements.tsv").exists()
    assert "no disagreement rate" in (out / "taxonomy_log.txt").read_text()


def test_the_amplicon_ceiling_is_stated(tmp_path, monkeypatch):
    assert run_main(monkeypatch, tmp_path, FakeRunner()) == 0
    log = (tmp_path / "out" / "taxonomy_log.txt").read_text()
    assert "property of the amplicon" in log


def test_mock_calls_are_listed(tmp_path, monkeypatch):
    targets = tmp_path / "targets.tsv"
    targets.write_text("reference_names\tlength\tsequence\n"
                       "r1;r2\t4\tACGT\n"
                       "r3\t\tNO REGION FOUND BETWEEN THE PRIMERS\n")
    assert run_main(monkeypatch, tmp_path, FakeRunner(),
                    extra=["--mock-targets", str(targets)]) == 0
    rows = list(csv.DictReader((tmp_path / "out" / "taxonomy_mock_calls.tsv").open(),
                               delimiter="\t"))
    assert len(rows) == 1 and rows[0]["reference_names"] == "r1;r2"
    assert rows[0]["feature"] == "f1" and rows[0]["gg2_taxon"] == GG
    assert "listed, not" in (tmp_path / "out" / "taxonomy_log.txt").read_text()


def test_a_targets_file_with_no_sequences_fires(tmp_path, monkeypatch):
    targets = tmp_path / "targets.tsv"
    targets.write_text("reference_names\tlength\tsequence\n")
    assert run_main(monkeypatch, tmp_path, FakeRunner(),
                    extra=["--mock-targets", str(targets)]) == 1
    assert "no target sequences" in (tmp_path / "out" / "taxonomy_log.txt").read_text()


def test_classify_failure_is_reported(tmp_path, monkeypatch):
    assert run_main(monkeypatch, tmp_path, FakeRunner(rc=1, write=False)) == 1
    log = (tmp_path / "out" / "taxonomy_log.txt").read_text()
    assert "Plugin error from feature-classifier" in log


def test_missing_rep_seqs_fires(tmp_path, monkeypatch):
    monkeypatch.setattr(t, "run_cmd", FakeRunner())
    _, specs = setup(tmp_path)
    assert t.main(["-r", str(tmp_path / "nope.qza"), "-o", str(tmp_path / "out"),
                   *specs]) == 1


@pytest.mark.parametrize("bad", [["--n-jobs", "0"], ["--n-jobs", "-2"],
                                 ["--confidence", "1.5"], ["--confidence", "-0.1"]])
def test_out_of_range_options_fire(tmp_path, monkeypatch, bad):
    assert run_main(monkeypatch, tmp_path, FakeRunner(), extra=bad) == 1


# ---- the no-hang rule ----------------------------------------------------

def test_a_stalled_command_is_killed():
    start = time.monotonic()
    with pytest.raises(t.TaxonomyError, match="timed out after 1 s"):
        t.run_cmd(["sh", "-c", "sleep 30"], timeout=1)
    assert time.monotonic() - start < 10
