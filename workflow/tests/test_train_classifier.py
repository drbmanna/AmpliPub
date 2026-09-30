"""Offline tests for scripts/train_classifier.py. Each guard is shown to fire.

QIIME 2 is replaced by a fake runner that writes the artifacts the real commands would
write, so the order of the two steps, the type checks on both inputs and on our own output,
and the provenance record are all exercised without QIIME 2 installed.
"""

import importlib.util
import pathlib
import sys
import zipfile

import pytest

HERE = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "train_classifier", HERE.parent / "scripts" / "train_classifier.py")
t = importlib.util.module_from_spec(spec)
sys.modules["train_classifier"] = t
spec.loader.exec_module(t)

FWD = "CCTACGGGNGGCWGCAG"        # 341F
REV = "GACTACHVGGGTATCTAATCC"    # 805R


def make_qza(path, kind, uuid="0000-1111"):
    """A minimal artifact: enough metadata.yaml for artifact_meta to read a type."""
    path = pathlib.Path(path)
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr(f"{uuid}/metadata.yaml", f"uuid: {uuid}\ntype: {kind}\nformat: whatever\n")
    return str(path)


class FakeQiime:
    """Stands in for t.run_cmd. Writes whatever --o-reads or --o-classifier names."""

    def __init__(self, reads_type=t.SEQ_TYPE, classifier_type=t.CLS_TYPE, rc=0, write=True):
        self.reads_type, self.classifier_type = reads_type, classifier_type
        self.rc, self.write = rc, write
        self.calls = []

    def __call__(self, cmd, timeout):
        self.calls.append(cmd)
        assert timeout > 0
        if self.write:
            for flag, kind in (("--o-reads", self.reads_type),
                               ("--o-classifier", self.classifier_type)):
                if flag in cmd:
                    make_qza(cmd[cmd.index(flag) + 1], kind)
        return self.rc, "", "Plugin error from feature-classifier\n" if self.rc else ""

    @property
    def subcommands(self):
        return [c[c.index("qiime") + 1:c.index("qiime") + 3] for c in self.calls]


def run_main(monkeypatch, tmp_path, runner, seqs=None, tax=None, extra=()):
    seqs = seqs or make_qza(tmp_path / "seqs.qza", t.SEQ_TYPE)
    tax = tax or make_qza(tmp_path / "tax.qza", t.TAX_TYPE)
    monkeypatch.setattr(t, "run_cmd", runner)
    return t.main(["--sequences", seqs, "--taxonomy", tax,
                   "--forward", FWD, "--reverse", REV,
                   "-o", str(tmp_path / "out"), "--name", "gg2_test_v3v4", *extra])


# Input guards, all before any compute

def test_missing_sequences_are_refused(monkeypatch, tmp_path, caplog):
    runner = FakeQiime()
    assert run_main(monkeypatch, tmp_path, runner,
                    seqs=str(tmp_path / "nope.qza")) == 1
    assert "not found or empty" in caplog.text
    assert runner.calls == []


def test_a_file_that_is_not_an_artifact_is_refused(monkeypatch, tmp_path, caplog):
    bad = tmp_path / "bad.qza"
    bad.write_bytes(b"not a zip")
    runner = FakeQiime()
    assert run_main(monkeypatch, tmp_path, runner, seqs=str(bad)) == 1
    assert "not a readable artifact" in caplog.text
    assert runner.calls == []


def test_sequences_of_the_wrong_type_are_refused(monkeypatch, tmp_path, caplog):
    wrong = make_qza(tmp_path / "wrong.qza", "FeatureTable[Frequency]")
    runner = FakeQiime()
    assert run_main(monkeypatch, tmp_path, runner, seqs=wrong) == 1
    assert "not a FeatureData[Sequence]" in caplog.text
    assert runner.calls == []


def test_taxonomy_of_the_wrong_type_is_refused(monkeypatch, tmp_path, caplog):
    wrong = make_qza(tmp_path / "wrong.qza", t.SEQ_TYPE)
    runner = FakeQiime()
    assert run_main(monkeypatch, tmp_path, runner, tax=wrong) == 1
    assert "not a FeatureData[Taxonomy]" in caplog.text
    assert runner.calls == []


def test_an_existing_classifier_is_never_silently_replaced(monkeypatch, tmp_path, caplog):
    out = tmp_path / "out"
    out.mkdir()
    make_qza(out / "gg2_test_v3v4.qza", t.CLS_TYPE)
    runner = FakeQiime()
    assert run_main(monkeypatch, tmp_path, runner) == 1
    assert "already exists" in caplog.text
    assert runner.calls == []


# Failure of the qiime calls themselves

def test_a_failing_qiime_call_is_fatal(monkeypatch, tmp_path, caplog):
    assert run_main(monkeypatch, tmp_path, FakeQiime(rc=1)) == 1
    assert "exited with code 1" in caplog.text


def test_a_step_that_writes_nothing_is_fatal(monkeypatch, tmp_path, caplog):
    assert run_main(monkeypatch, tmp_path, FakeQiime(write=False)) == 1
    assert "wrote no" in caplog.text


def test_our_own_output_is_checked_not_assumed(monkeypatch, tmp_path, caplog):
    """A file of the right size and the wrong type must fail here, not hours later."""
    runner = FakeQiime(classifier_type="FeatureData[Sequence]")
    assert run_main(monkeypatch, tmp_path, runner) == 1
    assert "not a TaxonomicClassifier" in caplog.text


# The happy path

def test_both_steps_run_in_order_with_the_given_primers(monkeypatch, tmp_path):
    runner = FakeQiime()
    assert run_main(monkeypatch, tmp_path, runner) == 0
    assert runner.subcommands == [["feature-classifier", "extract-reads"],
                                  ["feature-classifier", "fit-classifier-naive-bayes"]]
    extract = runner.calls[0]
    assert extract[extract.index("--p-f-primer") + 1] == FWD
    assert extract[extract.index("--p-r-primer") + 1] == REV
    fit = runner.calls[1]
    # Training must consume what extraction produced, not the untrimmed reference.
    assert fit[fit.index("--i-reference-reads") + 1] == extract[extract.index("--o-reads") + 1]


def test_provenance_records_the_primers_and_the_checksum(monkeypatch, tmp_path):
    assert run_main(monkeypatch, tmp_path, FakeQiime()) == 0
    prov = tmp_path / "out" / "gg2_test_v3v4_provenance.tsv"
    values = dict(line.rstrip("\n").split("\t", 1) for line in open(prov))
    assert values["forward_primer"] == FWD and values["reverse_primer"] == REV
    assert values["name"] == "gg2_test_v3v4"
    written = tmp_path / "out" / "gg2_test_v3v4.qza"
    assert values["classifier_sha256"] == t.sha256_of(str(written))
    assert len(values["classifier_sha256"]) == 64
    for key in ("reference_sequences_sha256", "reference_taxonomy_sha256", "command", "env"):
        assert values[key]


def test_the_intermediate_reads_go_unless_asked_for(monkeypatch, tmp_path):
    assert run_main(monkeypatch, tmp_path, FakeQiime()) == 0
    assert not (tmp_path / "out" / "gg2_test_v3v4_ref_reads.qza").exists()


def test_keep_reads_keeps_them(monkeypatch, tmp_path):
    assert run_main(monkeypatch, tmp_path, FakeQiime(), extra=("--keep-reads",)) == 0
    assert (tmp_path / "out" / "gg2_test_v3v4_ref_reads.qza").exists()


def test_the_config_block_to_paste_is_printed(monkeypatch, tmp_path, caplog):
    assert run_main(monkeypatch, tmp_path, FakeQiime()) == 0
    assert "put this in the config:" in caplog.text
    assert "name: gg2_test_v3v4" in caplog.text
