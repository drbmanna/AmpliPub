"""Offline tests for 09_tree.py. Each guard is shown to fire.

The two that matter: table features missing from the tree, which breaks phylogenetic
diversity much later, and a masked alignment that has collapsed, which builds a tree on
noise without any error.
"""

import csv
import importlib.util
import pathlib
import sys
import time
import zipfile

import pytest

HERE = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("tree", HERE.parent / "09_tree.py")
t = importlib.util.module_from_spec(spec)
sys.modules["tree"] = t
spec.loader.exec_module(t)

NEWICK = "((f1:0.1,f2:0.2)0.95:0.05,(f3:0.3,f4:0.1)0.88:0.02);"


def artifact(path, files):
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("uuid/metadata.yaml", "type: X\n")
        for name, content in files.items():
            zf.writestr(f"uuid/data/{name}", content)
    return str(path)


def alignment_fasta(ids, length):
    return "".join(f">{i}\n{'A' * length}\n" for i in ids)


# ---- Newick ---------------------------------------------------------------

def test_tips_are_extracted():
    assert t.newick_tips(NEWICK) == {"f1", "f2", "f3", "f4"}


def test_fasttree_support_values_are_not_tips():
    """The bug this parser was rewritten for: support values sit exactly where an
    internal label goes, and a naive parse reported them as tips."""
    assert t.newick_tips(NEWICK) == {"f1", "f2", "f3", "f4"}
    assert "0.95" not in t.newick_tips(NEWICK)
    assert "0.88" not in t.newick_tips(NEWICK)


def test_named_internal_nodes_are_not_tips():
    tree = "((a:0.1,b:0.2)clade_x:0.05,(c:0.3,d:0.1)clade_y:0.02)root;"
    assert t.newick_tips(tree) == {"a", "b", "c", "d"}


def test_a_tree_without_branch_lengths_still_parses():
    assert t.newick_tips("((a,b),(c,d));") == {"a", "b", "c", "d"}


def test_quoted_labels_are_unquoted():
    assert "a b" in t.newick_tips("(('a b':0.1,f2:0.2):0.05);")


def test_an_empty_tree_fires():
    with pytest.raises(t.TreeError, match="empty or not Newick"):
        t.newick_tips("   ")


def test_a_non_newick_string_fires():
    with pytest.raises(t.TreeError, match="empty or not Newick"):
        t.newick_tips("f1,f2,f3")


# ---- alignment ------------------------------------------------------------

def test_alignment_length_is_read():
    assert t.fasta_lengths(alignment_fasta(["a", "b"], 100)) == (2, 100)


def test_a_ragged_alignment_fires():
    with pytest.raises(t.TreeError, match="ragged"):
        t.fasta_lengths(">a\nAAAA\n>b\nAA\n")


def test_an_empty_alignment_fires():
    with pytest.raises(t.TreeError, match="holds no sequences"):
        t.fasta_lengths("\n")


def test_wrapped_alignment_lines_are_joined():
    assert t.fasta_lengths(">a\nAAAA\nAAAA\n>b\nCCCC\nCCCC\n") == (2, 8)


# ---- end to end -----------------------------------------------------------

class FakeRunner:
    def __init__(self, tips=NEWICK, aln_len=200, masked_len=150, n_aln=4, n_mask=4,
                 features=("f1", "f2", "f3", "f4"), rc=0, write=True):
        self.tips, self.aln_len, self.masked_len = tips, aln_len, masked_len
        self.n_aln, self.n_mask = n_aln, n_mask
        self.features = features
        self.rc, self.write = rc, write
        self.calls = []

    def __call__(self, cmd, timeout):
        self.calls.append(cmd)
        assert timeout > 0
        if not self.write:
            return self.rc, "", "Plugin error from phylogeny\n"
        if "align-to-tree-mafft-fasttree" in cmd:
            ids = [f"s{i}" for i in range(self.n_aln)]
            artifact(cmd[cmd.index("--o-alignment") + 1],
                     {"aligned-dna-sequences.fasta":
                      alignment_fasta(ids, self.aln_len)})
            artifact(cmd[cmd.index("--o-masked-alignment") + 1],
                     {"aligned-dna-sequences.fasta":
                      alignment_fasta([f"s{i}" for i in range(self.n_mask)],
                                      self.masked_len)})
            artifact(cmd[cmd.index("--o-tree") + 1], {"tree.nwk": self.tips})
            artifact(cmd[cmd.index("--o-rooted-tree") + 1], {"tree.nwk": self.tips})
        elif "export" in cmd:
            dest = pathlib.Path(cmd[cmd.index("--output-path") + 1])
            dest.mkdir(parents=True, exist_ok=True)
            (dest / "feature-table.biom").write_bytes(b"biom")
        else:  # biom convert
            lines = ["# Constructed from biom file", "#OTU ID\ts1"]
            lines += [f"{x}\t1.0" for x in self.features]
            pathlib.Path(cmd[cmd.index("-o") + 1]).write_text("\n".join(lines) + "\n")
        return self.rc, "", ""


def run_main(tmp_path, runner, extra=(), with_table=True, monkeypatch=None):
    (tmp_path / "seqs.qza").write_bytes(b"x")
    (tmp_path / "table.qza").write_bytes(b"x")
    monkeypatch.setattr(t, "run_cmd", runner)
    argv = ["-r", str(tmp_path / "seqs.qza"), "-o", str(tmp_path / "out"), *extra]
    if with_table:
        argv += ["-b", str(tmp_path / "table.qza")]
    return t.main(argv)


def test_end_to_end(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(), monkeypatch=monkeypatch) == 0
    out = tmp_path / "out"
    for name in ("alignment.qza", "masked_alignment.qza", "unrooted_tree.qza",
                 "rooted_tree.qza", "tree_tips.tsv"):
        assert (out / name).exists()
    log = (out / "tree_log.txt").read_text()
    assert "masking left 150 columns (75.0%)" in log
    assert "every one of the 4 table features is a tip" in log


def test_threads_and_mask_settings_are_explicit(tmp_path, monkeypatch):
    runner = FakeRunner()
    assert run_main(tmp_path, runner, extra=["--threads", "4"],
                    monkeypatch=monkeypatch) == 0
    cmd = runner.calls[0]
    assert cmd[cmd.index("--p-n-threads") + 1] == "4"
    assert cmd[cmd.index("--p-mask-max-gap-frequency") + 1] == "1.0"
    assert cmd[cmd.index("--p-mask-min-conservation") + 1] == "0.4"


def test_a_feature_missing_from_the_tree_is_fatal(tmp_path, monkeypatch):
    """The guard this stage exists for."""
    runner = FakeRunner(features=("f1", "f2", "f3", "f4", "f5"))
    assert run_main(tmp_path, runner, monkeypatch=monkeypatch) == 1
    log = (tmp_path / "out" / "tree_log.txt").read_text()
    assert "not tips in the tree" in log
    assert "build the tree from the filtered sequences instead" in log


def test_extra_tips_only_warn(tmp_path, monkeypatch):
    runner = FakeRunner(features=("f1", "f2"))
    assert run_main(tmp_path, runner, monkeypatch=monkeypatch) == 0
    assert "tip(s) are not in the table" in (tmp_path / "out" / "tree_log.txt").read_text()


def test_a_collapsed_masked_alignment_is_fatal(tmp_path, monkeypatch):
    runner = FakeRunner(aln_len=1000, masked_len=100)
    assert run_main(tmp_path, runner, monkeypatch=monkeypatch) == 1
    log = (tmp_path / "out" / "tree_log.txt").read_text()
    assert "masking left only 10.0% of the alignment" in log
    assert "would be noise" in log


def test_the_collapse_floor_can_be_lowered(tmp_path, monkeypatch):
    runner = FakeRunner(aln_len=1000, masked_len=100)
    assert run_main(tmp_path, runner, extra=["--min-masked-fraction", "0.05"],
                    monkeypatch=monkeypatch) == 0


def test_a_sequence_lost_in_masking_is_fatal(tmp_path, monkeypatch):
    runner = FakeRunner(n_aln=4, n_mask=3)
    assert run_main(tmp_path, runner, monkeypatch=monkeypatch) == 1
    assert "masked alignment holds 3" in (tmp_path / "out" / "tree_log.txt").read_text()


def test_running_without_a_table_warns_loudly(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(), with_table=False,
                    monkeypatch=monkeypatch) == 0
    log = (tmp_path / "out" / "tree_log.txt").read_text()
    assert "That check is the reason this stage exists" in log


def test_tips_table_records_both_sides(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(features=("f1", "f2")),
                    monkeypatch=monkeypatch) == 0
    rows = {r["id"]: r for r in
            csv.DictReader((tmp_path / "out" / "tree_tips.tsv").open(), delimiter="\t")}
    assert rows["f1"] == {"id": "f1", "in_tree": "yes", "in_table": "yes"}
    assert rows["f3"]["in_table"] == "no"


def test_a_qiime_failure_is_reported(tmp_path, monkeypatch):
    assert run_main(tmp_path, FakeRunner(rc=1, write=False), monkeypatch=monkeypatch) == 1
    assert "Plugin error from phylogeny" in (tmp_path / "out" / "tree_log.txt").read_text()


def test_a_missing_output_fires(tmp_path, monkeypatch):
    class NoTree(FakeRunner):
        def __call__(self, cmd, timeout):
            super().__call__(cmd, timeout)
            if "align-to-tree-mafft-fasttree" in cmd:
                pathlib.Path(cmd[cmd.index("--o-rooted-tree") + 1]).unlink()
            return 0, "", ""
    assert run_main(tmp_path, NoTree(), monkeypatch=monkeypatch) == 1
    assert "wrote no" in (tmp_path / "out" / "tree_log.txt").read_text()


@pytest.mark.parametrize("bad", [["--threads", "0"],
                                 ["--mask-min-conservation", "1.5"],
                                 ["--min-masked-fraction", "-0.1"]])
def test_out_of_range_options_fire(tmp_path, monkeypatch, bad):
    assert run_main(tmp_path, FakeRunner(), extra=bad, monkeypatch=monkeypatch) == 1


def test_missing_rep_seqs_fires(tmp_path, monkeypatch):
    monkeypatch.setattr(t, "run_cmd", FakeRunner())
    assert t.main(["-r", str(tmp_path / "nope.qza"), "-o", str(tmp_path / "out")]) == 1


# ---- the no-hang rule -----------------------------------------------------

def test_a_stalled_command_is_killed():
    start = time.monotonic()
    with pytest.raises(t.TreeError, match="timed out after 1 s"):
        t.run_cmd(["sh", "-c", "sleep 30"], timeout=1)
    assert time.monotonic() - start < 10
