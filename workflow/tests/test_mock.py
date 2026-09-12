"""Offline tests for 04_mock.py. Each guard is shown to fire.

The reference here is synthetic and built to exercise the awkward cases: two records
sharing one region, a record whose forward primer site carries a mismatch, and a record
with no reverse site at all. The real HMP mock reference is not in this repository, so
the check against it is recorded in the project notes rather than run here.
"""

import csv
import importlib.util
import pathlib
import sys
import time

import pytest

HERE = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("mock", HERE.parent / "04_mock.py")
m = importlib.util.module_from_spec(spec)
sys.modules["mock"] = m
spec.loader.exec_module(m)

def concrete(primer):
    """One real sequence satisfying an IUPAC primer, for building fixtures."""
    return "".join(m.IUPAC[c][0] for c in primer)


F = m.PRIMER_F                    # GTGCCAGCMGCCGCGGTAA, with an ambiguous M
F_SEQ = concrete(F)               # what a real template carries at that site
R_RC = m.revcomp(concrete(m.PRIMER_R))   # 806R site as it appears on the sense strand
FLANK = "AAGGTTCCAA"

V4_A = "ACGT" * 20                 # 80 bp
V4_B = "TTGCA" * 16                # 80 bp
V4_C = "GGCCT" * 16                # 80 bp, only reachable through a relaxed primer site


def record(v4, fwd=F_SEQ, rev_rc=R_RC, tail=True):
    return FLANK + fwd + v4 + (rev_rc if tail else "") + FLANK


def reference(tmp_path, extra=()):
    """Two records share V4_A; one needs a relaxed site; one has no reverse site."""
    relaxed_f = "A" + F_SEQ[1:]    # one mismatch at the first position
    records = [("r1", record(V4_A)), ("r2", record(V4_A)), ("r3", record(V4_B)),
               ("r4", record(V4_C, fwd=relaxed_f)), ("r5", record(V4_B, tail=False))]
    records.extend(extra)
    p = tmp_path / "ref.fasta"
    p.write_text("".join(f">{n}\n{s}\n" for n, s in records))
    return str(p)


# ---- sequence helpers ----------------------------------------------------

def test_revcomp_handles_iupac():
    assert m.revcomp("ACGT") == "ACGT"
    assert m.revcomp("GGACTACHVGGGTWTCTAAT") == "ATTAGAWACCCBDGTAGTCC"


def test_revcomp_rejects_junk():
    with pytest.raises(m.MockError, match="cannot complement"):
        m.revcomp("ACGZ")


def test_iupac_mismatches_are_counted():
    assert m.mismatches("ACGT", "ACGT") == 0
    assert m.mismatches("AMGT", "ACGT") == 0   # M is A or C
    assert m.mismatches("AMGT", "AGGT") == 1   # G is neither
    assert m.mismatches("NNNN", "ACGT") == 0


def test_unknown_primer_code_fires():
    with pytest.raises(m.MockError, match="not a IUPAC code"):
        m.mismatches("ACGZ", "ACGT")


def test_exact_match_wins_over_an_earlier_near_match():
    primer = "ACGTACGT"
    near = "ACGTACGA"                       # one mismatch, comes first
    seq = "TT" + near + "TT" + primer + "TT"
    assert m.find_primer(seq, primer, 1) == (12, 20)


def test_near_match_is_used_when_nothing_is_exact():
    primer = "ACGTACGT"
    seq = "TT" + "ACGTACGA" + "TT"
    assert m.find_primer(seq, primer, 1) == (2, 10)
    assert m.find_primer(seq, primer, 0) is None


def test_region_is_what_lies_between_the_primers(tmp_path):
    assert m.extract_region(record(V4_A), F, m.PRIMER_R, 0) == V4_A


def test_region_is_none_without_a_reverse_site():
    assert m.extract_region(record(V4_B, tail=False), F, m.PRIMER_R, 0) is None


def test_targets_collapse_and_report_what_is_missing(tmp_path):
    records = m.read_fasta(reference(tmp_path))
    targets, missing = m.build_targets(records, F, m.PRIMER_R, 0)
    assert len(targets) == 2                       # V4_A (shared) and V4_B
    assert targets[V4_A] == ["r1", "r2"]
    assert sorted(missing) == ["r4", "r5"]


def test_a_relaxed_primer_site_rescues_a_record(tmp_path):
    records = m.read_fasta(reference(tmp_path))
    targets, missing = m.build_targets(records, F, m.PRIMER_R, 1)
    assert len(targets) == 3 and V4_C in targets
    assert missing == ["r5"]                       # no reverse site at any budget


def test_no_region_anywhere_fires(tmp_path):
    p = tmp_path / "bad.fasta"
    p.write_text(">x\nACGTACGTACGTACGT\n")
    with pytest.raises(m.MockError, match="no reference record yielded"):
        m.build_targets(m.read_fasta(str(p)), F, m.PRIMER_R, 0)


# ---- FASTA reading -------------------------------------------------------

def test_duplicate_record_names_fire(tmp_path):
    p = tmp_path / "dup.fasta"
    p.write_text(">a\nACGT\n>a\nTGCA\n")
    with pytest.raises(m.MockError, match="more than once"):
        m.read_fasta(str(p))


def test_sequence_before_a_header_fires(tmp_path):
    p = tmp_path / "head.fasta"
    p.write_text("ACGT\n>a\nACGT\n")
    with pytest.raises(m.MockError, match="before the first header"):
        m.read_fasta(str(p))


def test_empty_fasta_fires(tmp_path):
    p = tmp_path / "empty.fasta"
    p.write_text("\n\n")
    with pytest.raises(m.MockError, match="no sequences"):
        m.read_fasta(str(p))


def test_wrapped_sequence_is_joined(tmp_path):
    p = tmp_path / "wrap.fasta"
    p.write_text(">a\nACGT\nTTTT\n")
    assert m.read_fasta(str(p)) == {"a": "ACGTTTTT"}


# ---- distances -----------------------------------------------------------

def test_hamming_needs_equal_lengths():
    with pytest.raises(m.MockError, match="equal lengths"):
        m.hamming("ACGT", "ACG")


def test_nearest_ignores_other_lengths():
    assert m.nearest_same_length("ACGT", ["ACGA", "ACG"]) == (1, 4)
    assert m.nearest_same_length("ACGT", ["ACG"]) is None


# ---- the exported table --------------------------------------------------

BIOM = ("# Constructed from biom file\n"
        "#OTU ID\tmock1\tmock2\n"
        "f1\t100.0\t0.0\n"
        "f2\t5.0\t50.0\n")


def test_biom_tsv_is_read_and_zeros_dropped():
    got = m.parse_biom_tsv(BIOM, "x")
    assert got == {"mock1": {"f1": 100.0, "f2": 5.0}, "mock2": {"f2": 50.0}}


def test_wrong_first_column_fires():
    with pytest.raises(m.MockError, match="expected '#OTU ID'"):
        m.parse_biom_tsv("feature\ts1\nf1\t1\n", "x")


def test_ragged_row_fires():
    with pytest.raises(m.MockError, match="cells for"):
        m.parse_biom_tsv("#OTU ID\ts1\ts2\nf1\t1\n", "x")


def test_negative_count_fires():
    with pytest.raises(m.MockError, match="negative"):
        m.parse_biom_tsv("#OTU ID\ts1\nf1\t-1\n", "x")


def test_non_numeric_count_fires():
    with pytest.raises(m.MockError, match="not a number"):
        m.parse_biom_tsv("#OTU ID\ts1\nf1\tmany\n", "x")


def test_table_without_samples_fires():
    with pytest.raises(m.MockError, match="no samples"):
        m.parse_biom_tsv("#OTU ID\nf1\n", "x")


# ---- scoring one sample --------------------------------------------------

def test_exact_and_spurious_are_separated():
    targets = {V4_A: ["r1"], V4_B: ["r3"]}
    off_by_two = "TT" + V4_A[2:]
    seqs = {"f1": V4_A, "f2": off_by_two, "f3": V4_A[:-1]}   # f3 is one base short
    counts = {"f1": 900.0, "f2": 100.0, "f3": 10.0}
    r = m.score_sample(counts, seqs, targets, [3])
    assert r["targets_recovered"] == 1 and r["targets_total"] == 2
    assert r["other_asvs"] == 2
    assert r["other_read_fraction"] == pytest.approx(110 / 1010)
    assert r["asvs_without_same_length_reference"] == 1   # f3
    # f2 is 2 mismatches over 80 bases in 100 reads; the 900 exact reads are zero-error
    # and belong in the denominator, so the rate is 200 / (8000 + 72000)
    assert r["within"][3]["mismatch_rate"] == pytest.approx(200 / 80000)
    assert r["within"][3]["mismatch_rate_variants_only"] == pytest.approx(2 / 80)
    assert r["unattributable_asvs"] == 1 and r["unattributable_reads"] == 10.0


def test_a_foreign_organism_is_excluded_from_the_error_rate():
    """The bug this metric was rebuilt to avoid: 40 mismatches is a different organism."""
    targets = {V4_A: ["r1"]}
    foreign = V4_B                                   # unrelated, far from V4_A
    one_off = "T" + V4_A[1:] if V4_A[0] != "T" else "A" + V4_A[1:]
    seqs = {"f1": V4_A, "f2": one_off, "f3": foreign}
    r = m.score_sample({"f1": 100.0, "f2": 10.0, "f3": 890.0}, seqs, targets, [1, 3, 10])
    assert m.hamming(foreign, V4_A) > 10             # the fixture really is far away
    assert r["within"][1]["asvs"] == 1 and r["within"][10]["asvs"] == 1
    # 1 error over 80 bases in 10 reads, against those reads plus the 100 exact ones
    assert r["within"][3]["mismatch_rate"] == pytest.approx(10 / (800 + 8000))
    assert r["unattributable_read_fraction"] == pytest.approx(0.89)


def test_cutoffs_are_parsed_and_checked():
    assert m.read_cutoffs("10,1,3,3") == [1, 3, 10]
    for bad in ("", "0,3", "-1", "a,b"):
        with pytest.raises(m.MockError):
            m.read_cutoffs(bad)


def test_a_feature_missing_from_the_sequences_fires():
    with pytest.raises(m.MockError, match="not in the sequences"):
        m.score_sample({"f9": 1.0}, {}, {V4_A: ["r1"]}, [3])


def test_two_references_sharing_a_region_count_once():
    targets = {V4_A: ["r1", "r2"]}
    r = m.score_sample({"f1": 10.0}, {"f1": V4_A}, targets, [3])
    assert r["targets_recovered"] == 1 and r["targets_total"] == 1


# ---- the sample list -----------------------------------------------------

def test_sample_list_from_a_comma_string():
    assert m.read_sample_list("mock1, mock2 ,mock1") == ["mock1", "mock2"]


def test_sample_list_from_a_file(tmp_path):
    p = tmp_path / "mocks.txt"
    p.write_text("# mocks\nmock1\n\nmock2 # note\n")
    assert m.read_sample_list(str(p)) == ["mock1", "mock2"]


def test_empty_sample_list_fires():
    with pytest.raises(m.MockError, match="no mock sample ids"):
        m.read_sample_list(" , ")


# ---- end to end ----------------------------------------------------------

class FakeRunner:
    """Stands in for m.run_cmd. Writes what qiime tools export and biom convert write."""

    def __init__(self, biom=BIOM, seqs=None, rc=0, write=True):
        self.biom = biom
        self.seqs = seqs if seqs is not None else {"f1": V4_A, "f2": V4_B}
        self.rc, self.write = rc, write
        self.calls = []

    def __call__(self, cmd, timeout):
        self.calls.append(cmd)
        assert timeout > 0
        if not self.write:
            return self.rc, "", "Plugin error\n"
        if "export" in cmd:
            dest = pathlib.Path(cmd[cmd.index("--output-path") + 1])
            dest.mkdir(parents=True, exist_ok=True)
            if dest.name == "table":
                (dest / "feature-table.biom").write_bytes(b"biom")
            else:
                (dest / "dna-sequences.fasta").write_text(
                    "".join(f">{k}\n{v}\n" for k, v in self.seqs.items()))
        else:  # biom convert
            pathlib.Path(cmd[cmd.index("-o") + 1]).write_text(self.biom)
        return self.rc, "", ""


def run_main(monkeypatch, tmp_path, runner, extra=(), samples="mock1,mock2"):
    for name in ("table.qza", "rep_seqs.qza"):
        (tmp_path / name).write_bytes(b"x")
    monkeypatch.setattr(m, "run_cmd", runner)
    return m.main(["-b", str(tmp_path / "table.qza"), "-r", str(tmp_path / "rep_seqs.qza"),
                   "-m", reference(tmp_path), "-s", samples,
                   "-o", str(tmp_path / "out"), *extra])


def test_end_to_end_writes_the_report(tmp_path, monkeypatch):
    assert run_main(monkeypatch, tmp_path, FakeRunner()) == 0
    out = tmp_path / "out"
    rows = {r["sample-id"]: r for r in
            csv.DictReader((out / "mock_summary.tsv").open(), delimiter="\t")}
    assert rows["mock1"]["reads"] == "105.0" and rows["mock1"]["targets_recovered"] == "2"
    assert rows["mock2"]["targets_recovered"] == "1"
    targets = (out / "targets.tsv").read_text()
    assert "r1;r2" in targets and "NO REGION FOUND" in targets
    missing = list(csv.DictReader((out / "mock_missing_targets.tsv").open(), delimiter="\t"))
    # mock2 holds only V4_B, so both of the other targets are missing from it
    assert sorted(r["reference_names"] for r in missing
                  if r["sample-id"] == "mock2") == ["r1;r2", "r4"]
    assert "does not pass or fail" in (out / "mock_log.txt").read_text()


def test_relaxed_site_is_reported_in_the_log(tmp_path, monkeypatch):
    assert run_main(monkeypatch, tmp_path, FakeRunner(),
                    extra=["--max-primer-mismatch", "1"]) == 0
    log = (tmp_path / "out" / "mock_log.txt").read_text()
    assert "only found with a relaxed primer site: r4" in log
    assert "no region between the primers in: r5" in log


def test_published_figures_are_quoted_not_applied(tmp_path, monkeypatch):
    assert run_main(monkeypatch, tmp_path, FakeRunner()) == 0
    log = (tmp_path / "out" / "mock_log.txt").read_text()
    assert "not as a threshold" in log and "Kozich" in log
    assert "not the same unit" in log


def test_low_depth_note_is_opt_in(tmp_path, monkeypatch):
    assert run_main(monkeypatch, tmp_path, FakeRunner()) == 0
    assert "below the" not in (tmp_path / "out" / "mock_log.txt").read_text()


def test_low_depth_note_fires_when_asked(tmp_path, monkeypatch):
    assert run_main(monkeypatch, tmp_path, FakeRunner(),
                    extra=["--low-depth-note", "200"]) == 0
    assert "below the 200 you asked" in (tmp_path / "out" / "mock_log.txt").read_text()


def test_unknown_mock_sample_fires(tmp_path, monkeypatch):
    assert run_main(monkeypatch, tmp_path, FakeRunner(), samples="mock9") == 1
    assert "not in the table" in (tmp_path / "out" / "mock_log.txt").read_text()


def test_export_failure_is_reported(tmp_path, monkeypatch):
    assert run_main(monkeypatch, tmp_path, FakeRunner(rc=1, write=False)) == 1
    assert "Plugin error" in (tmp_path / "out" / "mock_log.txt").read_text()


def test_missing_input_fires(tmp_path, monkeypatch):
    monkeypatch.setattr(m, "run_cmd", FakeRunner())
    assert m.main(["-b", str(tmp_path / "nope.qza"), "-r", str(tmp_path / "nope.qza"),
                   "-m", reference(tmp_path), "-s", "mock1",
                   "-o", str(tmp_path / "out")]) == 1


def test_negative_mismatch_budget_fires(tmp_path, monkeypatch):
    assert run_main(monkeypatch, tmp_path, FakeRunner(),
                    extra=["--max-primer-mismatch", "-1"]) == 1


# ---- the no-hang rule ----------------------------------------------------

def test_a_stalled_command_is_killed():
    start = time.monotonic()
    with pytest.raises(m.MockError, match="timed out after 1 s"):
        m.run_cmd(["sh", "-c", "sleep 30"], timeout=1)
    assert time.monotonic() - start < 10
