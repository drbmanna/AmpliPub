"""Offline tests for 06_resolution.py. Each guard is shown to fire.

The reference is synthetic and built so the right answer is known by construction: two
records share a short region but differ outside it, so a short amplicon collapses them
and a long one does not. That is the whole claim of this stage, so it is the first test.
"""

import csv
import importlib.util
import pathlib
import sys

import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
spec = importlib.util.spec_from_file_location("resolution", HERE.parent / "06_resolution.py")
r = importlib.util.module_from_spec(spec)
sys.modules["resolution"] = r
spec.loader.exec_module(r)

import amplicon_regions as ar  # noqa: E402

F = "GTGCCAGCMGCCGCGGTAA"        # 515F, with an ambiguous M
R = "GGACTACHVGGGTWTCTAAT"       # 806R


def concrete(primer):
    return "".join(ar.IUPAC[c][0] for c in primer)


F_SEQ = concrete(F)
R_RC = ar.revcomp(concrete(R))
SHARED = "ACGT" * 20             # the region two records share
UNIQUE = "TTGCA" * 16


def record(inner, tail_tag):
    """Same inner region, different sequence outside it."""
    return "AAAA" + F_SEQ + inner + R_RC + tail_tag


def reference(tmp_path, name="ref.fasta"):
    """s1 and s2 collapse inside the region and differ outside it. s3 is distinct.
    s4 has no reverse primer site at all."""
    records = [("s1", record(SHARED, "GGGGGGGGGG")),
               ("s2", record(SHARED, "CCCCCCCCCC")),
               ("s3", record(UNIQUE, "AAAAAAAAAA")),
               ("s4", "AAAA" + F_SEQ + UNIQUE + "TTTTTTTT")]
    p = tmp_path / name
    p.write_text("".join(f">{n}\n{s}\n" for n, s in records))
    return str(p)


def taxonomy(tmp_path, name="tax.tsv"):
    p = tmp_path / name
    p.write_text(
        "Feature ID\tTaxon\n"
        "s1\td__Bacteria; g__Staphylococcus; s__Staphylococcus aureus\n"
        "s2\td__Bacteria; g__Staphylococcus; s__Staphylococcus epidermidis\n"
        "s3\td__Bacteria; g__Escherichia; s__Escherichia coli\n"
        "s4\td__Bacteria; g__Bacillus; s__Bacillus cereus\n")
    return str(p)


# ---- the central claim ---------------------------------------------------

def test_a_short_region_collapses_what_full_length_separates(tmp_path):
    records = ar.read_fasta(reference(tmp_path))
    short = r.analyse(records, F_SEQ, concrete(R), 0, None)
    full = r.analyse(records, None, None, 0, None)
    # s1 and s2 are identical inside the region and different outside it
    assert short["collapse_groups"] == 1 and short["records_collapsed"] == 2
    assert full["collapse_groups"] == 0
    assert full["distinct_sequences"] == 4


def test_species_resolution_is_counted_per_rank(tmp_path):
    records = ar.read_fasta(reference(tmp_path))
    tax = r.read_taxonomy(taxonomy(tmp_path))
    got = r.analyse(records, F_SEQ, concrete(R), 0, tax)
    sp = got["ranks"]["species"]
    # S. aureus and S. epidermidis collapse; E. coli does not; B. cereus has no region
    assert sp["labels"] == 3 and sp["resolved"] == 1 and sp["unresolved"] == 2
    assert len(sp["ambiguity_sets"]) == 1
    assert sp["ambiguity_sets"][0]["labels"] == ["Staphylococcus aureus",
                                                 "Staphylococcus epidermidis"]
    # both belong to one genus, so genus is not ambiguous at all
    assert got["ranks"]["genus"]["unresolved"] == 0


def test_a_record_without_a_region_is_reported_not_dropped(tmp_path):
    records = ar.read_fasta(reference(tmp_path))
    got = r.analyse(records, F_SEQ, concrete(R), 0, None)
    assert got["no_region"] == 1 and got["no_region_names"] == ["s4"]
    assert got["with_region"] == 3


# ---- taxon strings -------------------------------------------------------

def test_prefixed_taxon_splits_by_rank():
    got = r.split_taxon("d__Bacteria; p__Bacillota; g__Blautia; s__Blautia wexlerae")
    assert got["domain"] == "Bacteria" and got["species"] == "Blautia wexlerae"
    assert "class" not in got


def test_placeholder_labels_do_not_count():
    got = r.split_taxon("d__Bacteria; p__; g__uncultured; s__Unassigned")
    assert got == {"domain": "Bacteria"}


def test_k_prefix_is_domain():
    assert r.split_taxon("k__Bacteria")["domain"] == "Bacteria"


def test_unprefixed_taxon_falls_back_to_position():
    got = r.split_taxon("Bacteria; Firmicutes")
    assert got == {"domain": "Bacteria", "phylum": "Firmicutes"}


# ---- inputs --------------------------------------------------------------

def test_taxonomy_header_and_types_row_are_skipped(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("Feature ID\tTaxon\n#q2:types\tcategorical\ns1\td__Bacteria\n")
    assert list(r.read_taxonomy(str(p))) == ["s1"]


def test_duplicate_taxonomy_row_fires(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("s1\td__Bacteria\ns1\td__Archaea\n")
    with pytest.raises(r.ResolutionError, match="more than once"):
        r.read_taxonomy(str(p))


def test_empty_taxonomy_fires(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("Feature ID\tTaxon\n")
    with pytest.raises(r.ResolutionError, match="no taxonomy rows"):
        r.read_taxonomy(str(p))


@pytest.mark.parametrize("bad", ["v4", "v4=ACGT", "=ACGT,ACGT", "v4=,ACGT", "v4=ACGT,"])
def test_bad_primer_spec_fires(bad):
    with pytest.raises(r.ResolutionError, match="NAME=FORWARD,REVERSE"):
        r.parse_primers([bad])


def test_repeated_region_name_fires():
    with pytest.raises(r.ResolutionError, match="more than once"):
        r.parse_primers(["v4=AAAA,CCCC", "v4=GGGG,TTTT"])


def test_primers_are_upper_cased():
    assert r.parse_primers(["v4=acgt,tgca"]) == {"v4": ("ACGT", "TGCA")}


# ---- end to end ----------------------------------------------------------

def run_main(tmp_path, extra):
    return r.main(["-m", reference(tmp_path), "-o", str(tmp_path / "out"), *extra])


def test_end_to_end_compares_two_regions(tmp_path):
    assert run_main(tmp_path, ["-p", f"v4={F},{R}", "--full-length",
                               "-t", taxonomy(tmp_path)]) == 0
    out = tmp_path / "out"
    rows = {x["region"]: x for x in
            csv.DictReader((out / "resolution_summary.tsv").open(), delimiter="\t")}
    assert rows["v4"]["species_resolved"] == "1"
    assert rows["full_length"]["species_resolved"] == "4"
    assert rows["v4"]["no_region"] == "1" and rows["full_length"]["no_region"] == "0"

    sets = list(csv.DictReader((out / "resolution_ambiguity_sets.tsv").open(), delimiter="\t"))
    species = [s for s in sets if s["rank"] == "species" and s["region"] == "v4"]
    assert len(species) == 1
    assert "Staphylococcus aureus" in species[0]["labels"]
    assert species[0]["records"] == "s1; s2"

    groups = list(csv.DictReader((out / "resolution_groups.tsv").open(), delimiter="\t"))
    assert [g["records"] for g in groups if g["region"] == "v4"] == ["s1; s2"]
    missing = list(csv.DictReader((out / "resolution_no_region.tsv").open(), delimiter="\t"))
    assert [m["record"] for m in missing] == ["s4"]


def test_it_works_without_taxonomy(tmp_path):
    assert run_main(tmp_path, ["-p", f"v4={F},{R}"]) == 0
    out = tmp_path / "out"
    assert not (out / "resolution_ambiguity_sets.tsv").exists()
    header = (out / "resolution_summary.tsv").read_text().splitlines()[0]
    assert "species_resolved" not in header
    assert (out / "resolution_groups.tsv").read_text().count("s1; s2") == 1


def test_the_caveat_is_always_stated(tmp_path):
    assert run_main(tmp_path, ["-p", f"v4={F},{R}", "-t", taxonomy(tmp_path)]) == 0
    log = (tmp_path / "out" / "resolution_log.txt").read_text()
    assert "property of the reference and the region, not of your samples" in log


def test_every_region_failing_is_fatal(tmp_path):
    assert run_main(tmp_path, ["-p", "bogus=TTTTTTTTTTTTTTTT,AAAAAAAAAAAAAAAA"]) == 1
    log = (tmp_path / "out" / "resolution_log.txt").read_text()
    assert "no region could be measured" in log


def test_one_failing_region_does_not_abort_the_comparison(tmp_path):
    """The whole point is comparing regions, so a reference that cannot answer for one
    of them must still report the others."""
    assert run_main(tmp_path, ["-p", f"v4={F},{R}",
                               "-p", "bogus=TTTTTTTTTTTTTTTT,AAAAAAAAAAAAAAAA",
                               "-t", taxonomy(tmp_path)]) == 0
    rows = {x["region"]: x for x in
            csv.DictReader((tmp_path / "out" / "resolution_summary.tsv").open(),
                           delimiter="\t")}
    # s1+s2 collapse to one sequence, s3 gives another, s4 yields no region
    assert rows["v4"]["distinct_sequences"] == "2"
    assert rows["bogus"]["distinct_sequences"] == "0"
    assert rows["bogus"]["no_region"] == "4"
    log = (tmp_path / "out" / "resolution_log.txt").read_text()
    assert "reported as zeroes, not as regions that resolve nothing" in log


def test_a_missing_forward_site_is_named_as_such(tmp_path):
    """A reference trimmed inside the amplicon has lost its forward site. Saying
    'no region found' would blame the region instead of the reference."""
    records = ar.read_fasta(reference(tmp_path))
    # keep the reverse site, remove the forward one
    trimmed = {n: s[s.index(R_RC):] for n, s in records.items() if R_RC in s}
    p = tmp_path / "trimmed.fasta"
    p.write_text("".join(f">{n}\n{s}\n" for n, s in trimmed.items()))
    assert r.main(["-m", str(p), "-o", str(tmp_path / "out2"),
                   "-p", f"v4={F},{R}"]) == 1
    log = (tmp_path / "out2" / "resolution_log.txt").read_text()
    assert "the forward site is absent everywhere" in log
    assert "trimmed to start inside the amplicon" in log


def test_primer_presence_counts_both_orientations(tmp_path):
    records = ar.read_fasta(reference(tmp_path))
    assert r.primer_presence(records, F, 0) == (4, 0)
    assert r.primer_presence(records, R, 0) == (0, 3)


def test_mismatched_taxonomy_ids_fire(tmp_path):
    p = tmp_path / "other.tsv"
    p.write_text("zz1\td__Bacteria\n")
    assert run_main(tmp_path, ["-p", f"v4={F},{R}", "-t", str(p)]) == 1
    assert "ids do not match" in (tmp_path / "out" / "resolution_log.txt").read_text()


def test_partial_taxonomy_warns_but_runs(tmp_path, caplog):
    p = tmp_path / "part.tsv"
    p.write_text("s1\td__Bacteria; s__Staphylococcus aureus\n")
    assert run_main(tmp_path, ["-p", f"v4={F},{R}", "-t", str(p)]) == 0
    assert "have no taxonomy row" in (tmp_path / "out" / "resolution_log.txt").read_text()


def test_no_region_asked_for_fires(tmp_path, caplog):
    """Argument checks fire before any output directory is made, so the message goes to
    the logger's fallback rather than to a file in a directory we should not create."""
    with caplog.at_level("ERROR", logger="resolution"):
        assert run_main(tmp_path, []) == 1
    assert "at least one --primers" in caplog.text
    assert not (tmp_path / "out").exists()


def test_missing_reference_fires(tmp_path):
    assert r.main(["-m", str(tmp_path / "nope.fasta"), "-o", str(tmp_path / "out"),
                   "--full-length"]) == 1


def test_negative_mismatch_budget_fires(tmp_path):
    assert run_main(tmp_path, ["--full-length", "--max-primer-mismatch", "-1"]) == 1
