"""Region extraction and the amplicon-length cross-check.

The shared module was covered only indirectly, through 04_mock.py and 06_resolution.py.
The cross-check added here exists to catch one specific mistake: primers for one region
configured alongside `quality.amplicon_len` for another. That mistake is invisible
downstream, because the overlap floor is amplicon_len + min overlap + margin, so a value
left too short passes trivially and DADA2 merges almost nothing after a full run.

Success criteria fixed before the tests were run: a matched pair passes, a region mix-up
of V4 against V3-V4 is refused, primers that match nothing are refused with a different
message, and natural within-region variation never fires the guard.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from amplicon_regions import (  # noqa: E402
    IUPAC, RegionError, check_amplicon_len, extract_region, length_distribution,
    read_fasta, region_lengths, revcomp,
)

FWD = "GTGCCAGCMGCCGCGGTAA"          # 515F, carries M
REV = "GGACTACHVGGGTWTCTAAT"         # 806R, carries H, V and W


def resolve(primer: str) -> str:
    """A concrete sequence the ambiguous primer matches: first base each code allows.

    A reference carries real bases, never IUPAC codes. Writing the codes into the
    reference is what a naive fixture does, and it makes every ambiguous position count
    as a mismatch.
    """
    return "".join(IUPAC[b][0] for b in primer)


def build(region: str, fwd: str = FWD, rev: str = REV,
          lead: str = "ACGTACGTAC", tail: str = "TGCATGCATG") -> str:
    """A reference sequence carrying both primer sites around a known region."""
    return lead + resolve(fwd) + region + revcomp(resolve(rev)) + tail


def fasta(tmp_path, seqs: dict[str, str]) -> str:
    path = tmp_path / "ref.fasta"
    path.write_text("".join(f">{n}\n{s}\n" for n, s in seqs.items()))
    return str(path)


def test_the_region_between_the_primers_is_recovered_exactly():
    region = "A" * 253
    assert extract_region(build(region), FWD, REV, 0) == region


def test_iupac_codes_in_the_primer_match_every_base_they_stand_for():
    # M is A or C. Both must match at zero mismatches, or a real primer whose site
    # carries the other base silently finds nothing.
    idx = FWD.index("M")
    for base in "AC":
        site = list(resolve(FWD))
        site[idx] = base
        seq = "TTTT" + "".join(site) + "GGGG" + revcomp(resolve(REV)) + "TTTT"
        assert extract_region(seq, FWD, REV, 0) == "GGGG"


def test_region_lengths_reports_what_was_found_and_what_was_not(tmp_path):
    good = build("A" * 253)
    trimmed = "A" * 100                       # no primer sites at all
    path = fasta(tmp_path, {"good": good, "also_good": build("A" * 255),
                            "trimmed": trimmed})
    records = read_fasta(path)
    lengths, missing = region_lengths(records, FWD, REV, 2)
    assert lengths == [253, 255]
    assert missing == ["trimmed"]


def test_a_matched_primer_pair_and_amplicon_len_pass():
    summary = check_amplicon_len([252, 253, 254], [], configured=253, tolerance=50)
    assert summary["median"] == 253
    assert summary["n_found"] == 3


def test_natural_within_region_variation_does_not_fire_the_guard():
    # V4 references vary by a few bases between taxa. That must never be an error.
    lengths = [248, 251, 253, 253, 254, 259]
    assert check_amplicon_len(lengths, [], configured=253, tolerance=50)["median"] == 253


def test_v4_length_against_v3v4_primers_is_refused():
    # The mistake the check exists for: primers cut ~465 bp, amplicon_len left at V4's 253.
    with pytest.raises(RegionError) as exc:
        check_amplicon_len([460, 465, 470], [], configured=253, tolerance=50)
    msg = str(exc.value)
    assert "different regions" in msg
    assert "465" in msg and "253" in msg


def test_the_guard_fires_symmetrically_when_amplicon_len_is_too_long():
    with pytest.raises(RegionError, match="different regions"):
        check_amplicon_len([253], [], configured=465, tolerance=50)


def test_primers_matching_nothing_are_refused_with_their_own_message():
    with pytest.raises(RegionError) as exc:
        check_amplicon_len([], ["a", "b"], configured=253, tolerance=50)
    msg = str(exc.value)
    assert "no reference yielded a region" in msg
    assert "different regions" not in msg      # not a length problem, a primer problem


def test_a_difference_of_exactly_the_primer_length_is_refused():
    """The real case: 341F/805R quoted as a 465 bp insert, region between them 427 bp.

    17 + 21 = 38 = 465 - 427, which sits inside a 50 bp tolerance and would otherwise pass.
    """
    with pytest.raises(RegionError, match="combined length of the two primers"):
        check_amplicon_len([427, 427, 428], [], 465, 50, primer_len=38)


def test_the_primer_length_check_does_not_fire_on_ordinary_variation():
    # Same 38 bp primers, but the difference is 10 bp, not 38.
    summary = check_amplicon_len([427, 427, 428], [], 437, 50, primer_len=38)
    assert summary["median"] == 427
    assert "primer_length_mismatch" not in summary


def test_the_primer_length_check_is_skipped_when_no_primer_length_is_given():
    # Older callers pass no primer_len; behaviour must be unchanged for them.
    summary = check_amplicon_len([427, 427, 428], [], 465, 50)
    assert summary["median"] == 427


def test_the_primer_length_check_only_fires_in_the_over_counting_direction():
    # A configured length SHORTER than the region by the primer length is a different
    # mistake and must not be reported as primers counted twice.
    summary = check_amplicon_len([427, 427, 428], [], 427 - 38, 50, primer_len=38)
    assert summary["median"] == 427


def test_the_boundary_of_the_tolerance_is_inclusive():
    # Exactly at tolerance passes; one base beyond does not.
    assert check_amplicon_len([303], [], configured=253, tolerance=50)["median"] == 303
    with pytest.raises(RegionError):
        check_amplicon_len([304], [], configured=253, tolerance=50)


def test_the_length_source_can_be_named_without_changing_the_old_wording():
    # 02_quality.py's message must not move; 03_dada2.py passes ASVs, not references.
    with pytest.raises(RegionError, match="the primers cut these references to a median"):
        check_amplicon_len([460, 465], [], configured=253, tolerance=50)
    with pytest.raises(RegionError, match="the denoised ASVs have a median"):
        check_amplicon_len([460, 465], [], configured=253, tolerance=50,
                           observed_as="the denoised ASVs have")


# ---- the tail the median cannot see --------------------------------------
#
# Criteria fixed before these were run: a catalogue sitting on the amplicon reports no
# tail, the PRJNA643648 shape (a correct median with 8.9% of ASVs at up to 1.8x) reports
# one, repeated lengths are reported as clusters, the cutoffs are inclusive, and an empty
# catalogue is refused rather than summarised as zero.

# The V4 arm's shape, scaled down: a median of exactly 253 bp and an oversize tail.
V4_WITH_TAIL = {f"in{i}": L for i, L in
                enumerate([253, 253, 253, 253, 252, 254, 248, 258, 251, 290], 1)}
V4_WITH_TAIL.update({"big1": 441, "big2": 441, "big3": 446, "big4": 456})


def test_a_catalogue_on_the_amplicon_has_no_tail():
    d = length_distribution({"a": 253, "b": 252, "c": 254}, 253, 50)
    assert d["n_oversize"] == 0 and d["n_undersize"] == 0
    assert d["oversize"] == [] and d["undersize"] == []
    assert d["median"] == 253 and d["clusters"] == []
    assert d["low_cutoff"] == 203 and d["high_cutoff"] == 303


def test_the_tail_is_found_although_the_median_is_exactly_right():
    """The whole point. check_amplicon_len passes this catalogue; this does not."""
    assert check_amplicon_len(sorted(V4_WITH_TAIL.values()), [], 253, 50)["median"] == 253
    d = length_distribution(V4_WITH_TAIL, 253, 50)
    assert d["median"] == 253                     # the median says nothing is wrong
    assert d["n_oversize"] == 4                   # and four ASVs are at up to 1.8x
    assert d["n_total"] == 14
    assert round(d["fraction_oversize"], 4) == round(4 / 14, 4)
    assert [L for _, L in d["oversize"]] == [456, 446, 441, 441]   # longest first
    assert d["n_undersize"] == 0


def test_repeated_off_length_values_are_reported_as_clusters():
    """17 at 441 bp is one amplified product; 17 scattered lengths are not."""
    d = length_distribution(V4_WITH_TAIL, 253, 50)
    assert d["clusters"] == [(441, 2)]            # only lengths seen more than once
    many = {f"a{i}": 441 for i in range(17)}
    many.update({f"b{i}": 446 for i in range(16)})
    many["in"] = 253
    d = length_distribution(many, 253, 50)
    assert d["clusters"] == [(441, 17), (446, 16)]  # commonest first


def test_both_cutoffs_are_inclusive():
    d = length_distribution({"hi": 303, "lo": 203, "over": 304, "under": 202}, 253, 50)
    assert [n for n, _ in d["oversize"]] == ["over"]
    assert [n for n, _ in d["undersize"]] == ["under"]


def test_a_zero_tolerance_is_allowed_and_reports_every_difference():
    d = length_distribution({"a": 253, "b": 254}, 253, 0)
    assert d["n_oversize"] == 1 and d["n_undersize"] == 0


def test_the_guard_fires_on_an_empty_catalogue_or_nonsense_settings():
    with pytest.raises(RegionError, match="empty catalogue"):
        length_distribution({}, 253, 50)
    with pytest.raises(RegionError, match="tolerance cannot be negative"):
        length_distribution({"a": 253}, 253, -1)
    with pytest.raises(RegionError, match="must be positive"):
        length_distribution({"a": 253}, 0, 50)
