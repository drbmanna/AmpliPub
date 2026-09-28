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
    IUPAC, RegionError, check_amplicon_len, extract_region, read_fasta, region_lengths,
    revcomp,
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


def test_the_boundary_of_the_tolerance_is_inclusive():
    # Exactly at tolerance passes; one base beyond does not.
    assert check_amplicon_len([303], [], configured=253, tolerance=50)["median"] == 303
    with pytest.raises(RegionError):
        check_amplicon_len([304], [], configured=253, tolerance=50)
