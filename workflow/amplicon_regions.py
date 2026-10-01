#!/usr/bin/env python3
"""Find the region a primer pair amplifies, and group sequences that collapse in it.

Shared by 04_mock.py, which asks what a mock community's reference looks like over the
sequenced region, and 06_resolution.py, which asks the same question of a whole reference
database: what can this region actually tell apart, and what does it not. 02_quality.py
and 03_dada2.py use the length checks here to hold `quality.amplicon_len`, the primer
pair, and the lengths actually observed to the same region.

Everything here is IUPAC-aware, because primers carry ambiguity codes and treating them
as literal bases silently finds nothing.

Standard library only.
"""

from __future__ import annotations

import statistics

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T",
         "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A",
              "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
              "B": "V", "V": "B", "D": "H", "H": "D", "N": "N"}


class RegionError(RuntimeError):
    """A check failed. The message says which one and why."""


def revcomp(seq: str) -> str:
    try:
        return "".join(COMPLEMENT[b] for b in reversed(seq.upper()))
    except KeyError as exc:
        raise RegionError(f"cannot complement base {exc.args[0]!r} in {seq!r}") from exc


def read_fasta(path: str) -> dict[str, str]:
    """Read a FASTA file. Duplicate names are an error, not a silent overwrite."""
    records: dict[str, str] = {}
    name = None
    parts: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(parts).upper()
                name = line[1:].split()[0]
                if name in records:
                    raise RegionError(f"{path}: {name} appears more than once")
                parts = []
            else:
                if name is None:
                    raise RegionError(f"{path}: sequence data before the first header")
                parts.append(line)
    if name is not None:
        records[name] = "".join(parts).upper()
    if not records:
        raise RegionError(f"{path}: no sequences")
    return records


def mismatches(primer: str, window: str) -> int:
    """Count positions where the window does not satisfy the primer's IUPAC code."""
    n = 0
    for p, b in zip(primer, window):
        allowed = IUPAC.get(p)
        if allowed is None:
            raise RegionError(f"{p!r} is not a IUPAC code in primer {primer!r}")
        if b not in allowed:
            n += 1
    return n


def find_primer(seq: str, primer: str, max_mismatch: int) -> tuple[int, int] | None:
    """Leftmost window matching the primer within max_mismatch. Returns (start, end).

    An exact match anywhere wins over an earlier near match, because a near match is a
    fallback for references whose primer site carries a real mismatch, not a licence to
    cut at the first thing that looks close.
    """
    width = len(primer)
    best = None
    for i in range(len(seq) - width + 1):
        m = mismatches(primer, seq[i:i + width])
        if m == 0:
            return i, i + width
        if m <= max_mismatch and best is None:
            best = (i, i + width)
    return best


def extract_region(seq: str, fwd: str, rev: str, max_mismatch: int) -> str | None:
    """The sequence between the forward primer and the reverse primer's complement.

    The reverse primer is searched for as its reverse complement, which is how it appears
    on the strand the reference is written on. Returns None if either site is absent.
    """
    f = find_primer(seq, fwd, max_mismatch)
    if f is None:
        return None
    r = find_primer(seq[f[1]:], revcomp(rev), max_mismatch)
    if r is None:
        return None
    return seq[f[1]:f[1] + r[0]]


def group_by_region(records: dict[str, str], fwd: str, rev: str,
                    max_mismatch: int) -> tuple[dict[str, list[str]], list[str]]:
    """Map each distinct region sequence to every record that produces it.

    The grouping is the point: two records sharing one group are indistinguishable over
    this region, so no method can separate them from an amplicon of it. Also returns the
    records where no region was found at this mismatch budget.
    """
    groups: dict[str, list[str]] = {}
    missing = []
    for name, seq in records.items():
        region = extract_region(seq, fwd, rev, max_mismatch)
        if region is None or not region:
            missing.append(name)
            continue
        groups.setdefault(region, []).append(name)
    if not groups:
        raise RegionError("no reference record yielded a region between the primers. "
                          "Check the primers and the reference orientation")
    return groups, missing


def region_lengths(records: dict[str, str], fwd: str, rev: str,
                   max_mismatch: int) -> tuple[list[int], list[str]]:
    """Length of the region each reference yields between the primers.

    Returns the sorted lengths, one per record that yielded a region, and the names of
    the records where a primer site was absent. A reference trimmed to start inside the
    amplicon has lost its primer site and lands in the second list; that is normal for
    some databases and is reported rather than treated as an error.
    """
    lengths: list[int] = []
    missing: list[str] = []
    for name, seq in records.items():
        region = extract_region(seq, fwd, rev, max_mismatch)
        if region is None or not region:
            missing.append(name)
            continue
        lengths.append(len(region))
    return sorted(lengths), missing


def check_amplicon_len(lengths: list[int], missing: list[str], configured: int,
                       tolerance: int, primer_len: int | None = None,
                       observed_as: str = "the primers cut these references to"
                       ) -> dict[str, object]:
    """Refuse a configured amplicon length the observed lengths do not support.

    `quality.amplicon_len` and the primer pair are two independent settings that have to
    describe the same region. Nothing downstream catches them disagreeing: the overlap
    floor is amplicon + min overlap + margin, so a length left too short passes the check
    trivially, DADA2 runs to completion, and almost nothing merges. The cost is a full
    denoising run, and the symptom appears only as a flagged read loss afterwards.

    The tolerance exists to catch a region mix-up (V4's 253 bp against V3-V4's ~465), not
    to police natural variation. Within one region references vary by tens of bases;
    between regions they differ by hundreds. No published source sets this number, so it
    is a chosen value stated plainly, not a standard, and it is a CLI flag.

    `lengths` do not have to come from a reference. 03_dada2.py passes the lengths of the
    denoised ASVs, which asks the same question of the data instead of a database, so
    `observed_as` names where the lengths came from and the default keeps 02_quality.py's
    wording unchanged. The median is all this function reads; the tail it cannot see is
    `length_distribution`'s job.
    """
    if not lengths:
        raise RegionError(
            f"no reference yielded a region between the primers ({len(missing)} checked). "
            "The primers do not match this reference, or it is in the other orientation. "
            "Nothing can be concluded about the amplicon length.")
    median = int(statistics.median(lengths))
    summary = {"n_found": len(lengths), "n_missing": len(missing), "median": median,
               "min": lengths[0], "max": lengths[-1], "configured": configured,
               "tolerance": tolerance}
    # A difference of exactly the primer length is not variation between taxa, it is the
    # primers counted twice. Published inserts are usually quoted WITH the primers (341F
    # and 805R are given as 465 bp, while the region between them is 427), and amplicon_len
    # here means the region WITHOUT them. Caught on PRJNA643648, where 465 - 427 = 38 was
    # exactly len(341F) + len(805R) and sat inside a 50 bp tolerance, so the run proceeded
    # on a length that was wrong in a specific, diagnosable way.
    if primer_len and configured - median == primer_len:
        summary["primer_length_mismatch"] = True
        raise RegionError(
            f"configured amplicon_len {configured} bp is exactly {primer_len} bp more than "
            f"the {median} bp the primers cut from these references, and {primer_len} bp is "
            "the combined length of the two primers. amplicon_len is the region WITHOUT the "
            "primers, but published insert sizes usually include them. Use "
            f"{median} bp. This is not length variation between taxa: the difference matches "
            "the primers exactly.")
    if abs(median - configured) > tolerance:
        raise RegionError(
            f"configured amplicon_len {configured} bp, but {observed_as} "
            f"a median of {median} bp (range {lengths[0]}-{lengths[-1]}, n = {len(lengths)}), "
            f"a difference of {abs(median - configured)} bp above the {tolerance} bp tolerance. "
            "The primer pair and amplicon_len describe different regions. Fix whichever is "
            "wrong before denoising: the overlap floor is built from amplicon_len, so a wrong "
            "value passes the overlap check and DADA2 then merges almost nothing.")
    return summary


def length_distribution(lengths_by_id: dict[str, int], configured: int,
                        tolerance: int) -> dict[str, object]:
    """Group sequences by how far their length sits from the configured amplicon.

    Separate from `check_amplicon_len` on purpose, and both are needed. That function
    reads the median and asks whether the catalogue as a whole describes the configured
    region. **The median cannot see the tail.** On PRJNA643648's V4 arm the median was
    253 bp, exactly the configured value, while 72 of 806 ASVs (8.9%) sat at 294 to
    456 bp, up to 1.8x the amplicon, and every check in the pipeline passed in silence.
    This function is what would have said so.

    `tolerance` is bp either side of `configured`. Outside it a sequence is reported as
    oversize or undersize. No published source sets the point at which an ASV is too long
    for its amplicon: six sources were checked for the denoising criteria in 03_dada2.py
    and none of them sets a length threshold either, so the caller passes the value in as
    a chosen number and 03_dada2.py exposes it as a CLI flag. It is not a standard.

    Nothing here fails. Being off-length is evidence about a sequence, not a verdict on
    it: chimeras, off-target amplification and carryover between runs all produce
    off-length ASVs, and so does real length variation in a few taxa. Which one it is
    takes work this function does not do, so it reports and leaves the judgement.

    Returns the counts, shares and cutoffs, the off-length ids with their lengths, and the
    lengths that repeat across more than one ASV, commonest first. A cluster of identical
    off-length sequences is the shape that distinguishes one amplified product from
    scattered noise: that V4 arm had 17 ASVs at 441 bp and 16 at 446 bp.
    """
    if tolerance < 0:
        raise RegionError(f"tolerance cannot be negative, got {tolerance}")
    if configured <= 0:
        raise RegionError(f"configured amplicon length must be positive, got {configured}")
    if not lengths_by_id:
        raise RegionError("no sequences to measure. Nothing can be said about the length "
                          "distribution of an empty catalogue")
    low, high = configured - tolerance, configured + tolerance
    oversize = sorted(((n, L) for n, L in lengths_by_id.items() if L > high),
                      key=lambda x: (-x[1], x[0]))
    undersize = sorted(((n, L) for n, L in lengths_by_id.items() if L < low),
                       key=lambda x: (x[1], x[0]))
    lengths = sorted(lengths_by_id.values())
    counts: dict[int, int] = {}
    for _, L in oversize + undersize:
        counts[L] = counts.get(L, 0) + 1
    total = len(lengths_by_id)
    return {
        "n_total": total,
        "median": int(statistics.median(lengths)),
        "min": lengths[0],
        "max": lengths[-1],
        "configured": configured,
        "tolerance": tolerance,
        "low_cutoff": low,
        "high_cutoff": high,
        "oversize": oversize,
        "undersize": undersize,
        "n_oversize": len(oversize),
        "n_undersize": len(undersize),
        "fraction_oversize": len(oversize) / total,
        "fraction_undersize": len(undersize) / total,
        "clusters": sorted(((L, n) for L, n in counts.items() if n > 1),
                           key=lambda x: (-x[1], -x[0])),
    }


def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        raise RegionError(f"hamming needs equal lengths, got {len(a)} and {len(b)}")
    return sum(1 for x, y in zip(a, b) if x != y)


def nearest_same_length(seq: str, targets: list[str]) -> tuple[int, int] | None:
    """Fewest differing positions against any target of the same length, and its length."""
    same = [t for t in targets if len(t) == len(seq)]
    if not same:
        return None
    return min(hamming(seq, t) for t in same), len(seq)
