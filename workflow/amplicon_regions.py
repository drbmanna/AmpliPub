#!/usr/bin/env python3
"""Find the region a primer pair amplifies, and group sequences that collapse in it.

Shared by 04_mock.py, which asks what a mock community's reference looks like over the
sequenced region, and 06_resolution.py, which asks the same question of a whole reference
database: what can this region actually tell apart, and what does it not.

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
                       tolerance: int) -> dict[str, object]:
    """Refuse a configured amplicon length the references do not support.

    `quality.amplicon_len` and the primer pair are two independent settings that have to
    describe the same region. Nothing downstream catches them disagreeing: the overlap
    floor is amplicon + min overlap + margin, so a length left too short passes the check
    trivially, DADA2 runs to completion, and almost nothing merges. The cost is a full
    denoising run, and the symptom appears only as a flagged read loss afterwards.

    The tolerance exists to catch a region mix-up (V4's 253 bp against V3-V4's ~465), not
    to police natural variation. Within one region references vary by tens of bases;
    between regions they differ by hundreds. No published source sets this number, so it
    is a chosen value stated plainly, not a standard, and it is a CLI flag.
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
    if abs(median - configured) > tolerance:
        raise RegionError(
            f"configured amplicon_len {configured} bp, but the primers cut these references "
            f"to a median of {median} bp (range {lengths[0]}-{lengths[-1]}, n = {len(lengths)}), "
            f"a difference of {abs(median - configured)} bp above the {tolerance} bp tolerance. "
            "The primer pair and amplicon_len describe different regions. Fix whichever is "
            "wrong before denoising: the overlap floor is built from amplicon_len, so a wrong "
            "value passes the overlap check and DADA2 then merges almost nothing.")
    return summary


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
