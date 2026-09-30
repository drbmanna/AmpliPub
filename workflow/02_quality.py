#!/usr/bin/env python3
"""Propose DADA2 truncation lengths from read quality, and refuse lengths that cannot merge.

Give it the paired-end reads from 01_primers.py (trimmed.qza). The script

  1. runs qiime demux summarize, which draws a random subset of reads,
  2. reads the median quality at every position of R1 and R2,
  3. truncates each read just before the first position whose median falls below Q30,
  4. checks that the truncated reads still overlap enough to merge,
  5. writes the lengths for 03_dada2 and the quality profile behind them,
  6. logs the command, versions and results.

The rule and the overlap floor were fixed before the full data were seen. The floor is
amplicon length + DADA2 minimum overlap + a margin for length variation; for V4 that is
253 + 12 + 20 = 285 bp. If trunc-len-f + trunc-len-r falls below it, the reads would not
merge and DADA2 would return a nearly empty table without raising an error, so the
script stops instead.

qiime demux summarize has no seed, so a rerun can shift a median slightly. The exact
table used is kept in quality.qzv and quality_profile.tsv.

Standard library only. Needs Python 3.8 or later, Linux or WSL, and conda with the
QIIME 2 environment (see workflow/setup_envs.sh).

Example:
    python workflow/02_quality.py -i ~/research/baxter2016/q2/primers/trimmed.qza \\
        -o ~/research/baxter2016/q2/quality
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import platform
import shlex
import signal
import subprocess
import sys
import zipfile
from datetime import datetime, timezone

__version__ = "0.1.0"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from amplicon_regions import (  # noqa: E402
    RegionError, check_amplicon_len, read_fasta, region_lengths,
)

MIN_Q = 30          # median quality a position must reach to be kept
AMPLICON_LEN = 253  # V4, 515F to 806R, primers excluded
MIN_OVERLAP = 12    # qiime dada2 denoise-paired --p-min-overlap default in 2025.7
MARGIN = 20         # allowance for V4 length variation between taxa
N_SAMPLED = 10000   # qiime demux summarize --p-n default
LEN_TOLERANCE = 50  # bp between amplicon_len and what the primers cut from a reference.
#                     Chosen, not published: within one region references vary by tens of
#                     bases, between regions by hundreds, so this sits in the gap.
MIN_REACH = 0.75    # least fraction of sampled reads that must survive truncation.
#                     Adapted from nf-core/ampliseq's trunc_rmin default, but see
#                     check_reach: theirs picks a cutoff, ours checks one. A chosen value.
#                     Measured on the validated Baxter V4 run: 100.0% at both cutoffs, so
#                     this floor does not move that result.

log = logging.getLogger("quality")


class QualityError(RuntimeError):
    """A check failed. The message says which one and why."""


def run_cmd(cmd: list[str], timeout: int) -> tuple[int, str, str]:
    """Run a command with stdin closed and a hard time limit.

    The command gets its own process group, so a timeout also kills what it started.
    """
    log.info("run: %s", shlex.join(cmd))
    try:
        proc = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True, start_new_session=True)
    except FileNotFoundError as exc:
        raise QualityError(f"cannot start {cmd[0]}: {exc}") from exc
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.communicate()
        raise QualityError(f"timed out after {timeout} s and was killed: {shlex.join(cmd)}")
    return proc.returncode, out, err


def _tail(text: str, n: int = 10) -> str:
    return "\n".join(text.strip().splitlines()[-n:])


def parse_summary(text: str, name: str) -> dict[str, list[float]]:
    """Read a demux summarize seven-number table: one column per position, one row per stat."""
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        raise QualityError(f"{name} is empty")
    positions = lines[0].split("\t")[1:]
    if positions != [str(i) for i in range(1, len(positions) + 1)]:
        raise QualityError(f"{name}: positions are not 1..N in order")
    rows: dict[str, list[float]] = {}
    for ln in lines[1:]:
        label, *values = ln.split("\t")
        if len(values) != len(positions):
            raise QualityError(f"{name}: row {label!r} has {len(values)} values for "
                               f"{len(positions)} positions")
        rows[label] = [float(v) for v in values]
    for label in ("count", "50%"):
        if label not in rows:
            raise QualityError(f"{name}: no {label!r} row")
    return rows


def read_profiles(qzv: str) -> tuple[dict, dict]:
    try:
        with zipfile.ZipFile(qzv) as zf:
            out = []
            for direction in ("forward", "reverse"):
                suffix = f"/data/{direction}-seven-number-summaries.tsv"
                hits = [n for n in zf.namelist() if n.endswith(suffix)]
                if len(hits) != 1:
                    raise QualityError(f"{qzv}: expected one {direction} quality table, "
                                       f"found {len(hits)}. Is the input paired-end?")
                out.append(parse_summary(zf.read(hits[0]).decode("utf-8"), hits[0]))
    except zipfile.BadZipFile as exc:
        raise QualityError(f"{qzv} is not a readable visualization: {exc}") from exc
    return out[0], out[1]


def trunc_len(medians: list[float], min_q: float, name: str) -> int:
    """Keep positions up to, not including, the first whose median is below min_q."""
    for i, m in enumerate(medians):
        if m < min_q:
            if i == 0:
                raise QualityError(f"{name}: median quality is below Q{min_q:g} at position 1")
            return i
    return len(medians)


def check_overlap(f: int, r: int, amplicon: int, min_overlap: int, margin: int) -> int:
    """Return the expected overlap. Stop if the reads could not merge."""
    need = amplicon + min_overlap + margin
    if f + r < need:
        raise QualityError(
            f"trunc-len-f {f} + trunc-len-r {r} = {f + r} bp, below the {need} bp needed "
            f"(amplicon {amplicon} + min overlap {min_overlap} + margin {margin}). The reads "
            "would not merge. Do not lower the quality threshold to get past this without a "
            "written reason")
    return f + r - amplicon


def reach(counts: list[float], length: int) -> float:
    """Percent of sampled reads at least `length` long. DADA2 discards shorter reads."""
    return 100.0 * counts[length - 1] / counts[0] if counts[0] else 0.0


def retention_cap(counts: list[float], min_reach: float) -> int:
    """Highest position still reached by at least `min_reach` of the sampled reads.

    The quality profile runs to the longest read in the sample, not to the length most
    reads have. Where a minority of reads is longer, every statistic past that point is
    computed from that minority.
    """
    total = counts[0] if counts else 0.0
    if not total:
        return 0
    ok = [i + 1 for i, c in enumerate(counts) if c / total >= min_reach]
    return max(ok) if ok else 0


def apply_retention_cap(name: str, chosen: int, cap: int, counts: list[float],
                        min_reach: float) -> int:
    """Hold truncation to a position the reads actually reach, and say so when it binds.

    Found on real data (PRJNA643648 V4O, 2026-09-28): cutadapt was run keeping untrimmed
    pairs, so most reads lost a 19 bp primer and a minority kept full length. The median
    quality at position 281 was Q35, comfortably above the threshold, but it was the median
    of the 1,002 reads out of 10,000 that still existed there. Truncating at 281 would have
    silently discarded 90% of R1: no error, ordinary-looking quality, a table built from a
    tenth of the data.

    So this is not a choice between a quality guarantee and a retention one. A median over
    10% of the reads is not a statistic about the sample, and the fix is to stop reading the
    profile where it stops describing the data.
    """
    if cap and chosen > cap:
        log.warning(
            "%s: truncation capped at %d, down from the %d the quality rule chose. Only "
            "%.1f%% of reads reach %d, under the required %.0f%%; %.1f%% reach %d. Past "
            "position %d the profile describes a minority of long reads, not the sample.",
            name, cap, chosen, reach(counts, chosen), chosen, min_reach * 100,
            reach(counts, cap), cap, cap)
        return cap
    return chosen


def check_reach(profiles: list[tuple[str, list[float], int]],
                min_reach: float) -> dict[str, float]:
    """Require that truncation keeps at least `min_reach` of the sampled reads.

    Returns the surviving fraction per read, keyed by name.

    DADA2 discards every read shorter than the truncation position, and nothing downstream
    treats that as an error. A cutoff that keeps a third of the reads produces a run that
    completes, a table that looks ordinary, and diversity estimates drawn from a third of
    the data. That is the failure class worth guarding: plausible, not loud. So this is a
    floor, not a log line.

    This is the backstop, not the main mechanism: apply_retention_cap already holds the
    cutoff to a position the reads reach. This fires when even that is not enough, for
    instance when the sample is so short that no position retains min_reach.

    The threshold matches nf-core/ampliseq's `trunc_rmin` default of 0.75, and after seeing
    what it caught on real data their approach is right: an earlier version of this file
    argued that capping would trade a quality guarantee for a retention one. That was wrong.
    A median computed over 10% of the reads is not a statement about the sample, so there is
    no guarantee being given up. The number is still ours to justify, not theirs.

    `reach()` is a percentage and `min_reach` is a fraction; converting here keeps that
    difference in one place.
    """
    kept = {name: reach(counts, length) / 100.0 for name, counts, length in profiles}
    lengths = {name: length for name, _, length in profiles}
    below = {name: frac for name, frac in kept.items() if frac < min_reach}
    if below:
        detail = "; ".join(
            f"{name} keeps {frac:.1%} at position {lengths[name]}"
            for name, frac in sorted(below.items()))
        raise QualityError(
            f"truncation discards too much: {detail}, below the required {min_reach:.0%}. "
            "DADA2 drops every read shorter than the cutoff without reporting it, so the "
            "run would finish on a fraction of the data. Either the reads are shorter than "
            "the profile suggests, or the quality threshold is cutting too late. Lower "
            "--min-reach only with a written reason")
    return kept


def write_outputs(outdir: str, fwd: dict, rev: dict, f: int, r: int, overlap: int,
                  args, kept: dict[str, float] | None = None,
                  from_quality: tuple[int, int] | None = None,
                  caps: tuple[int, int] | None = None) -> None:
    with open(os.path.join(outdir, "quality_profile.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["position", "median_f", "median_r", "count_f", "count_r"])
        for i in range(max(len(fwd["50%"]), len(rev["50%"]))):
            get = lambda rows, k: rows[k][i] if i < len(rows[k]) else ""
            w.writerow([i + 1, get(fwd, "50%"), get(rev, "50%"),
                        get(fwd, "count"), get(rev, "count")])
    with open(os.path.join(outdir, "trunc_len.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        rows = [("trunc_len_f", f), ("trunc_len_r", r), ("min_q", args.min_q),
                ("amplicon_len", args.amplicon_len), ("min_overlap", args.min_overlap),
                ("margin", args.margin), ("expected_overlap", overlap),
                ("n_sampled", args.n), ("min_reach", args.min_reach)]
        # What the run actually retained, not only what it required. A reader checking
        # whether a table was built from most of the data should not have to rerun anything.
        if kept is not None:
            rows += [("reach_f", round(kept["R1"], 6)), ("reach_r", round(kept["R2"], 6))]
        # Both the position the quality rule proposed and the retention limit, so a reader
        # can see whether the cap bound without rerunning anything. trunc_len_f and
        # trunc_len_r above stay the values actually used, which is what 03_dada2.py reads.
        if from_quality is not None:
            rows += [("trunc_len_f_from_quality", from_quality[0]),
                     ("trunc_len_r_from_quality", from_quality[1])]
        if caps is not None:
            rows += [("retention_cap_f", caps[0]), ("retention_cap_r", caps[1])]
        for key, value in rows:
            w.writerow([key, value])


def setup_logging(outdir: str | None) -> None:
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(fmt)
    log.addHandler(console)
    if outdir:
        fh = logging.FileHandler(os.path.join(outdir, "quality_log.txt"), mode="a")
        fh.setFormatter(fmt)
        log.addHandler(fh)


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Propose DADA2 truncation lengths from read quality.")
    p.add_argument("-i", "--input", required=True, help="paired-end reads, e.g. trimmed.qza")
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("--min-q", type=float, default=MIN_Q,
                   help=f"median quality a position must reach (default {MIN_Q})")
    p.add_argument("--amplicon-len", type=int, default=AMPLICON_LEN,
                   help=f"amplicon length without primers (default {AMPLICON_LEN}, V4)")
    p.add_argument("--min-overlap", type=int, default=MIN_OVERLAP,
                   help=f"DADA2 minimum overlap (default {MIN_OVERLAP}, as in QIIME 2 2025.7)")
    p.add_argument("--margin", type=int, default=MARGIN,
                   help=f"extra overlap for amplicon length variation (default {MARGIN})")
    p.add_argument("--reference", help="reference FASTA to check --amplicon-len against, "
                                       "using the primers the reads were trimmed with")
    p.add_argument("--forward", help="forward primer, required with --reference")
    p.add_argument("--reverse", help="reverse primer, required with --reference")
    p.add_argument("--max-mismatch", type=int, default=2,
                   help="mismatches allowed per primer site in the reference (default 2)")
    p.add_argument("--len-tolerance", type=int, default=LEN_TOLERANCE,
                   help=f"bp that --amplicon-len may differ from the reference median "
                        f"(default {LEN_TOLERANCE}, a chosen value, not a published one)")
    p.add_argument("--min-reach", type=float, default=MIN_REACH,
                   help=f"least fraction of sampled reads that must survive truncation "
                        f"(default {MIN_REACH}, a chosen value, not a published one)")
    p.add_argument("--n", type=int, default=N_SAMPLED,
                   help=f"reads sampled for the quality profile (default {N_SAMPLED})")
    p.add_argument("--env", default="amplipub-qiime2-2025.7", help="conda env with QIIME 2")
    p.add_argument("--timeout", type=int, default=3600,
                   help="seconds before demux summarize is killed (default 3600)")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p.parse_args(argv)


def check_amplicon_against_reference(args) -> dict[str, object] | None:
    """Cross-check --amplicon-len against what the primers cut from a reference.

    Optional, because many runs have no reference FASTA to hand. Where there is none the
    configured length stays unverified, and this says so rather than implying it was
    checked: an unverified number that looks checked is worse than one known to be
    unchecked.
    """
    if not args.reference:
        log.info("amplicon_len %d bp is NOT verified against a reference: none given. "
                 "Pass --reference with --forward and --reverse to have it checked "
                 "(a wrong value passes the overlap floor and costs a full denoising run)",
                 args.amplicon_len)
        return None
    if not (args.forward and args.reverse):
        raise QualityError("--reference needs --forward and --reverse: the amplicon length "
                           "is only defined by a primer pair.")
    path = os.path.abspath(os.path.expanduser(args.reference))
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise QualityError(f"reference FASTA not found or empty: {path}")
    try:
        records = read_fasta(path)
        lengths, missing = region_lengths(records, args.forward, args.reverse,
                                          args.max_mismatch)
        summary = check_amplicon_len(lengths, missing, args.amplicon_len,
                                     args.len_tolerance,
                                     primer_len=len(args.forward) + len(args.reverse))
    except RegionError as exc:
        raise QualityError(str(exc)) from exc
    log.info("amplicon_len %d bp agrees with the reference: median %d bp "
             "(range %d-%d, n = %d, tolerance %d)",
             args.amplicon_len, summary["median"], summary["min"], summary["max"],
             summary["n_found"], summary["tolerance"])
    if summary["n_missing"]:
        log.info("%d of %d references yielded no region at %d mismatches, so the median is "
                 "taken over the rest; a database trimmed to the amplicon has lost its "
                 "primer sites and does this legitimately",
                 summary["n_missing"], summary["n_found"] + summary["n_missing"],
                 args.max_mismatch)
    return summary


def run(args) -> None:
    qza = os.path.abspath(os.path.expanduser(args.input))
    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    if not os.path.isfile(qza) or os.path.getsize(qza) == 0:
        raise QualityError(f"input artifact not found or empty: {qza}")
    os.makedirs(outdir, exist_ok=True)
    setup_logging(outdir)
    log.info("=" * 70)
    log.info("quality %s | Python %s | %s", __version__, platform.python_version(), platform.platform())
    log.info("command: %s", shlex.join([sys.executable] + sys.argv))
    log.info("started: %s", datetime.now(timezone.utc).isoformat(timespec="seconds"))
    log.info("rule: truncate before the first position with median quality < Q%g; "
             "floor %d + %d + %d = %d bp", args.min_q, args.amplicon_len, args.min_overlap,
             args.margin, args.amplicon_len + args.min_overlap + args.margin)

    # Before the expensive call, because this is the cheap check that decides whether
    # spending it is worth anything. amplicon_len and the primer pair are independent
    # settings describing one region; when they disagree the overlap floor is built from
    # the wrong number, passes, and DADA2 merges almost nothing.
    check_amplicon_against_reference(args)

    qzv = os.path.join(outdir, "quality.qzv")
    rc, _, err = run_cmd(["conda", "run", "-n", args.env, "qiime", "demux", "summarize",
                          "--i-data", qza, "--p-n", str(args.n), "--o-visualization", qzv],
                         args.timeout)
    if rc != 0:
        raise QualityError(f"qiime demux summarize exited with code {rc}:\n{_tail(err)}")
    if not os.path.isfile(qzv):
        raise QualityError(f"qiime demux summarize finished but wrote no {qzv}")

    fwd, rev = read_profiles(qzv)
    f_q = trunc_len(fwd["50%"], args.min_q, "R1")
    r_q = trunc_len(rev["50%"], args.min_q, "R2")
    log.info("R1: %d positions, median below Q%g first at %s -> trunc-len-f %d",
             len(fwd["50%"]), args.min_q, f_q + 1 if f_q < len(fwd["50%"]) else "none", f_q)
    log.info("R2: %d positions, median below Q%g first at %s -> trunc-len-r %d",
             len(rev["50%"]), args.min_q, r_q + 1 if r_q < len(rev["50%"]) else "none", r_q)
    # The quality rule can land past the point where most reads end, because the profile
    # runs to the longest read rather than the typical one. Hold it back to where the reads
    # actually are, before the overlap check, so the overlap is computed on real lengths.
    f_cap = retention_cap(fwd["count"], args.min_reach)
    r_cap = retention_cap(rev["count"], args.min_reach)
    f = apply_retention_cap("R1", f_q, f_cap, fwd["count"], args.min_reach)
    r = apply_retention_cap("R2", r_q, r_cap, rev["count"], args.min_reach)
    overlap = check_overlap(f, r, args.amplicon_len, args.min_overlap, args.margin)
    log.info("expected overlap %d bp (DADA2 needs %d; floor with margin %d)",
             overlap, args.min_overlap, args.min_overlap + args.margin)
    for name, rows, length in (("R1", fwd, f), ("R2", rev, r)):
        log.info("%s: %.1f%% of sampled reads reach position %d (DADA2 discards shorter reads)",
                 name, reach(rows["count"], length), length)
    # Enforced, not just reported: see check_reach. Both reads are measured before either
    # can fail, so the log shows the whole picture rather than stopping at the first one.
    kept = check_reach([("R1", fwd["count"], f), ("R2", rev["count"], r)], args.min_reach)
    log.info("retention floor %.0f%% met: R1 %.1f%%, R2 %.1f%%",
             args.min_reach * 100, kept["R1"] * 100, kept["R2"] * 100)
    write_outputs(outdir, fwd, rev, f, r, overlap, args, kept, (f_q, r_q), (f_cap, r_cap))
    log.info("done: --p-trunc-len-f %d --p-trunc-len-r %d  (%s)", f, r,
             os.path.join(outdir, "trunc_len.tsv"))


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        run(args)
    except QualityError as exc:
        if not log.handlers:
            setup_logging(None)
        log.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
