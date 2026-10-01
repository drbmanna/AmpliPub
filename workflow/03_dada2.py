#!/usr/bin/env python3
"""Denoise paired reads with DADA2, then judge the run against pre-registered criteria.

Give it the trimmed reads from 01_primers.py and the lengths from 02_quality.py. The
script

  1. reads trunc_len.tsv and repeats the overlap check before spending any compute,
  2. runs qiime dada2 denoise-paired with every parameter stated explicitly,
  3. reads the denoising stats back out of the artifact,
  4. measures the length of every ASV against `quality.amplicon_len`,
  5. applies the criteria below and writes them next to the results.

Why the criteria look like this. Six sources were checked on 2026-09-12: the DADA2
tutorial, the mothur MiSeq SOP, Kozich et al. 2013 (AEM 79:5112), Callahan et al. 2016
(Nat Methods 13:581), Estaki et al. 2020 (Curr Protoc Bioinformatics 70:e100), and the
QIIME 2 denoising tutorial. **None of them sets a pass or fail threshold.** They report
what they observed and leave the judgement to the reader. So this script hard-fails only
where the output is wrong rather than merely poor, and everything else is flagged, with
the source quoted in the message.

Hard fail (structural, the table cannot be believed):
  - truncation lengths that leave too little overlap to merge (the 02_quality check,
    repeated here because this is the boundary where the compute is spent),
  - a library that comes out with zero reads,
  - an empty feature table,
  - an ASV catalogue whose median length does not agree with amplicon_len. 02_quality.py
    runs the same check against a reference FASTA, but that is optional and many runs have
    none, so this repeats it against the data, where it cannot be skipped. A median far
    from amplicon_len means the configured region is not the region that was sequenced,
    and the overlap floor, the resolution claims and the region-matched classifier are all
    built on it. This fires after the denoise, so the compute is already spent and
    Snakemake will remove the artifacts: that is the same cost the zero-read and
    empty-table checks already carry at this boundary.

Flag and report (sourced, qualitative in the original):
  - more than half the reads lost at any single step outside filtering. DADA2 tutorial:
    "Outside of filtering, there should no step in which a majority of reads are lost."
    (sic, benjjneb.github.io/dada2/tutorial.html, fetched 2026-09-12),
  - more than half the merged reads lost to chimeras. Same source: "If most of your
    reads were removed as chimeric, upstream processing may need to be revisited."

Flag and report (no source exists at all, our own check):
  - ASVs further than --asv-len-tolerance from amplicon_len, in either direction. None of
    the six sources sets a length threshold for a denoised ASV, so the tolerance is a
    chosen value exposed as a flag, in the same way 02_quality.py exposes --len-tolerance
    and --min-reach. Written because nothing here caught it: on PRJNA643648's V4 arm the
    median ASV was 253 bp, exactly the configured amplicon, while 72 of 806 ASVs (8.9%)
    ran from 294 to 456 bp with 17 at 441 bp and 16 at 446 bp, and every existing check
    passed in silence. Being off-length is evidence, not a verdict. Chimeras, off-target
    amplification, carryover between runs on the same BioProject and genuine variation in
    a few taxa all produce it, and separating them is not this stage's job, so the ASVs
    are named in asv_length_flags.tsv and the judgement is left to whoever reads it.

No numeric threshold here is presented as a standard, because none exists. "More than
half" is the tutorial's own "majority" and "most", made countable.

One semantic trap, read out of run_dada.R in this build rather than assumed: in the
paired track table `denoised` is denoisedF, the **forward** reads only. The drop from
denoised to merged therefore carries both reverse denoising and merging, and the
per-step loss is reported with that caveat rather than as a clean merge rate.

Standard library only. Needs Python 3.8 or later, Linux or WSL, and conda with the
QIIME 2 environment (see workflow/setup_envs.sh).

Example:
    python workflow/03_dada2.py -i ~/research/baxter2016/q2/primers/trimmed.qza \\
        -t ~/research/baxter2016/q2/quality/trunc_len.tsv \\
        -o ~/research/baxter2016/q2/dada2 --threads 12
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

__version__ = "0.2.0"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from amplicon_regions import (  # noqa: E402
    RegionError, check_amplicon_len, length_distribution, read_fasta,
)

MAJORITY = 0.5  # the tutorial's "majority" and "most", made countable
ASV_LEN_TOLERANCE = 50  # bp an ASV may differ from amplicon_len before it is reported.
#                         The same chosen value and the same reasoning as 02_quality.py's
#                         LEN_TOLERANCE: within one region lengths vary by tens of bases,
#                         between regions by hundreds, so this sits in the gap. No
#                         published source sets it. On a 253 bp V4 amplicon it reports
#                         anything outside 203-303 bp, which leaves the 268-293 bp primer
#                         readthrough band 01_primers.py already handles alone and reports
#                         the group above it.

TUTORIAL = "benjjneb.github.io/dada2/tutorial.html, fetched 2026-09-12"
QUOTE_STEP = ('DADA2 tutorial: "Outside of filtering, there should no step in which a '
              'majority of reads are lost." (sic)')
QUOTE_CHIMERA = ('DADA2 tutorial: "If most of your reads were removed as chimeric, '
                 'upstream processing may need to be revisited."')

# Paired-end stats columns, read from q2_dada2/_denoise.py in the 2025.7 build.
COUNT_COLS = ["input", "filtered", "denoised", "merged", "non-chimeric"]

# Steps the tutorial rule applies to: everything after filtering.
STEPS = [("filtered", "denoised", "forward reads denoised"),
         ("denoised", "merged", "reverse denoising and merging together"),
         ("merged", "non-chimeric", "chimera removal")]

log = logging.getLogger("dada2")


class Dada2Error(RuntimeError):
    """A hard check failed. The message says which one and why."""


def run_cmd(cmd: list[str], timeout: int) -> tuple[int, str, str]:
    """Run a command with stdin closed and a hard time limit.

    The command gets its own process group, so a timeout also kills what it started.
    """
    log.info("run: %s", shlex.join(cmd))
    try:
        proc = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True, start_new_session=True)
    except FileNotFoundError as exc:
        raise Dada2Error(f"cannot start {cmd[0]}: {exc}") from exc
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.communicate()
        raise Dada2Error(f"timed out after {timeout} s and was killed: {shlex.join(cmd)}")
    return proc.returncode, out, err


def _tail(text: str, n: int = 15) -> str:
    return "\n".join(text.strip().splitlines()[-n:])


def read_trunc_len(path: str) -> dict[str, int]:
    """Read the two-column file 02_quality.py wrote. Every key it needs must be there."""
    values: dict[str, int] = {}
    with open(path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) != 2:
                raise Dada2Error(f"{path}: expected two columns, got {len(row)}: {row}")
            try:
                values[row[0]] = int(float(row[1]))
            except ValueError as exc:
                raise Dada2Error(f"{path}: {row[0]} is not a number: {row[1]!r}") from exc
    needed = ["trunc_len_f", "trunc_len_r", "amplicon_len", "min_overlap", "margin"]
    missing = [k for k in needed if k not in values]
    if missing:
        raise Dada2Error(f"{path} has no {', '.join(missing)}. Was it written by 02_quality.py?")
    for key in ("trunc_len_f", "trunc_len_r"):
        if values[key] <= 0:
            raise Dada2Error(f"{path}: {key} is {values[key]}")
    return values


def check_overlap(v: dict[str, int]) -> int:
    """Repeat the 02_quality check here, where the compute is about to be spent."""
    f, r = v["trunc_len_f"], v["trunc_len_r"]
    need = v["amplicon_len"] + v["min_overlap"] + v["margin"]
    if f + r < need:
        raise Dada2Error(
            f"trunc-len-f {f} + trunc-len-r {r} = {f + r} bp, below the {need} bp needed "
            f"(amplicon {v['amplicon_len']} + min overlap {v['min_overlap']} + margin "
            f"{v['margin']}). The reads would not merge and DADA2 would return a nearly "
            "empty table without raising an error")
    return f + r - v["amplicon_len"]


def read_stats(qza: str) -> list[dict[str, object]]:
    """Pull the per-sample track table out of a DADA2 stats artifact."""
    try:
        with zipfile.ZipFile(qza) as zf:
            hits = [n for n in zf.namelist() if n.endswith("/data/stats.tsv")]
            if len(hits) != 1:
                raise Dada2Error(f"{qza}: expected one stats table, found {len(hits)}")
            text = zf.read(hits[0]).decode("utf-8")
    except zipfile.BadZipFile as exc:
        raise Dada2Error(f"{qza} is not a readable artifact: {exc}") from exc
    return parse_stats(text, qza)


def parse_stats(text: str, name: str) -> list[dict[str, object]]:
    """Parse stats.tsv. The '#q2:types' directive row is metadata, not a sample."""
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    if not rows:
        raise Dada2Error(f"{name}: no rows")
    first = list(rows[0].keys())
    if first[0] != "sample-id":
        raise Dada2Error(f"{name}: first column is {first[0]!r}, expected 'sample-id'")
    missing = [c for c in COUNT_COLS if c not in first]
    if missing:
        raise Dada2Error(
            f"{name}: no {', '.join(missing)} column. Single-end runs have no 'merged' "
            "column; this script is for paired-end data")
    out = []
    for row in rows:
        if str(row["sample-id"]).startswith("#q2:types"):
            continue
        rec: dict[str, object] = {"sample-id": row["sample-id"]}
        for col in COUNT_COLS:
            try:
                rec[col] = int(float(row[col]))
            except (TypeError, ValueError) as exc:
                raise Dada2Error(f"{name}: {row['sample-id']} has a non-numeric "
                                 f"{col}: {row[col]!r}") from exc
        for col in COUNT_COLS:
            if rec[col] < 0:
                raise Dada2Error(f"{name}: {row['sample-id']} has a negative {col}")
        out.append(rec)
    if not out:
        raise Dada2Error(f"{name}: the table holds no samples")
    return out


def lost_fraction(before: int, after: int) -> float:
    """Fraction of reads lost between two steps. Zero in means nothing was lost."""
    return 0.0 if before <= 0 else 1.0 - (after / before)


def judge(stats: list[dict[str, object]], majority: float) -> tuple[list[dict], dict]:
    """Apply the criteria. Returns the per-sample flags and the pooled totals."""
    totals = {c: sum(int(s[c]) for s in stats) for c in COUNT_COLS}
    flags = []
    for s in stats:
        for before, after, note in STEPS:
            frac = lost_fraction(int(s[before]), int(s[after]))
            if frac > majority:
                flags.append({"sample-id": s["sample-id"], "step": f"{before}->{after}",
                              "lost_fraction": round(frac, 4),
                              "reads_before": s[before], "reads_after": s[after],
                              "note": note})
    return flags, totals


def zero_read_samples(stats: list[dict[str, object]]) -> list[str]:
    return [str(s["sample-id"]) for s in stats if int(s["non-chimeric"]) == 0]


def write_stats(path: str, stats: list[dict[str, object]]) -> None:
    """Per-sample counts plus the loss at each step, so nobody has to recompute them."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample-id"] + COUNT_COLS
                   + [f"lost_{b}_to_{a}" for b, a, _ in STEPS] + ["retained_of_input"])
        for s in stats:
            losses = [round(lost_fraction(int(s[b]), int(s[a])), 4) for b, a, _ in STEPS]
            kept = 0.0 if int(s["input"]) == 0 else int(s["non-chimeric"]) / int(s["input"])
            w.writerow([s["sample-id"]] + [s[c] for c in COUNT_COLS] + losses
                       + [round(kept, 4)])


def write_flags(path: str, flags: list[dict]) -> None:
    cols = ["sample-id", "step", "lost_fraction", "reads_before", "reads_after", "note"]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for f in flags:
            w.writerow([f[c] for c in cols])


# ---- ASV lengths ---------------------------------------------------------

def export_rep_seqs(qza: str, dest: str, env: str, timeout: int) -> str:
    """Export rep_seqs.qza and return the FASTA inside it."""
    rc, _, err = run_cmd(["conda", "run", "-n", env, "qiime", "tools", "export",
                          "--input-path", qza, "--output-path", dest], timeout)
    if rc != 0:
        raise Dada2Error(f"qiime tools export failed on {qza} with code {rc}:\n{_tail(err)}")
    path = os.path.join(dest, "dna-sequences.fasta")
    if not os.path.isfile(path):
        raise Dada2Error(f"export finished but wrote no {path}")
    return path


def asv_lengths(fasta: str) -> dict[str, int]:
    """Length of every denoised ASV, keyed by feature id."""
    try:
        records = read_fasta(fasta)
    except RegionError as exc:
        raise Dada2Error(f"cannot read the exported ASVs: {exc}") from exc
    return {name: len(seq) for name, seq in records.items()}


def check_asv_lengths(lengths: dict[str, int], configured: int,
                      tolerance: int) -> dict[str, object]:
    """The median must agree with amplicon_len; the tails are measured and returned.

    Two questions, and the second is the one that was being missed. The median says
    whether the catalogue describes the configured region at all, and disagreeing with it
    is structural, so it fails. The tails say what else is in the catalogue, and the
    median is blind to them: a run can sit exactly on 253 bp and still carry 8.9% of its
    ASVs at up to 1.8x that. Those are reported, not failed.
    """
    try:
        check_amplicon_len(sorted(lengths.values()), [], configured, tolerance,
                           observed_as="the denoised ASVs have")
        return length_distribution(lengths, configured, tolerance)
    except RegionError as exc:
        raise Dada2Error(
            f"the {len(lengths)} denoised ASVs do not describe the configured region. "
            f"{exc} Here the reads are the evidence, not a reference database, so this is "
            "amplicon_len disagreeing with what was actually sequenced") from exc


def write_asv_lengths(path: str, lengths: dict[str, int], dist: dict[str, object]) -> None:
    """Every ASV with its length and how far it sits from amplicon_len."""
    low, high = int(dist["low_cutoff"]), int(dist["high_cutoff"])
    configured = int(dist["configured"])
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["feature-id", "length", "difference_from_amplicon_len", "class"])
        for name, L in sorted(lengths.items(), key=lambda x: (-x[1], x[0])):
            kind = "oversize" if L > high else "undersize" if L < low else "in_range"
            w.writerow([name, L, L - configured, kind])


def write_length_flags(path: str, dist: dict[str, object]) -> None:
    """Only the ASVs outside the tolerance, longest first, then shortest first."""
    with open(path, "w", newline="") as fh:
        fh.write(f"# ASVs further than {dist['tolerance']} bp from amplicon_len "
                 f"{dist['configured']} bp, so outside {dist['low_cutoff']}-"
                 f"{dist['high_cutoff']} bp.\n")
        fh.write("# The tolerance is a chosen value, not a published threshold.\n"
                 "# Being off-length is evidence about an ASV, not a verdict on it:\n"
                 "# chimeras, off-target amplification and carryover between runs all\n"
                 "# look like this, and so does real variation in a few taxa.\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["feature-id", "length", "difference_from_amplicon_len", "class"])
        for name, L in dist["oversize"]:
            w.writerow([name, L, L - int(dist["configured"]), "oversize"])
        for name, L in dist["undersize"]:
            w.writerow([name, L, L - int(dist["configured"]), "undersize"])


def report_lengths(dist: dict[str, object]) -> None:
    """Say what the length distribution is, and shout only if something is outside it."""
    log.info("%d ASVs, median %d bp against amplicon_len %d (range %d-%d, tolerance "
             "+/-%d bp, a chosen value)", dist["n_total"], dist["median"],
             dist["configured"], dist["min"], dist["max"], dist["tolerance"])
    if not dist["n_oversize"] and not dist["n_undersize"]:
        log.info("every ASV is within %d bp of amplicon_len", dist["tolerance"])
        return
    for kind, cut, cmp_ in (("oversize", dist["high_cutoff"], "above"),
                            ("undersize", dist["low_cutoff"], "below")):
        n = int(dist[f"n_{kind}"])
        if not n:
            continue
        rows = dist[kind]
        log.warning("%d of %d ASVs (%.1f%%) are %s, %s %d bp: %d-%d bp. See "
                    "asv_length_flags.tsv", n, dist["n_total"],
                    100 * float(dist[f"fraction_{kind}"]), kind, cmp_, cut,
                    min(L for _, L in rows), max(L for _, L in rows))
    if dist["clusters"]:
        top = ", ".join(f"{n} at {L} bp" for L, n in dist["clusters"][:5])
        log.warning("off-length ASVs cluster at identical lengths (%s), which is the "
                    "shape of one amplified product rather than scattered noise. "
                    "Chimeras, off-target amplification and carryover between runs are "
                    "all candidates; nothing here tells them apart", top)
    log.warning("no published source sets a length threshold for a denoised ASV, so the "
                "%d bp tolerance is our own choice and is recorded in criteria.tsv. "
                "These ASVs are flagged, not removed", dist["tolerance"])


def write_criteria(path: str, args, v: dict[str, int], overlap: int) -> None:
    """Write down what was being tested for, in the same directory as the result."""
    with open(path, "w", newline="") as fh:
        fh.write("# Criteria applied by 03_dada2.py. Pre-registered before the run.\n")
        fh.write("# No published source sets a pass/fail threshold; six were checked on\n")
        fh.write("# 2026-09-12 and all report observed numbers only. The two flag rules\n")
        fh.write(f"# below are the DADA2 tutorial's own wording ({TUTORIAL}).\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["criterion", "action", "value", "source"])
        w.writerow(["overlap after truncation", "hard fail",
                    f"trunc_f+trunc_r >= {v['amplicon_len'] + v['min_overlap'] + v['margin']} bp",
                    "arithmetic; min overlap is the QIIME 2 2025.7 default"])
        w.writerow(["library with zero reads", "hard fail" if not args.allow_zero_read_samples
                    else "flag (override passed)", "non-chimeric > 0", "our choice"])
        w.writerow(["empty feature table", "hard fail", "total non-chimeric > 0", "our choice"])
        w.writerow(["reads lost at one step outside filtering", "flag",
                    f"<= {args.majority:.0%}", QUOTE_STEP])
        w.writerow(["reads lost to chimeras", "flag", f"<= {args.majority:.0%}", QUOTE_CHIMERA])
        w.writerow(["median ASV length against amplicon_len", "hard fail",
                    f"|median - {v['amplicon_len']}| <= {args.asv_len_tolerance} bp",
                    "arithmetic; the tolerance is our choice, no source sets one"])
        w.writerow(["ASVs outside amplicon_len +/- tolerance", "flag",
                    f"{v['amplicon_len'] - args.asv_len_tolerance}-"
                    f"{v['amplicon_len'] + args.asv_len_tolerance} bp",
                    "our choice; none of the six sources sets a length threshold"])
        w.writerow(["asv_len_tolerance", "parameter", args.asv_len_tolerance,
                    "a chosen value, not a published one"])
        w.writerow(["trunc_len_f", "parameter", v["trunc_len_f"], "02_quality.py"])
        w.writerow(["trunc_len_r", "parameter", v["trunc_len_r"], "02_quality.py"])
        w.writerow(["expected_overlap", "derived", overlap, "02_quality.py"])
        w.writerow(["n_threads", "parameter", args.threads, "set explicitly; QIIME 2 default is 1"])


def setup_logging(outdir: str | None) -> None:
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(fmt)
    log.addHandler(console)
    if outdir:
        fh = logging.FileHandler(os.path.join(outdir, "dada2_log.txt"), mode="a")
        fh.setFormatter(fmt)
        log.addHandler(fh)


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Denoise paired reads with DADA2 and judge the run against "
                    "pre-registered criteria.")
    p.add_argument("-i", "--input", required=True, help="paired-end reads, e.g. trimmed.qza")
    p.add_argument("-t", "--trunc-len", required=True,
                   help="trunc_len.tsv written by 02_quality.py")
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("--threads", type=int, default=1,
                   help="DADA2 threads. QIIME 2's own default is 1, so this is always "
                        "passed explicitly (default 1)")
    p.add_argument("--majority", type=float, default=MAJORITY,
                   help=f"fraction that counts as the tutorial's 'majority' (default {MAJORITY})")
    p.add_argument("--asv-len-tolerance", type=int, default=ASV_LEN_TOLERANCE,
                   help=f"bp an ASV may differ from amplicon_len before it is reported, "
                        f"and bp the median may differ before the run fails "
                        f"(default {ASV_LEN_TOLERANCE}, a chosen value, not a published one)")
    p.add_argument("--allow-zero-read-samples", action="store_true",
                   help="downgrade the zero-read library check from hard fail to a flag. "
                        "A deliberate choice, recorded in criteria.tsv")
    p.add_argument("--env", default="amplipub-qiime2-2025.7", help="conda env with QIIME 2")
    p.add_argument("--timeout", type=int, default=86400,
                   help="seconds before denoise-paired is killed (default 86400)")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p.parse_args(argv)


def run(args) -> None:
    qza = os.path.abspath(os.path.expanduser(args.input))
    trunc = os.path.abspath(os.path.expanduser(args.trunc_len))
    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    if not os.path.isfile(qza) or os.path.getsize(qza) == 0:
        raise Dada2Error(f"input artifact not found or empty: {qza}")
    if not os.path.isfile(trunc) or os.path.getsize(trunc) == 0:
        raise Dada2Error(f"trunc_len.tsv not found or empty: {trunc}")
    if args.threads < 1:
        raise Dada2Error(f"--threads must be at least 1, got {args.threads}")
    if not 0 < args.majority < 1:
        raise Dada2Error(f"--majority must be between 0 and 1, got {args.majority}")
    if args.asv_len_tolerance < 0:
        raise Dada2Error("--asv-len-tolerance cannot be negative, got "
                         f"{args.asv_len_tolerance}")
    os.makedirs(outdir, exist_ok=True)
    setup_logging(outdir)
    log.info("=" * 70)
    log.info("dada2 %s | Python %s | %s", __version__, platform.python_version(),
             platform.platform())
    log.info("command: %s", shlex.join([sys.executable] + sys.argv))
    log.info("started: %s", datetime.now(timezone.utc).isoformat(timespec="seconds"))

    v = read_trunc_len(trunc)
    overlap = check_overlap(v)
    log.info("trunc-len-f %d, trunc-len-r %d, expected overlap %d bp",
             v["trunc_len_f"], v["trunc_len_r"], overlap)
    log.info("criteria: hard fail on no overlap, a zero-read library%s, an empty table, or "
             "a median ASV length more than %d bp from amplicon_len %d; flag above %.0f%% "
             "loss at any step after filtering and any ASV outside that tolerance. No "
             "published source sets a threshold, see criteria.tsv", " (overridden to a flag)"
             if args.allow_zero_read_samples else "", args.asv_len_tolerance,
             v["amplicon_len"], args.majority * 100)

    table = os.path.join(outdir, "table.qza")
    rep_seqs = os.path.join(outdir, "rep_seqs.qza")
    stats_qza = os.path.join(outdir, "denoising_stats.qza")
    rc, _, err = run_cmd(
        ["conda", "run", "-n", args.env, "qiime", "dada2", "denoise-paired",
         "--i-demultiplexed-seqs", qza,
         "--p-trunc-len-f", str(v["trunc_len_f"]),
         "--p-trunc-len-r", str(v["trunc_len_r"]),
         "--p-min-overlap", str(v["min_overlap"]),
         "--p-n-threads", str(args.threads),
         "--o-table", table, "--o-representative-sequences", rep_seqs,
         "--o-denoising-stats", stats_qza],
        args.timeout)
    if rc != 0:
        raise Dada2Error(f"qiime dada2 denoise-paired exited with code {rc}:\n{_tail(err)}")
    for path in (table, rep_seqs, stats_qza):
        if not os.path.isfile(path):
            raise Dada2Error(f"denoise-paired finished but wrote no {path}")

    stats = read_stats(stats_qza)
    flags, totals = judge(stats, args.majority)
    write_stats(os.path.join(outdir, "dada2_stats.tsv"), stats)
    write_flags(os.path.join(outdir, "dada2_flags.tsv"), flags)
    write_criteria(os.path.join(outdir, "criteria.tsv"), args, v, overlap)

    kept = 0.0 if totals["input"] == 0 else 100.0 * totals["non-chimeric"] / totals["input"]
    log.info("%d samples, %d reads in, %d non-chimeric (%.1f%% retained overall)",
             len(stats), totals["input"], totals["non-chimeric"], kept)
    for before, after, note in STEPS:
        log.info("  %s -> %s: %.1f%% lost pooled (%s)", before, after,
                 100 * lost_fraction(totals[before], totals[after]), note)
    log.info("note: in the paired track table 'denoised' is denoisedF, the forward reads "
             "only, so denoised -> merged carries reverse denoising as well as merging")

    empty = zero_read_samples(stats)
    if totals["non-chimeric"] == 0:
        raise Dada2Error("the feature table is empty: no sample kept a single read. "
                         "Check the truncation lengths and the primer stage before rerunning")
    if empty and not args.allow_zero_read_samples:
        raise Dada2Error(
            f"{len(empty)} librar{'y' if len(empty) == 1 else 'ies'} came out with zero "
            f"reads: {', '.join(empty[:10])}{' ...' if len(empty) > 10 else ''}. Decide "
            "whether to drop them and say so in the methods, then rerun with "
            "--allow-zero-read-samples")
    if empty:
        log.warning("%d librar%s came out with zero reads and the hard check was "
                    "overridden: %s", len(empty), "y" if len(empty) == 1 else "ies",
                    ", ".join(empty[:10]))

    # The ASV catalogue against amplicon_len. The read counts above say nothing about it:
    # a run can retain every read and still denoise to a region that is not the one
    # configured, or carry a tail of ASVs at twice the amplicon.
    lengths = asv_lengths(export_rep_seqs(rep_seqs, os.path.join(outdir, "seqs"),
                                          args.env, args.timeout))
    dist = check_asv_lengths(lengths, v["amplicon_len"], args.asv_len_tolerance)
    write_asv_lengths(os.path.join(outdir, "asv_lengths.tsv"), lengths, dist)
    write_length_flags(os.path.join(outdir, "asv_length_flags.tsv"), dist)
    report_lengths(dist)

    pooled = [(b, a) for b, a, _ in STEPS if lost_fraction(totals[b], totals[a]) > args.majority]
    if pooled:
        log.warning("pooled across all samples, more than %.0f%% of reads are lost at: %s. %s",
                    args.majority * 100, ", ".join(f"{b}->{a}" for b, a in pooled), QUOTE_STEP)
    if flags:
        per_step = {}
        for f in flags:
            per_step[f["step"]] = per_step.get(f["step"], 0) + 1
        log.warning("%d sample-level flag(s) across %d sample(s): %s. See dada2_flags.tsv",
                    len(flags), len({f["sample-id"] for f in flags}),
                    ", ".join(f"{k} in {n}" for k, n in sorted(per_step.items())))
        if any(f["step"] == "merged->non-chimeric" for f in flags):
            log.warning("%s", QUOTE_CHIMERA)
    else:
        log.info("no sample crossed the %.0f%% loss flag at any step after filtering",
                 args.majority * 100)
    log.info("done: %s, %s, %s, %s", table, rep_seqs,
             os.path.join(outdir, "dada2_stats.tsv"),
             os.path.join(outdir, "asv_lengths.tsv"))


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        run(args)
    except Dada2Error as exc:
        if not log.handlers:
            setup_logging(None)
        log.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
