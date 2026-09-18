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

MIN_Q = 30          # median quality a position must reach to be kept
AMPLICON_LEN = 253  # V4, 515F to 806R, primers excluded
MIN_OVERLAP = 12    # qiime dada2 denoise-paired --p-min-overlap default in 2025.7
MARGIN = 20         # allowance for V4 length variation between taxa
N_SAMPLED = 10000   # qiime demux summarize --p-n default

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


def write_outputs(outdir: str, fwd: dict, rev: dict, f: int, r: int, overlap: int,
                  args) -> None:
    with open(os.path.join(outdir, "quality_profile.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["position", "median_f", "median_r", "count_f", "count_r"])
        for i in range(max(len(fwd["50%"]), len(rev["50%"]))):
            get = lambda rows, k: rows[k][i] if i < len(rows[k]) else ""
            w.writerow([i + 1, get(fwd, "50%"), get(rev, "50%"),
                        get(fwd, "count"), get(rev, "count")])
    with open(os.path.join(outdir, "trunc_len.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for key, value in (("trunc_len_f", f), ("trunc_len_r", r), ("min_q", args.min_q),
                           ("amplicon_len", args.amplicon_len), ("min_overlap", args.min_overlap),
                           ("margin", args.margin), ("expected_overlap", overlap),
                           ("n_sampled", args.n)):
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
    p.add_argument("--n", type=int, default=N_SAMPLED,
                   help=f"reads sampled for the quality profile (default {N_SAMPLED})")
    p.add_argument("--env", default="amplipub-qiime2-2025.7", help="conda env with QIIME 2")
    p.add_argument("--timeout", type=int, default=3600,
                   help="seconds before demux summarize is killed (default 3600)")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p.parse_args(argv)


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

    qzv = os.path.join(outdir, "quality.qzv")
    rc, _, err = run_cmd(["conda", "run", "-n", args.env, "qiime", "demux", "summarize",
                          "--i-data", qza, "--p-n", str(args.n), "--o-visualization", qzv],
                         args.timeout)
    if rc != 0:
        raise QualityError(f"qiime demux summarize exited with code {rc}:\n{_tail(err)}")
    if not os.path.isfile(qzv):
        raise QualityError(f"qiime demux summarize finished but wrote no {qzv}")

    fwd, rev = read_profiles(qzv)
    f = trunc_len(fwd["50%"], args.min_q, "R1")
    r = trunc_len(rev["50%"], args.min_q, "R2")
    log.info("R1: %d positions, median below Q%g first at %s -> trunc-len-f %d",
             len(fwd["50%"]), args.min_q, f + 1 if f < len(fwd["50%"]) else "none", f)
    log.info("R2: %d positions, median below Q%g first at %s -> trunc-len-r %d",
             len(rev["50%"]), args.min_q, r + 1 if r < len(rev["50%"]) else "none", r)
    overlap = check_overlap(f, r, args.amplicon_len, args.min_overlap, args.margin)
    log.info("expected overlap %d bp (DADA2 needs %d; floor with margin %d)",
             overlap, args.min_overlap, args.min_overlap + args.margin)
    for name, rows, length in (("R1", fwd, f), ("R2", rev, r)):
        log.info("%s: %.1f%% of sampled reads reach position %d (DADA2 discards shorter reads)",
                 name, reach(rows["count"], length), length)
    write_outputs(outdir, fwd, rev, f, r, overlap, args)
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
