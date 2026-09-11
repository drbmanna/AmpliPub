#!/usr/bin/env python3
"""Remove primers with q2-cutadapt and report how many reads carried them.

Give it the paired-end demux.qza. The script

  1. reads the sample list from the artifact's MANIFEST,
  2. runs qiime cutadapt trim-paired with both primers anchored at the read start,
  3. parses cutadapt's report for every sample,
  4. checks one report came back per sample and that no reads were lost,
  5. writes a per-sample table of reads in, reads with primer, and reads out,
  6. runs 00_qc_raw.py on the trimmed reads, only if any read was trimmed,
  7. logs the command, versions and results.

Cutadapt always runs. When the reads carry no primers (for example, the submitter
removed them) nothing is cut and the reads come out unchanged, which the report shows.

Primers are anchored (^) so only a full primer at position 1 is cut. Unanchored 5'
primers also match 3 or more bases of the primer's end at the read start, which can clip
primer-free reads. Untrimmed reads are kept unless --discard-untrimmed is given, because
discarding would drop every read of primer-free data.

Standard library only. Needs Python 3.8 or later, Linux or WSL, and conda with the
QIIME 2 environment (see workflow/setup_envs.sh).

Example:
    python workflow/01_primers.py -i ~/research/baxter2016/q2/demux.qza \\
        -o ~/research/baxter2016/q2/primers
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import platform
import re
import shlex
import signal
import subprocess
import sys
import zipfile
from datetime import datetime, timezone

__version__ = "0.1.0"

PRIMER_F = "GTGCCAGCMGCCGCGGTAA"   # 515F, V4
PRIMER_R = "GGACTACHVGGGTWTCTAAT"  # 806R, V4
PRIMER_RE = re.compile(r"^[ACGTRYSWKMBDHVN]{10,}$")
# QIIME 2 prints its own "Command: ..." lines into the same output, so a report can start
# mid-line ("Command: This is cutadapt 5.1 ..."). Split on the marker wherever it is.
REPORT_MARK = "This is cutadapt "
FIELDS = {
    "pairs_in": r"Total read pairs processed:\s+([\d,]+)",
    "r1_with_primer": r"Read 1 with adapter:\s+([\d,]+)",
    "r2_with_primer": r"Read 2 with adapter:\s+([\d,]+)",
    "pairs_out": r"Pairs written \(passing filters\):\s+([\d,]+)",
}
ABSENT_BELOW = 1.0   # percent of reads with a primer, in every sample
PRESENT_ABOVE = 90.0
QC_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "00_qc_raw.py")

log = logging.getLogger("primers")


class PrimerError(RuntimeError):
    """A check failed. The message says which one and why."""


def validate_primer(seq: str) -> str:
    seq = seq.upper()
    if not PRIMER_RE.match(seq):
        raise PrimerError(f"primer {seq!r} must be at least 10 IUPAC bases (ACGTRYSWKMBDHVN)")
    return seq


def read_manifest(qza: str) -> dict[str, str]:
    """Map each forward-read file name in a paired demux artifact to its sample ID."""
    if not os.path.isfile(qza) or os.path.getsize(qza) == 0:
        raise PrimerError(f"input artifact not found or empty: {qza}")
    try:
        with zipfile.ZipFile(qza) as zf:
            hits = [n for n in zf.namelist() if n.endswith("/data/MANIFEST")]
            if len(hits) != 1:
                raise PrimerError(f"{qza}: expected one data/MANIFEST, found {len(hits)}. "
                                  "Is this a demultiplexed paired-end artifact?")
            text = zf.read(hits[0]).decode("utf-8")
    except zipfile.BadZipFile as exc:
        raise PrimerError(f"{qza} is not a QIIME 2 artifact: {exc}") from exc
    lines = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    rows = list(csv.DictReader(lines))
    forward = {r["filename"]: r["sample-id"] for r in rows if r.get("direction") == "forward"}
    reverse = {r["sample-id"] for r in rows if r.get("direction") == "reverse"}
    if not forward:
        raise PrimerError(f"{qza}: MANIFEST lists no forward reads")
    unpaired = sorted(set(forward.values()) ^ reverse)
    if unpaired:
        raise PrimerError(f"{len(unpaired)} sample(s) lack a forward or reverse file, "
                          f"e.g. {unpaired[0]}")
    return forward


def run_cmd(cmd: list[str], timeout: int) -> tuple[int, str, str]:
    """Run a command with stdin closed and a hard time limit.

    The command gets its own process group, so a timeout also kills what it started.
    """
    log.info("run: %s", shlex.join(cmd))
    try:
        proc = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True, start_new_session=True)
    except FileNotFoundError as exc:
        raise PrimerError(f"cannot start {cmd[0]}: {exc}") from exc
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.communicate()
        raise PrimerError(f"timed out after {timeout} s and was killed: {shlex.join(cmd)}")
    return proc.returncode, out, err


def _tail(text: str, n: int = 10) -> str:
    return "\n".join(text.strip().splitlines()[-n:])


def cutadapt_cmd(qza: str, out_qza: str, env: str, primer_f: str, primer_r: str,
                 discard: bool, cores: int) -> list[str]:
    return ["conda", "run", "-n", env, "qiime", "cutadapt", "trim-paired",
            "--i-demultiplexed-sequences", qza,
            "--p-front-f", "^" + primer_f, "--p-front-r", "^" + primer_r,
            "--p-discard-untrimmed" if discard else "--p-no-discard-untrimmed",
            "--p-cores", str(cores), "--o-trimmed-sequences", out_qza, "--verbose"]


def _count(block: str, field: str, sample: str) -> int:
    m = re.search(FIELDS[field], block)
    if not m:
        raise PrimerError(f"cutadapt report for {sample} has no {field!r} line")
    return int(m.group(1).replace(",", ""))


def parse_reports(text: str, r1_to_sample: dict[str, str]) -> tuple[list[dict], str]:
    """One row per sample from cutadapt's reports. Returns rows and the cutadapt version."""
    blocks = text.split(REPORT_MARK)[1:]
    if not blocks:
        raise PrimerError("no cutadapt reports in the output")
    version = blocks[0].split()[0]
    rows, seen = [], set()
    for block in blocks:
        m = re.search(r"Command line parameters: (.+)", block)
        if not m:
            raise PrimerError("a cutadapt report has no command line")
        r1 = os.path.basename(shlex.split(m.group(1))[-2])
        if r1 not in r1_to_sample:
            raise PrimerError(f"cutadapt report for {r1}, which is not in the MANIFEST")
        sample = r1_to_sample[r1]
        if sample in seen:
            raise PrimerError(f"two cutadapt reports for sample {sample}")
        seen.add(sample)
        row = {"sample": sample, **{f: _count(block, f, sample) for f in FIELDS}}
        rows.append(row)
    missing = sorted(set(r1_to_sample.values()) - seen)
    if missing:
        raise PrimerError(f"{len(missing)} of {len(r1_to_sample)} samples have no cutadapt "
                          f"report, e.g. {missing[0]}")
    return sorted(rows, key=lambda r: r["sample"]), version


def check_counts(rows: list[dict], discard: bool) -> None:
    for r in rows:
        if r["pairs_out"] > r["pairs_in"]:
            raise PrimerError(f"{r['sample']}: {r['pairs_out']} pairs out from {r['pairs_in']} in")
    if not discard:
        lost = [r for r in rows if r["pairs_out"] != r["pairs_in"]]
        if lost:
            r = lost[0]
            raise PrimerError(f"{len(lost)} sample(s) lost reads although untrimmed reads are "
                              f"kept, e.g. {r['sample']}: {r['pairs_in']} in, {r['pairs_out']} out")


def pct(part: int, whole: int) -> float:
    return 100.0 * part / whole if whole else 0.0


def write_summary(path: str, rows: list[dict]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "pairs_in", "r1_with_primer", "r2_with_primer", "pairs_out",
                    "pct_r1_with_primer", "pct_r2_with_primer", "pct_pairs_kept"])
        for r in rows:
            w.writerow([r["sample"], r["pairs_in"], r["r1_with_primer"], r["r2_with_primer"],
                        r["pairs_out"], f"{pct(r['r1_with_primer'], r['pairs_in']):.2f}",
                        f"{pct(r['r2_with_primer'], r['pairs_in']):.2f}",
                        f"{pct(r['pairs_out'], r['pairs_in']):.2f}"])


def interpret(rows: list[dict]) -> str:
    """Say in one line whether the reads carried primers."""
    shares = [pct(min(r["r1_with_primer"], r["r2_with_primer"]), r["pairs_in"]) for r in rows]
    top = [pct(max(r["r1_with_primer"], r["r2_with_primer"]), r["pairs_in"]) for r in rows]
    if max(top) < ABSENT_BELOW:
        return (f"primers absent: under {ABSENT_BELOW:g}% of reads carried a primer in every "
                "sample, so the reads are effectively unchanged")
    if min(shares) > PRESENT_ABOVE:
        return f"primers present: over {PRESENT_ABOVE:g}% of both reads trimmed in every sample"
    return ("MIXED: the share of reads with a primer ranges from "
            f"{min(shares):.1f}% to {max(top):.1f}% across samples. Check primer sequences "
            "and whether some libraries were already trimmed")


def check_output(out_qza: str, samples: set[str]) -> None:
    got = set(read_manifest(out_qza).values())
    if got != samples:
        diff = sorted(got ^ samples)
        raise PrimerError(f"trimmed artifact has a different sample set, e.g. {diff[0]}")


def qc_trimmed(out_qza: str, outdir: str, args) -> str:
    """Export the trimmed reads and run 00_qc_raw.py on them. Returns the report path."""
    fastq_dir = os.path.join(outdir, "trimmed_fastq")
    rc, _, err = run_cmd(["conda", "run", "-n", args.env, "qiime", "tools", "export",
                          "--input-path", out_qza, "--output-path", fastq_dir], args.timeout)
    if rc != 0:
        raise PrimerError(f"qiime tools export exited with code {rc}:\n{_tail(err)}")
    qc_dir = os.path.join(outdir, "qc_trimmed")
    # 00_qc_raw applies its own limit to FastQC and to MultiQC.
    rc, _, err = run_cmd([sys.executable, QC_SCRIPT, "-i", fastq_dir, "-o", qc_dir,
                          "--env", args.qc_env, "--timeout", str(args.timeout)],
                         2 * args.timeout + 600)
    if rc != 0:
        raise PrimerError(f"trimmed-read QC failed, see {qc_dir}/qc_log.txt:\n{_tail(err)}")
    return os.path.join(qc_dir, "multiqc_report.html")


def setup_logging(outdir: str | None) -> None:
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(fmt)
    log.addHandler(console)
    if outdir:
        fh = logging.FileHandler(os.path.join(outdir, "primers_log.txt"), mode="a")
        fh.setFormatter(fmt)
        log.addHandler(fh)


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Remove primers with q2-cutadapt and report how many reads carried them.")
    p.add_argument("-i", "--input", required=True, help="paired-end demux.qza")
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("--primer-f", default=PRIMER_F, help=f"forward primer (default 515F {PRIMER_F})")
    p.add_argument("--primer-r", default=PRIMER_R, help=f"reverse primer (default 806R {PRIMER_R})")
    p.add_argument("--discard-untrimmed", action="store_true",
                   help="drop pairs where a primer was not found (only for reads that carry primers)")
    p.add_argument("--cores", type=int, default=4, help="CPU cores for cutadapt (default 4)")
    p.add_argument("--env", default="qiime2-amplicon-2025.7", help="conda env with QIIME 2")
    p.add_argument("--qc-env", default="amplipub-qc",
                   help="conda env with FastQC and MultiQC, for QC of trimmed reads")
    p.add_argument("--timeout", type=int, default=7200,
                   help="seconds before cutadapt is killed (default 7200)")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p.parse_args(argv)


def run(args) -> None:
    qza = os.path.abspath(os.path.expanduser(args.input))
    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    primer_f, primer_r = validate_primer(args.primer_f), validate_primer(args.primer_r)
    r1_to_sample = read_manifest(qza)
    os.makedirs(outdir, exist_ok=True)
    setup_logging(outdir)
    log.info("=" * 70)
    log.info("primers %s | Python %s | %s", __version__, platform.python_version(), platform.platform())
    log.info("command: %s", shlex.join([sys.executable] + sys.argv))
    log.info("started: %s", datetime.now(timezone.utc).isoformat(timespec="seconds"))
    log.info("input: %s, %d samples", qza, len(r1_to_sample))
    log.info("primers: forward ^%s, reverse ^%s; untrimmed pairs %s", primer_f, primer_r,
             "discarded" if args.discard_untrimmed else "kept")

    out_qza = os.path.join(outdir, "trimmed.qza")
    cmd = cutadapt_cmd(qza, out_qza, args.env, primer_f, primer_r, args.discard_untrimmed,
                       args.cores)
    rc, out, err = run_cmd(cmd, args.timeout)
    with open(os.path.join(outdir, "cutadapt_report.log"), "w") as fh:
        fh.write(out)
        fh.write("\n# ---- stderr ----\n")
        fh.write(err)
    if rc != 0:
        raise PrimerError(f"qiime cutadapt exited with code {rc}:\n{_tail(err)}")
    if not os.path.isfile(out_qza) or os.path.getsize(out_qza) == 0:
        raise PrimerError(f"qiime cutadapt finished but wrote no {out_qza}")

    rows, version = parse_reports(out + "\n" + err, r1_to_sample)
    check_counts(rows, args.discard_untrimmed)
    check_output(out_qza, set(r1_to_sample.values()))
    write_summary(os.path.join(outdir, "primer_summary.tsv"), rows)
    pairs_in = sum(r["pairs_in"] for r in rows)
    log.info("cutadapt %s: %d samples, %d pairs in, %d out; R1 with primer %d (%.2f%%), "
             "R2 with primer %d (%.2f%%)", version, len(rows), pairs_in,
             sum(r["pairs_out"] for r in rows),
             sum(r["r1_with_primer"] for r in rows), pct(sum(r["r1_with_primer"] for r in rows), pairs_in),
             sum(r["r2_with_primer"] for r in rows), pct(sum(r["r2_with_primer"] for r in rows), pairs_in))
    verdict = interpret(rows)
    (log.warning if verdict.startswith("MIXED") else log.info)("%s", verdict)
    trimmed = sum(r["r1_with_primer"] + r["r2_with_primer"] for r in rows)
    if trimmed == 0:
        log.info("no reads trimmed: the reads are unchanged, so the raw-read QC applies. "
                 "Trimmed-read QC skipped")
    else:
        log.info("%d reads trimmed: running QC on the trimmed reads", trimmed)
        log.info("trimmed-read QC: %s", qc_trimmed(out_qza, outdir, args))
    log.info("done: %s", out_qza)


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        run(args)
    except PrimerError as exc:
        if not log.handlers:
            setup_logging(None)
        log.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
