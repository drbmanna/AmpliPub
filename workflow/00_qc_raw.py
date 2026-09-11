#!/usr/bin/env python3
"""Run FastQC on every raw FASTQ and combine the reports with MultiQC.

Give it the FASTQ directory written by 00_fetch_sra.py. The script

  1. finds the FASTQ files and checks none are empty,
  2. runs FastQC on all of them in the QC conda environment,
  3. checks one FastQC report came back per input file,
  4. checks R1 and R2 of every pair hold the same number of reads,
  5. writes a per-file table of read counts and FastQC module results,
  6. runs MultiQC over the FastQC reports,
  7. logs the command, versions and results.

Four FastQC modules flag nearly every amplicon library, because all reads start with
the same primer and come from a mixed community. They are marked as expected in the
outputs. The modules worth acting on are quality, N content, length and adapters.

Standard library only. Needs Python 3.8 or later, Linux or WSL, and conda with the
amplipub-qc environment (see workflow/setup_envs.sh).

Example:
    python workflow/00_qc_raw.py -i ~/research/baxter2016/raw/fastq \\
        -o ~/research/baxter2016/qc_raw
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
from collections import defaultdict
from datetime import datetime, timezone

__version__ = "0.1.0"

FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
# Fail on nearly every amplicon library by design: identical primer starts, a mixed
# community instead of one genome, and many identical reads from abundant taxa.
EXPECTED_FOR_AMPLICONS = (
    "Per base sequence content",
    "Per sequence GC content",
    "Sequence Duplication Levels",
    "Overrepresented sequences",
)
STATUSES = ("PASS", "WARN", "FAIL")
PAIR_RE = re.compile(r"^(?P<run>.+)_(?P<mate>[12])$")
TOTAL_RE = re.compile(r"^Total Sequences\t(\d+)\s*$", re.MULTILINE)
MULTIQC_COMMENT = (
    "Per base sequence content, per sequence GC content, sequence duplication and "
    "overrepresented sequences are expected to fail for amplicon libraries. Act on "
    "quality, N content, length and adapter content."
)

log = logging.getLogger("qc_raw")


class QCError(RuntimeError):
    """A check failed. The message says which one and why."""


def stem(name: str) -> str:
    for suffix in FASTQ_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def find_fastqs(indir: str) -> list[str]:
    if not os.path.isdir(indir):
        raise QCError(f"input directory not found: {indir}")
    files = sorted(os.path.join(indir, n) for n in os.listdir(indir)
                   if n.endswith(FASTQ_SUFFIXES) and os.path.isfile(os.path.join(indir, n)))
    # FastQC given no files opens its GUI and waits for ever instead of failing.
    if not files:
        raise QCError(f"no FASTQ files ({', '.join(FASTQ_SUFFIXES)}) in {indir}")
    empty = [f for f in files if os.path.getsize(f) == 0]
    if empty:
        raise QCError(f"{len(empty)} empty FASTQ file(s), e.g. {os.path.basename(empty[0])}")
    seen: dict[str, str] = {}
    for f in files:
        s = stem(os.path.basename(f))
        if s in seen:
            raise QCError(f"{os.path.basename(seen[s])} and {os.path.basename(f)} would write "
                          "the same FastQC report")
        seen[s] = f
    return files


def run_cmd(cmd: list[str], timeout: int) -> tuple[int, str, str]:
    """Run a command with stdin closed and a hard time limit.

    The command gets its own process group, so a timeout also kills what it started
    (conda run starts the tool, FastQC starts java).
    """
    log.info("run: %s", shlex.join(cmd))
    try:
        proc = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True, start_new_session=True)
    except FileNotFoundError as exc:
        raise QCError(f"cannot start {cmd[0]}: {exc}") from exc
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.communicate()
        raise QCError(f"timed out after {timeout} s and was killed: {shlex.join(cmd)}")
    return proc.returncode, out, err


def _tail(text: str, n: int = 10) -> str:
    return "\n".join(text.strip().splitlines()[-n:])


def tool_version(env: str, tool: str) -> str:
    rc, out, err = run_cmd(["conda", "run", "-n", env, tool, "--version"], timeout=180)
    if rc != 0 or not out.strip():
        raise QCError(f"{tool} not usable in conda env {env}. Run workflow/setup_envs.sh.\n"
                      f"{_tail(err)}")
    return out.strip().splitlines()[0]


def run_fastqc(files: list[str], outdir: str, env: str, threads: int, timeout: int) -> None:
    cmd = ["conda", "run", "-n", env, "fastqc", "-q", "--noextract",
           "-t", str(threads), "-o", outdir, *files]
    rc, _, err = run_cmd(cmd, timeout)
    if rc != 0:
        raise QCError(f"FastQC exited with code {rc}:\n{_tail(err)}")


def _member(zf: zipfile.ZipFile, name: str) -> str:
    hits = [n for n in zf.namelist() if n.endswith("/" + name)]
    if len(hits) != 1:
        raise QCError(f"{zf.filename}: expected one {name}, found {len(hits)}")
    return zf.read(hits[0]).decode("utf-8")


def read_report(zip_path: str) -> dict:
    try:
        with zipfile.ZipFile(zip_path) as zf:
            summary = _member(zf, "summary.txt")
            data = _member(zf, "fastqc_data.txt")
    except zipfile.BadZipFile as exc:
        raise QCError(f"{zip_path} is not a readable zip: {exc}") from exc
    modules: dict[str, str] = {}
    filename = None
    for line in summary.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) != 3 or parts[0] not in STATUSES:
            raise QCError(f"{zip_path}: unexpected summary line {line!r}")
        modules[parts[1]] = parts[0]
        filename = parts[2]
    if not modules:
        raise QCError(f"{zip_path}: summary.txt is empty")
    m = TOTAL_RE.search(data)
    if not m:
        raise QCError(f"{zip_path}: no 'Total Sequences' in fastqc_data.txt")
    return {"filename": filename, "total_sequences": int(m.group(1)), "modules": modules}


def collect_reports(files: list[str], fastqc_dir: str) -> list[dict]:
    reports, missing = [], []
    for f in files:
        name = os.path.basename(f)
        zip_path = os.path.join(fastqc_dir, stem(name) + "_fastqc.zip")
        if not os.path.isfile(zip_path):
            missing.append(name)
            continue
        report = read_report(zip_path)
        if report["filename"] != name:
            raise QCError(f"{zip_path} describes {report['filename']}, expected {name}")
        report["file"] = name
        reports.append(report)
    if missing:
        raise QCError(f"FastQC finished but {len(missing)} of {len(files)} reports are "
                      f"missing, e.g. {missing[0]}")
    return reports


def check_pairs(reports: list[dict]) -> int:
    """R1 and R2 of a run must hold the same number of reads. Returns the pair count."""
    mates: dict[str, dict[str, int]] = defaultdict(dict)
    for r in reports:
        m = PAIR_RE.match(stem(r["file"]))
        if m:
            mates[m.group("run")][m.group("mate")] = r["total_sequences"]
    lone = sorted(run for run, d in mates.items() if len(d) != 2)
    if lone:
        raise QCError(f"{len(lone)} run(s) have only one mate, e.g. {lone[0]}")
    bad = sorted(run for run, d in mates.items() if d["1"] != d["2"])
    if bad:
        d = mates[bad[0]]
        raise QCError(f"{len(bad)} pair(s) have different read counts in R1 and R2, "
                      f"e.g. {bad[0]}: {d['1']} vs {d['2']}")
    return len(mates)


def flagged(report: dict, expected: bool) -> list[str]:
    """Non-PASS modules, either the ones expected for amplicons or the ones to check."""
    return [m for m, s in report["modules"].items()
            if s != "PASS" and (m in EXPECTED_FOR_AMPLICONS) == expected]


def module_order(reports: list[dict]) -> list[str]:
    order: list[str] = []
    for r in reports:
        for m in r["modules"]:
            if m not in order:
                order.append(m)
    return order


def write_summary(path: str, reports: list[dict]) -> None:
    modules = module_order(reports)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["file", "total_sequences", *modules, "to_check", "expected_for_amplicons"])
        for r in reports:
            w.writerow([r["file"], r["total_sequences"],
                        *(r["modules"].get(m, "") for m in modules),
                        ";".join(flagged(r, expected=False)),
                        ";".join(flagged(r, expected=True))])


def report_modules(reports: list[dict]) -> int:
    """Log module results. Returns how many files have a module worth checking."""
    for m in module_order(reports):
        counts = {s: sum(r["modules"].get(m) == s for r in reports) for s in ("WARN", "FAIL")}
        if not any(counts.values()):
            continue
        note = "expected for amplicons" if m in EXPECTED_FOR_AMPLICONS else "CHECK"
        log.info("%-30s WARN %4d  FAIL %4d  (%s)", m, counts["WARN"], counts["FAIL"], note)
    to_check = [r for r in reports if flagged(r, expected=False)]
    if to_check:
        log.warning("%d of %d files have a module worth checking; see to_check in "
                    "qc_summary.tsv", len(to_check), len(reports))
    return len(to_check)


def run_multiqc(fastqc_dir: str, outdir: str, env: str, timeout: int) -> str:
    report = os.path.join(outdir, "multiqc_report.html")
    cmd = ["conda", "run", "-n", env, "multiqc", "-f", "-q", "-m", "fastqc",
           "-i", "Raw read QC", "-b", MULTIQC_COMMENT,
           "-n", "multiqc_report.html", "-o", outdir, fastqc_dir]
    rc, _, err = run_cmd(cmd, timeout)
    if rc != 0:
        raise QCError(f"MultiQC exited with code {rc}:\n{_tail(err)}")
    if not os.path.isfile(report) or os.path.getsize(report) == 0:
        raise QCError(f"MultiQC finished but wrote no report at {report}")
    return report


def setup_logging(outdir: str | None) -> None:
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(fmt)
    log.addHandler(console)
    if outdir:
        fh = logging.FileHandler(os.path.join(outdir, "qc_log.txt"), mode="a")
        fh.setFormatter(fmt)
        log.addHandler(fh)


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Run FastQC on every raw FASTQ and combine the reports with MultiQC.")
    p.add_argument("-i", "--input", required=True, help="directory holding the FASTQ files")
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("--env", default="amplipub-qc", help="conda env with FastQC and MultiQC")
    p.add_argument("--threads", type=int, default=4, help="files FastQC processes at once "
                   "(default 4, about 250 MB RAM each)")
    p.add_argument("--timeout", type=int, default=7200,
                   help="seconds before FastQC or MultiQC is killed (default 7200)")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p.parse_args(argv)


def run(args) -> None:
    indir = os.path.abspath(os.path.expanduser(args.input))
    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    files = find_fastqs(indir)
    fastqc_dir = os.path.join(outdir, "fastqc")
    os.makedirs(fastqc_dir, exist_ok=True)
    setup_logging(outdir)
    log.info("=" * 70)
    log.info("qc_raw %s | Python %s | %s", __version__, platform.python_version(), platform.platform())
    log.info("command: %s", shlex.join([sys.executable] + sys.argv))
    log.info("started: %s", datetime.now(timezone.utc).isoformat(timespec="seconds"))
    log.info("input: %d FASTQ files, %.2f GB in %s", len(files),
             sum(os.path.getsize(f) for f in files) / 1e9, indir)
    log.info("versions: %s | %s", tool_version(args.env, "fastqc"), tool_version(args.env, "multiqc"))

    run_fastqc(files, fastqc_dir, args.env, args.threads, args.timeout)
    reports = collect_reports(files, fastqc_dir)
    n_pairs = check_pairs(reports)
    summary = os.path.join(outdir, "qc_summary.tsv")
    write_summary(summary, reports)
    reads = [r["total_sequences"] for r in reports]
    log.info("reads per file: min %d, median %d, max %d; %d complete pairs",
             min(reads), sorted(reads)[len(reads) // 2], max(reads), n_pairs)
    report_modules(reports)
    report = run_multiqc(fastqc_dir, outdir, args.env, args.timeout)
    log.info("done: %d files. Table: %s  Report: %s", len(reports), summary, report)


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        run(args)
    except QCError as exc:
        if not log.handlers:
            setup_logging(None)
        log.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
