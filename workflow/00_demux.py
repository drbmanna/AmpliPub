#!/usr/bin/env python3
"""Demultiplex EMP-protocol reads, and account for every read that goes missing.

Give it a directory holding `forward.fastq.gz`, `reverse.fastq.gz` and
`barcodes.fastq.gz`, plus a metadata file with a barcode column. The script

  1. checks the barcode column before spending any compute on it,
  2. imports the reads and demultiplexes them,
  3. accounts for every read: assigned plus unassigned must equal the input,
  4. reports the barcode error corrections, which QIIME 2 emits and nobody reads,
  5. diagnoses the barcode orientation when assignment comes out low.

**The orientation trap.** `--p-rev-comp-barcodes` and `--p-rev-comp-mapping-barcodes`
both default to False, and getting either wrong does not raise an error. It assigns
almost nothing, and the run continues to a near-empty feature table that looks like a
failed experiment rather than a wrong flag. This is the most common catastrophic failure
in 16S processing. So when the assigned fraction comes out below `--min-assigned`, this
script tries the other three combinations, reports what each one would assign, and stops
with the answer instead of handing back an empty table.

**Read accounting is the point of the stage.** A demultiplexing step that does not tell
you how many reads it threw away is not a quality control step. The unassigned fraction
is a real measurement: it carries barcode quality, index hopping and contamination from
other libraries on the same run.

**Already demultiplexed?** Most public data, including anything from SRA, arrives as one
FASTQ pair per sample. There is nothing to demultiplex, and this stage is not for you:
import with a manifest instead, as `00_fetch_sra.py` prints at the end of its log.

Standard library only. Needs Python 3.8 or later, Linux or WSL, and conda with the
QIIME 2 environment (see workflow/setup_envs.sh).

Example:
    python workflow/00_demux.py -i ~/research/run1/emp -m ~/research/run1/metadata.tsv \\
        --barcode-column barcode-sequence -o ~/research/run1/demux
"""

from __future__ import annotations

import argparse
import csv
import gzip
import itertools
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

# EMP-protocol imports require exactly these file names in the input directory.
EMP_FILES = ("forward.fastq.gz", "reverse.fastq.gz", "barcodes.fastq.gz")
MIN_ASSIGNED = 0.5  # our choice, not a standard; below this the orientation is checked

log = logging.getLogger("demux")


class DemuxError(RuntimeError):
    """A check failed. The message says which one and why."""


def run_cmd(cmd: list[str], timeout: int) -> tuple[int, str, str]:
    """Run a command with stdin closed and a hard time limit."""
    log.info("run: %s", shlex.join(cmd))
    try:
        proc = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True, start_new_session=True)
    except FileNotFoundError as exc:
        raise DemuxError(f"cannot start {cmd[0]}: {exc}") from exc
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.communicate()
        raise DemuxError(f"timed out after {timeout} s and was killed: {shlex.join(cmd)}")
    return proc.returncode, out, err


def _tail(text: str, n: int = 15) -> str:
    return "\n".join(text.strip().splitlines()[-n:])


# ---- inputs --------------------------------------------------------------

def check_input_dir(path: str) -> None:
    """EMP import needs three files with exact names. Say which one is missing."""
    if not os.path.isdir(path):
        raise DemuxError(f"not a directory: {path}")
    missing = [f for f in EMP_FILES if not os.path.isfile(os.path.join(path, f))]
    if missing:
        found = sorted(os.listdir(path))[:8]
        raise DemuxError(
            f"{path} is missing {', '.join(missing)}. An EMP import needs exactly "
            f"{', '.join(EMP_FILES)}. Found: {', '.join(found) or 'nothing'}. If your "
            "data is already one FASTQ pair per sample, it is demultiplexed and this "
            "stage is not the one you want; import with a manifest instead")
    empty = [f for f in EMP_FILES if os.path.getsize(os.path.join(path, f)) == 0]
    if empty:
        raise DemuxError(f"{path}: {', '.join(empty)} is empty")


def read_barcodes(path: str, column: str) -> dict[str, str]:
    """Read sample id to barcode. Duplicates and ragged lengths are fatal here.

    Two samples sharing a barcode cannot be told apart, and mixed barcode lengths mean
    the wrong column or a malformed file. Both produce a plausible-looking result.
    """
    text = open(path, "rb").read().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    rows = list(csv.DictReader(text.decode("utf-8").splitlines(), delimiter="\t"))
    if not rows:
        raise DemuxError(f"{path}: no rows")
    if column not in rows[0]:
        raise DemuxError(f"{path} has no column {column!r}. Columns: "
                         f"{', '.join(list(rows[0])[:12])}")
    id_col = list(rows[0])[0]
    out: dict[str, str] = {}
    for row in rows:
        sample = (row[id_col] or "").strip()
        if not sample or sample.startswith("#"):
            continue
        barcode = (row[column] or "").strip().upper()
        if not barcode:
            raise DemuxError(f"{path}: {sample} has no barcode in {column!r}")
        if sample in out:
            raise DemuxError(f"{path}: sample {sample} appears more than once")
        out[sample] = barcode
    if not out:
        raise DemuxError(f"{path}: no samples with a barcode")
    seen: dict[str, str] = {}
    for sample, barcode in out.items():
        if barcode in seen:
            raise DemuxError(
                f"{path}: {sample} and {seen[barcode]} share the barcode {barcode}. "
                "They cannot be told apart")
        seen[barcode] = sample
    lengths = {len(b) for b in out.values()}
    if len(lengths) > 1:
        raise DemuxError(f"{path}: barcodes have mixed lengths {sorted(lengths)}. "
                         f"Is {column!r} the right column?")
    bad = {b for b in out.values() if set(b) - set("ACGTN")}
    if bad:
        raise DemuxError(f"{path}: barcode {sorted(bad)[0]} is not a DNA sequence. "
                         f"Is {column!r} the right column?")
    return out


def count_reads(path: str) -> int:
    """Reads in a gzipped FASTQ. The honest denominator for read accounting."""
    n = 0
    with gzip.open(path, "rb") as fh:
        for n, _ in enumerate(fh, 1):
            pass
    if n % 4:
        raise DemuxError(f"{path}: {n} lines is not a whole number of FASTQ records")
    return n // 4


# ---- artifacts -----------------------------------------------------------

def read_tsv_from_artifact(qza: str, suffix: str) -> str:
    try:
        with zipfile.ZipFile(qza) as zf:
            hits = [n for n in zf.namelist() if n.endswith(suffix)]
            if len(hits) != 1:
                raise DemuxError(f"{qza}: expected one {suffix}, found {len(hits)}")
            return zf.read(hits[0]).decode("utf-8")
    except zipfile.BadZipFile as exc:
        raise DemuxError(f"{qza} is not a readable artifact: {exc}") from exc


def parse_per_sample_counts(text: str, name: str) -> dict[str, int]:
    """Read `qiime demux tabulate-read-counts` style output into {sample: reads}."""
    rows = list(csv.reader(text.splitlines(), delimiter="\t"))
    if len(rows) < 2:
        raise DemuxError(f"{name}: no counts")
    header = [c.strip() for c in rows[0]]
    try:
        col = next(i for i, c in enumerate(header)
                   if c.lower().replace("_", " ") in ("forward sequence count",
                                                      "sequence count", "count",
                                                      "reads"))
    except StopIteration:
        raise DemuxError(f"{name}: no read count column. Columns: {', '.join(header)}")
    out: dict[str, int] = {}
    for row in rows[1:]:
        if not row or not row[0].strip() or row[0].startswith("#"):
            continue
        try:
            out[row[0].strip()] = int(float(row[col]))
        except (IndexError, ValueError) as exc:
            raise DemuxError(f"{name}: {row[0]} has no usable count") from exc
    if not out:
        raise DemuxError(f"{name}: no samples")
    return out


def parse_error_corrections(text: str) -> dict[str, int]:
    """Summarise the ErrorCorrectionDetails table QIIME 2 emits and nobody reads."""
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    summary = {"records": 0, "corrected": 0, "uncorrectable": 0}
    for row in rows:
        first = (row.get("id") or row.get("") or "").strip()
        if first.startswith("#q2:types"):
            continue
        summary["records"] += 1
        errors = row.get("errors") or row.get("error") or ""
        sample = (row.get("sample") or "").strip()
        try:
            n_err = float(errors) if errors not in ("", None) else 0.0
        except ValueError:
            n_err = 0.0
        if not sample:
            summary["uncorrectable"] += 1
        elif n_err > 0:
            summary["corrected"] += 1
    return summary


# ---- demultiplexing ------------------------------------------------------

def import_emp(indir: str, outdir: str, env: str, timeout: int) -> str:
    qza = os.path.join(outdir, "emp.qza")
    rc, _, err = run_cmd(["conda", "run", "-n", env, "qiime", "tools", "import",
                          "--type", "EMPPairedEndSequences",
                          "--input-path", indir, "--output-path", qza], timeout)
    if rc != 0:
        raise DemuxError(f"qiime tools import exited with code {rc}:\n{_tail(err)}")
    if not os.path.isfile(qza):
        raise DemuxError(f"import finished but wrote no {qza}")
    return qza


def demultiplex(emp: str, metadata: str, column: str, rev_barcodes: bool,
                rev_mapping: bool, golay: bool, outdir: str, tag: str, env: str,
                timeout: int) -> tuple[str, str]:
    """One demultiplexing run. Every flag is passed explicitly, including the defaults."""
    seqs = os.path.join(outdir, f"per_sample_sequences{tag}.qza")
    details = os.path.join(outdir, f"error_correction{tag}.qza")
    cmd = ["conda", "run", "-n", env, "qiime", "demux", "emp-paired",
           "--i-seqs", emp,
           "--m-barcodes-file", metadata, "--m-barcodes-column", column,
           "--p-golay-error-correction" if golay else "--p-no-golay-error-correction",
           "--p-rev-comp-barcodes" if rev_barcodes else "--p-no-rev-comp-barcodes",
           "--p-rev-comp-mapping-barcodes" if rev_mapping
           else "--p-no-rev-comp-mapping-barcodes",
           "--o-per-sample-sequences", seqs,
           "--o-error-correction-details", details]
    rc, _, err = run_cmd(cmd, timeout)
    if rc != 0:
        raise DemuxError(f"qiime demux emp-paired exited with code {rc}:\n{_tail(err)}")
    for path in (seqs, details):
        if not os.path.isfile(path):
            raise DemuxError(f"emp-paired finished but wrote no {path}")
    return seqs, details


def count_assigned(seqs: str, outdir: str, tag: str, env: str,
                   timeout: int) -> dict[str, int]:
    counts = os.path.join(outdir, f"read_counts{tag}.qzv")
    rc, _, err = run_cmd(["conda", "run", "-n", env, "qiime", "demux",
                          "tabulate-read-counts", "--i-sequences", seqs,
                          "--o-visualization", counts], timeout)
    if rc != 0:
        raise DemuxError(f"tabulate-read-counts exited with code {rc}:\n{_tail(err)}")
    text = read_tsv_from_artifact(counts, "/data/metadata.tsv")
    return parse_per_sample_counts(text, counts)


def write_counts(path: str, counts: dict[str, int], barcodes: dict[str, str],
                 total: int) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample-id", "barcode", "reads", "fraction_of_input"])
        for sample in sorted(counts, key=lambda s: -counts[s]):
            w.writerow([sample, barcodes.get(sample, ""), counts[sample],
                        round(counts[sample] / total, 6) if total else 0.0])


def setup_logging(outdir: str | None) -> None:
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(fmt)
    log.addHandler(console)
    if outdir:
        fh = logging.FileHandler(os.path.join(outdir, "demux_log.txt"), mode="a")
        fh.setFormatter(fmt)
        log.addHandler(fh)


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Demultiplex EMP-protocol paired reads and account for every read.")
    p.add_argument("-i", "--input", required=True,
                   help=f"directory holding {', '.join(EMP_FILES)}")
    p.add_argument("-m", "--metadata", required=True, help="sample metadata file")
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("--barcode-column", default="barcode-sequence",
                   help="metadata column holding the barcodes "
                        "(default barcode-sequence)")
    p.add_argument("--rev-comp-barcodes", action="store_true",
                   help="reverse complement the barcode reads")
    p.add_argument("--rev-comp-mapping-barcodes", action="store_true",
                   help="reverse complement the barcodes in the metadata")
    p.add_argument("--no-golay", action="store_true",
                   help="switch off 12nt Golay error correction, which QIIME 2 has on "
                        "by default. Only correct for 12nt Golay barcodes")
    p.add_argument("--min-assigned", type=float, default=MIN_ASSIGNED,
                   help=f"fraction of reads that must be assigned before the barcode "
                        f"orientation is treated as wrong (default {MIN_ASSIGNED}, our "
                        f"choice, not a standard)")
    p.add_argument("--no-orientation-check", action="store_true",
                   help="do not try the other barcode orientations on low assignment")
    p.add_argument("--env", default="amplipub-qiime2-2025.7", help="conda env with QIIME 2")
    p.add_argument("--timeout", type=int, default=86400,
                   help="seconds before a step is killed (default 86400)")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p.parse_args(argv)


def run(args) -> None:
    indir = os.path.abspath(os.path.expanduser(args.input))
    metadata = os.path.abspath(os.path.expanduser(args.metadata))
    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    check_input_dir(indir)
    if not os.path.isfile(metadata) or os.path.getsize(metadata) == 0:
        raise DemuxError(f"metadata not found or empty: {metadata}")
    if not 0 < args.min_assigned <= 1:
        raise DemuxError(f"--min-assigned must be above 0 and at most 1, got "
                         f"{args.min_assigned}")
    os.makedirs(outdir, exist_ok=True)
    setup_logging(outdir)
    log.info("=" * 70)
    log.info("demux %s | Python %s | %s", __version__, platform.python_version(),
             platform.platform())
    log.info("command: %s", shlex.join([sys.executable] + sys.argv))
    log.info("started: %s", datetime.now(timezone.utc).isoformat(timespec="seconds"))

    barcodes = read_barcodes(metadata, args.barcode_column)
    length = len(next(iter(barcodes.values())))
    log.info("%d samples with unique %d nt barcodes in %r", len(barcodes), length,
             args.barcode_column)
    golay = not args.no_golay
    if golay and length != 12:
        log.warning("Golay error correction is on but the barcodes are %d nt, not 12. "
                    "Golay correction is for 12nt barcodes; pass --no-golay unless you "
                    "know these are Golay", length)

    total = count_reads(os.path.join(indir, "barcodes.fastq.gz"))
    log.info("input: %d read pairs in %s", total, indir)

    emp = import_emp(indir, outdir, args.env, args.timeout)
    seqs, details = demultiplex(emp, metadata, args.barcode_column,
                                args.rev_comp_barcodes, args.rev_comp_mapping_barcodes,
                                golay, outdir, "", args.env, args.timeout)
    counts = count_assigned(seqs, outdir, "", args.env, args.timeout)
    assigned = sum(counts.values())
    unassigned = total - assigned
    if unassigned < 0:
        raise DemuxError(
            f"{assigned} reads were assigned but only {total} were read from "
            "barcodes.fastq.gz. The accounting cannot be trusted; check that the three "
            "EMP files come from the same run")
    log.info("assigned %d of %d read pairs (%.2f%%); %d unassigned (%.2f%%)",
             assigned, total, 100 * assigned / total if total else 0.0,
             unassigned, 100 * unassigned / total if total else 0.0)

    fraction = assigned / total if total else 0.0
    if fraction < args.min_assigned and not args.no_orientation_check:
        log.warning("only %.2f%% of reads were assigned, below the %.0f%% you asked for. "
                    "Neither rev-comp flag errors when it is wrong, it just assigns "
                    "almost nothing, so the other orientations are being tried",
                    100 * fraction, 100 * args.min_assigned)
        tried = [(args.rev_comp_barcodes, args.rev_comp_mapping_barcodes, assigned)]
        for rb, rm in itertools.product((False, True), repeat=2):
            if (rb, rm) == (args.rev_comp_barcodes, args.rev_comp_mapping_barcodes):
                continue
            tag = f"_rb{int(rb)}_rm{int(rm)}"
            try:
                alt, _ = demultiplex(emp, metadata, args.barcode_column, rb, rm, golay,
                                     outdir, tag, args.env, args.timeout)
                got = sum(count_assigned(alt, outdir, tag, args.env,
                                         args.timeout).values())
            except DemuxError as exc:
                log.warning("  rev-comp-barcodes=%s rev-comp-mapping-barcodes=%s "
                            "failed: %s", rb, rm, exc)
                continue
            tried.append((rb, rm, got))
            log.warning("  rev-comp-barcodes=%s rev-comp-mapping-barcodes=%s would "
                        "assign %d reads (%.2f%%)", rb, rm, got,
                        100 * got / total if total else 0.0)
        best = max(tried, key=lambda t: t[2])
        if best[2] > assigned:
            raise DemuxError(
                f"the barcode orientation looks wrong. rev-comp-barcodes={best[0]} and "
                f"rev-comp-mapping-barcodes={best[1]} assigns {best[2]} reads "
                f"({100 * best[2] / total:.2f}%) against {assigned} "
                f"({100 * fraction:.2f}%) for the flags you gave. Rerun with those "
                "flags. Nothing downstream is written, because a near-empty table is "
                "the failure this check exists to prevent")
        raise DemuxError(
            f"only {100 * fraction:.2f}% of reads were assigned and no barcode "
            "orientation does better. The barcodes, the metadata column or the files "
            "themselves are the problem, not the orientation")

    corrections = parse_error_corrections(read_tsv_from_artifact(details,
                                                                "/data/metadata.tsv"))
    if corrections["records"]:
        log.info("barcode error correction: %d record(s), %d corrected, %d with no "
                 "sample assigned", corrections["records"], corrections["corrected"],
                 corrections["uncorrectable"])

    write_counts(os.path.join(outdir, "demux_counts.tsv"), counts, barcodes, total)
    empty = sorted(s for s in barcodes if counts.get(s, 0) == 0)
    if empty:
        log.warning("%d sample(s) in the metadata got zero reads: %s%s", len(empty),
                    ", ".join(empty[:10]), " ..." if len(empty) > 10 else "")
    extra = sorted(set(counts) - set(barcodes))
    if extra:
        log.warning("%d demultiplexed sample(s) are not in the metadata: %s",
                    len(extra), ", ".join(extra[:10]))
    values = sorted(counts.values())
    log.info("per-sample reads: min %d, median %d, max %d",
             values[0], values[len(values) // 2], values[-1])
    log.info("the unassigned fraction is a measurement, not waste: it carries barcode "
             "quality, index hopping and contamination from other libraries on the run")
    log.info("done: %s  (%s)", seqs, os.path.join(outdir, "demux_counts.tsv"))


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        run(args)
    except DemuxError as exc:
        if not log.handlers:
            setup_logging(None)
        log.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
