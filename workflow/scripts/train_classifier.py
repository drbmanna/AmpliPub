#!/usr/bin/env python3
"""Train a region-specific naive Bayes classifier for one primer pair.

Why this is a preparation tool and not a pipeline stage. The workflow's `classifier` rule
takes a classifier that already exists and verifies its sha256, so a run always states
exactly which artifact produced its taxonomy. Training inside the pipeline would replace
that pinned artifact with a derived binary rebuilt per run, and would put a slow,
memory-hungry step in the path of every analysis. nf-core/ampliseq does train per run,
because Nextflow caches the result; our contract is different on purpose. So this builds a
classifier once, prints its checksum, and you paste that into the config.

Why train at all, rather than using a full-length classifier for every region. Measured on
Baxter's 8891 V4 ASVs with the two GG2 2024.09 classifiers, same reference and same sklearn,
at confidence 0.7: the full-length classifier named a genus for 68.0% of ASVs against 77.4%
for the region-specific one, 9.4 points fewer, and 94.0% against 97.0% of reads. Where both
named a genus they agreed 98.1% of the time, so the cost is not wrong calls, it is calls not
made: a model trained on whole sequences produces flatter posteriors on a short fragment and
fewer clear the confidence threshold. That measurement is for 253 bp of V4. A longer region
carries more information and may close the gap, so compare rather than assume.

Both steps are QIIME 2 commands, run in the configured environment:
  qiime feature-classifier extract-reads
  qiime feature-classifier fit-classifier-naive-bayes

Standard library only.
"""

from __future__ import annotations

import argparse
import hashlib
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

# Type each input must have. Checked before any compute, because fitting a classifier on the
# wrong artifact wastes tens of minutes and then fails with a message about something else.
SEQ_TYPE = "FeatureData[Sequence]"
TAX_TYPE = "FeatureData[Taxonomy]"
CLS_TYPE = "TaxonomicClassifier"

log = logging.getLogger("train_classifier")


class TrainError(RuntimeError):
    """A check failed. The message says which one and why."""


def run_cmd(cmd: list[str], timeout: int) -> tuple[int, str, str]:
    """Run a command with stdin closed and a hard time limit."""
    log.info("run: %s", shlex.join(cmd))
    try:
        proc = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True, start_new_session=True)
    except FileNotFoundError as exc:
        raise TrainError(f"cannot start {cmd[0]}: {exc}") from exc
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except (AttributeError, ProcessLookupError, PermissionError):
            proc.kill()
        proc.wait()
        raise TrainError(f"timed out after {timeout}s: {shlex.join(cmd)}") from None
    return proc.returncode, out, err


def _tail(text: str, n: int = 15) -> str:
    lines = [line for line in (text or "").splitlines() if line.strip()]
    return "\n".join(lines[-n:])


def sha256_of(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def artifact_meta(path: str) -> dict[str, str]:
    """Read type and uuid out of an artifact's metadata.yaml, without QIIME 2."""
    try:
        with zipfile.ZipFile(path) as zf:
            hits = [n for n in zf.namelist()
                    if n.endswith("/metadata.yaml") and n.count("/") == 2]
            if not hits:
                hits = [n for n in zf.namelist() if n.endswith("/metadata.yaml")]
            if not hits:
                raise TrainError(f"{path}: no metadata.yaml, is this a QIIME 2 artifact?")
            text = zf.read(sorted(hits)[0]).decode("utf-8")
    except zipfile.BadZipFile as exc:
        raise TrainError(f"{path} is not a readable artifact: {exc}") from exc
    meta = {}
    for line in text.splitlines():
        if ":" in line:
            key, _, value = line.partition(":")
            meta[key.strip()] = value.strip().strip("'\"")
    return meta


def check_artifact(label: str, path: str, expected: str) -> dict[str, str]:
    """Refuse anything that is not the artifact type this step needs."""
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise TrainError(f"{label}: not found or empty: {path}")
    meta = artifact_meta(path)
    kind = meta.get("type", "")
    if expected not in kind:
        raise TrainError(
            f"{label} is a {kind or 'unknown type'}, not a {expected}. Check {path}")
    log.info("%s: %s, uuid %s, %d MB", label, kind, meta.get("uuid", "?"),
             round(os.path.getsize(path) / 1e6))
    return meta


def qiime(args, subcommand: list[str]) -> None:
    rc, _, err = run_cmd(["conda", "run", "-n", args.env, "qiime", *subcommand], args.timeout)
    if rc != 0:
        raise TrainError(f"qiime {' '.join(subcommand[:2])} exited with code {rc}:\n"
                         f"{_tail(err)}")


def write_provenance(path: str, rows: list[tuple[str, object]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        for key, value in rows:
            fh.write(f"{key}\t{value}\n")


def setup_logging(outdir: str | None) -> None:
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(fmt)
    log.addHandler(console)
    if outdir:
        fh = logging.FileHandler(os.path.join(outdir, "train_classifier_log.txt"), mode="a")
        fh.setFormatter(fmt)
        log.addHandler(fh)


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--sequences", required=True,
                   help=f"reference sequences, a {SEQ_TYPE} artifact")
    p.add_argument("--taxonomy", required=True,
                   help=f"reference taxonomy, a {TAX_TYPE} artifact")
    p.add_argument("--forward", required=True, help="forward primer")
    p.add_argument("--reverse", required=True, help="reverse primer")
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("--name", required=True,
                   help="name for the output, e.g. gg2_2024.09_v3v4. Becomes the filename "
                        "and the classifier name to put in the config")
    p.add_argument("--n-jobs", type=int, default=1,
                   help="jobs for extract-reads (default 1)")
    p.add_argument("--keep-reads", action="store_true",
                   help="keep the extracted per-region reads artifact, which is the input "
                        "to training and is useful for checking what the primers cut")
    p.add_argument("--env", default="amplipub-qiime2-2025.7", help="conda env with QIIME 2")
    p.add_argument("--timeout", type=int, default=86400,
                   help="seconds before a qiime call is killed (default 86400)")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p.parse_args(argv)


def run(args) -> None:
    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    os.makedirs(outdir, exist_ok=True)
    setup_logging(outdir)
    log.info("=" * 70)
    log.info("train_classifier %s | Python %s | %s", __version__,
             platform.python_version(), platform.platform())
    log.info("command: %s", shlex.join([sys.executable] + sys.argv))
    log.info("started: %s", datetime.now(timezone.utc).isoformat(timespec="seconds"))

    seqs = os.path.abspath(os.path.expanduser(args.sequences))
    tax = os.path.abspath(os.path.expanduser(args.taxonomy))
    seq_meta = check_artifact("sequences", seqs, SEQ_TYPE)
    tax_meta = check_artifact("taxonomy", tax, TAX_TYPE)

    reads = os.path.join(outdir, f"{args.name}_ref_reads.qza")
    out = os.path.join(outdir, f"{args.name}.qza")
    for path in (reads, out):
        if os.path.exists(path):
            raise TrainError(f"REFUSING: {path} already exists. Move it or choose another "
                             "--name, so a classifier is never silently replaced")

    log.info("extracting the %s / %s region from the reference", args.forward, args.reverse)
    qiime(args, ["feature-classifier", "extract-reads",
                 "--i-sequences", seqs,
                 "--p-f-primer", args.forward,
                 "--p-r-primer", args.reverse,
                 "--p-n-jobs", str(args.n_jobs),
                 "--o-reads", reads])
    if not os.path.isfile(reads):
        raise TrainError(f"extract-reads finished but wrote no {reads}")
    check_artifact("extracted reads", reads, SEQ_TYPE)

    log.info("fitting the classifier, this is the slow step")
    qiime(args, ["feature-classifier", "fit-classifier-naive-bayes",
                 "--i-reference-reads", reads,
                 "--i-reference-taxonomy", tax,
                 "--o-classifier", out])
    if not os.path.isfile(out):
        raise TrainError(f"fit-classifier-naive-bayes finished but wrote no {out}")
    # The point of checking our own output: a file of the right size and the wrong type
    # would be found only by the run that tried to use it, hours later.
    out_meta = check_artifact("classifier", out, CLS_TYPE)

    digest = sha256_of(out)
    write_provenance(os.path.join(outdir, f"{args.name}_provenance.tsv"), [
        ("name", args.name),
        ("classifier_path", out),
        ("classifier_sha256", digest),
        ("classifier_uuid", out_meta.get("uuid", "")),
        ("forward_primer", args.forward),
        ("reverse_primer", args.reverse),
        ("reference_sequences", seqs),
        ("reference_sequences_sha256", sha256_of(seqs)),
        ("reference_sequences_uuid", seq_meta.get("uuid", "")),
        ("reference_taxonomy", tax),
        ("reference_taxonomy_sha256", sha256_of(tax)),
        ("reference_taxonomy_uuid", tax_meta.get("uuid", "")),
        ("env", args.env),
        ("command", shlex.join([sys.executable] + sys.argv)),
        ("finished", datetime.now(timezone.utc).isoformat(timespec="seconds")),
    ])

    if not args.keep_reads:
        os.remove(reads)
        log.info("removed the intermediate reads artifact (--keep-reads to keep it)")

    log.info("done: %s", out)
    log.info("put this in the config:")
    log.info("  classifier:")
    log.info("    name: %s", args.name)
    log.info("    path: %s", out)
    log.info("    sha256: %s", digest)


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        run(args)
    except TrainError as exc:
        if not log.handlers:
            setup_logging(None)
        log.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
