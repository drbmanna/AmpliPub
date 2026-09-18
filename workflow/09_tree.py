#!/usr/bin/env python3
"""Build the phylogeny, and prove its tips are the features in the table.

`qiime phylogeny align-to-tree-mafft-fasttree` does four things in one call: align with
MAFFT, mask the alignment, build a tree with FastTree, and midpoint root it. This stage
runs it with every parameter stated, then checks the two things that silently go wrong.

**The tips must match the table.** A phylogenetic diversity metric needs every feature in
the table to be a tip in the tree. Filter the table after building the tree, or build the
tree from unfiltered sequences, and UniFrac either fails much later with an opaque message
or quietly computes on a subset. This stage compares the two sets and says exactly which
features are missing.

**The masking step can remove most of the alignment.** `--p-mask-max-gap-frequency` and
`--p-mask-min-conservation` decide how much survives. On badly aligned input the masked
alignment can collapse to a fraction of its length, and the tree is then built on almost
nothing without any error. The length before and after masking is reported so that
collapse is visible rather than invisible.

Standard library only. Needs Python 3.8 or later, Linux or WSL, and conda with the
QIIME 2 environment (see workflow/setup_envs.sh).

Example:
    python workflow/09_tree.py -r ~/research/baxter2016/q2/filtered/rep_seqs_filtered.qza \\
        -b ~/research/baxter2016/q2/filtered/table_filtered.qza \\
        -o ~/research/baxter2016/q2/tree --threads 4
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

MASK_MAX_GAP = 1.0      # qiime default
MASK_MIN_CONSERVATION = 0.4  # qiime default
MIN_MASKED_FRACTION = 0.25   # our choice: below this the alignment has collapsed

log = logging.getLogger("tree")


class TreeError(RuntimeError):
    """A check failed. The message says which one and why."""


def run_cmd(cmd: list[str], timeout: int) -> tuple[int, str, str]:
    """Run a command with stdin closed and a hard time limit."""
    log.info("run: %s", shlex.join(cmd))
    try:
        proc = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True, start_new_session=True)
    except FileNotFoundError as exc:
        raise TreeError(f"cannot start {cmd[0]}: {exc}") from exc
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.communicate()
        raise TreeError(f"timed out after {timeout} s and was killed: {shlex.join(cmd)}")
    return proc.returncode, out, err


def _tail(text: str, n: int = 15) -> str:
    return "\n".join(text.strip().splitlines()[-n:])


def read_from_artifact(qza: str, suffix: str) -> str:
    try:
        with zipfile.ZipFile(qza) as zf:
            hits = [n for n in zf.namelist() if n.endswith(suffix)]
            if len(hits) != 1:
                raise TreeError(f"{qza}: expected one {suffix}, found {len(hits)}")
            return zf.read(hits[0]).decode("utf-8")
    except zipfile.BadZipFile as exc:
        raise TreeError(f"{qza} is not a readable artifact: {exc}") from exc


def newick_tips(text: str) -> set[str]:
    """Tip labels from a Newick tree, excluding internal node labels.

    A label that follows a closing parenthesis belongs to the internal node that
    parenthesis closed. A label anywhere else is a leaf. FastTree writes support values
    exactly where an internal label goes, so a parser that takes every token reports
    them as tips: on this project's tree that turned 654 real tips into 973 and invented
    319 "extra tips" with names like 0.000 and 0.093. Position is the only thing that
    separates the two, so position is what this checks.
    """
    text = text.strip()
    if not text or "(" not in text:
        raise TreeError("the tree is empty or not Newick")
    tips: set[str] = set()
    buf: list[str] = []
    starts_after_close = False
    previous = ""

    def flush() -> None:
        token = "".join(buf).strip()
        buf.clear()
        if not token or starts_after_close:
            return
        label = token.split(":")[0].strip().strip("'\"")
        if label:
            tips.add(label)

    for char in text:
        if char in "(),;":
            flush()
            previous = char
            starts_after_close = char == ")"
            continue
        if not buf and not char.isspace():
            # remember whether this label began right after a ')'
            starts_after_close = previous == ")"
        buf.append(char)
    flush()
    if not tips:
        raise TreeError("no tip labels found in the tree")
    return tips


def fasta_lengths(text: str) -> tuple[int, int]:
    """Number of records and the aligned length, which must be the same for all."""
    lengths = []
    current = 0
    n = 0
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if n:
                lengths.append(current)
            n += 1
            current = 0
        else:
            current += len(line)
    if n:
        lengths.append(current)
    if not lengths:
        raise TreeError("the alignment holds no sequences")
    if len(set(lengths)) > 1:
        raise TreeError(f"the alignment is ragged: lengths {sorted(set(lengths))[:5]}")
    return n, lengths[0]


def table_features(qza: str, outdir: str, env: str, timeout: int) -> set[str]:
    """Feature ids in the table, which the tree tips must cover."""
    dest = os.path.join(outdir, "table_export")
    rc, _, err = run_cmd(["conda", "run", "-n", env, "qiime", "tools", "export",
                          "--input-path", qza, "--output-path", dest], timeout)
    if rc != 0:
        raise TreeError(f"export failed on {qza} with code {rc}:\n{_tail(err)}")
    biom = os.path.join(dest, "feature-table.biom")
    if not os.path.isfile(biom):
        raise TreeError(f"export finished but wrote no {biom}")
    tsv = os.path.join(dest, "feature-table.tsv")
    rc, _, err = run_cmd(["conda", "run", "-n", env, "biom", "convert", "-i", biom,
                          "-o", tsv, "--to-tsv"], timeout)
    if rc != 0:
        raise TreeError(f"biom convert exited with code {rc}:\n{_tail(err)}")
    features = set()
    with open(tsv) as fh:
        for i, line in enumerate(fh):
            if line.startswith("# Const") or not line.strip():
                continue
            first = line.split("\t")[0]
            if first.startswith("#"):
                continue
            features.add(first.strip())
    if not features:
        raise TreeError(f"{tsv}: no features")
    return features


def write_tips(path: str, tips: set[str], features: set[str] | None) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["id", "in_tree", "in_table"])
        for value in sorted(tips | (features or set())):
            w.writerow([value, "yes" if value in tips else "no",
                        "" if features is None else
                        ("yes" if value in features else "no")])


def setup_logging(outdir: str | None) -> None:
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(fmt)
    log.addHandler(console)
    if outdir:
        fh = logging.FileHandler(os.path.join(outdir, "tree_log.txt"), mode="a")
        fh.setFormatter(fmt)
        log.addHandler(fh)


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Align, mask, build and root a tree, then prove its tips match "
                    "the table.")
    p.add_argument("-r", "--rep-seqs", required=True, help="ASV sequences")
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("-b", "--table", help="feature table, to check the tips cover it")
    p.add_argument("--threads", type=int, default=1,
                   help="MAFFT and FastTree threads. QIIME 2's own default is 1, so "
                        "this is always passed explicitly (default 1)")
    p.add_argument("--mask-max-gap-frequency", type=float, default=MASK_MAX_GAP,
                   help=f"QIIME 2 default {MASK_MAX_GAP}")
    p.add_argument("--mask-min-conservation", type=float, default=MASK_MIN_CONSERVATION,
                   help=f"QIIME 2 default {MASK_MIN_CONSERVATION}")
    p.add_argument("--min-masked-fraction", type=float, default=MIN_MASKED_FRACTION,
                   help="fail if masking leaves less than this fraction of the "
                        f"alignment (default {MIN_MASKED_FRACTION}, our choice)")
    p.add_argument("--env", default="amplipub-qiime2-2025.7", help="conda env with QIIME 2")
    p.add_argument("--timeout", type=int, default=86400,
                   help="seconds before a step is killed (default 86400)")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p.parse_args(argv)


def run(args) -> None:
    rep_seqs = os.path.abspath(os.path.expanduser(args.rep_seqs))
    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    if not os.path.isfile(rep_seqs) or os.path.getsize(rep_seqs) == 0:
        raise TreeError(f"rep-seqs not found or empty: {rep_seqs}")
    if args.threads < 1:
        raise TreeError(f"--threads must be at least 1, got {args.threads}")
    for name in ("mask_max_gap_frequency", "mask_min_conservation",
                 "min_masked_fraction"):
        value = getattr(args, name)
        if not 0 <= value <= 1:
            raise TreeError(f"--{name.replace('_', '-')} must be between 0 and 1, "
                            f"got {value}")
    table = None
    if args.table:
        table = os.path.abspath(os.path.expanduser(args.table))
        if not os.path.isfile(table) or os.path.getsize(table) == 0:
            raise TreeError(f"table not found or empty: {table}")
    os.makedirs(outdir, exist_ok=True)
    setup_logging(outdir)
    log.info("=" * 70)
    log.info("tree %s | Python %s | %s", __version__, platform.python_version(),
             platform.platform())
    log.info("command: %s", shlex.join([sys.executable] + sys.argv))
    log.info("started: %s", datetime.now(timezone.utc).isoformat(timespec="seconds"))

    alignment = os.path.join(outdir, "alignment.qza")
    masked = os.path.join(outdir, "masked_alignment.qza")
    unrooted = os.path.join(outdir, "unrooted_tree.qza")
    rooted = os.path.join(outdir, "rooted_tree.qza")
    rc, _, err = run_cmd(
        ["conda", "run", "-n", args.env, "qiime", "phylogeny",
         "align-to-tree-mafft-fasttree", "--i-sequences", rep_seqs,
         "--p-n-threads", str(args.threads),
         "--p-mask-max-gap-frequency", str(args.mask_max_gap_frequency),
         "--p-mask-min-conservation", str(args.mask_min_conservation),
         "--o-alignment", alignment, "--o-masked-alignment", masked,
         "--o-tree", unrooted, "--o-rooted-tree", rooted], args.timeout)
    if rc != 0:
        raise TreeError(f"align-to-tree-mafft-fasttree exited with code {rc}:"
                        f"\n{_tail(err)}")
    for path in (alignment, masked, unrooted, rooted):
        if not os.path.isfile(path):
            raise TreeError(f"the pipeline finished but wrote no {path}")

    n_aln, len_aln = fasta_lengths(read_from_artifact(alignment,
                                                      "/data/aligned-dna-sequences.fasta"))
    n_mask, len_mask = fasta_lengths(read_from_artifact(masked,
                                                        "/data/aligned-dna-sequences.fasta"))
    kept = len_mask / len_aln if len_aln else 0.0
    log.info("aligned %d sequences to %d columns; masking left %d columns (%.1f%%)",
             n_aln, len_aln, len_mask, 100 * kept)
    if n_aln != n_mask:
        raise TreeError(f"the alignment holds {n_aln} sequences but the masked "
                        f"alignment holds {n_mask}")
    if kept < args.min_masked_fraction:
        raise TreeError(
            f"masking left only {100 * kept:.1f}% of the alignment ({len_mask} of "
            f"{len_aln} columns), below the {100 * args.min_masked_fraction:.0f}% floor. "
            "A tree built on this would be noise, and nothing about that raises an "
            "error. Check the sequences are all the same region and orientation")

    tips = newick_tips(read_from_artifact(rooted, "/data/tree.nwk"))
    log.info("rooted tree: %d tip label(s)", len(tips))

    features = None
    if table:
        features = table_features(table, outdir, args.env, args.timeout)
        missing = sorted(features - tips)
        if missing:
            raise TreeError(
                f"{len(missing)} feature(s) in the table are not tips in the tree, for "
                f"example {missing[0]}. A phylogenetic diversity metric needs every "
                "feature to be a tip. This happens when the table is filtered after the "
                "tree is built; build the tree from the filtered sequences instead")
        extra = sorted(tips - features)
        if extra:
            log.warning("%d tip(s) are not in the table: %s%s. Harmless for UniFrac, "
                        "but it means the tree was built from a larger set of sequences",
                        len(extra), ", ".join(extra[:5]),
                        " ..." if len(extra) > 5 else "")
        log.info("every one of the %d table features is a tip in the tree",
                 len(features))
    else:
        log.warning("no --table given, so nothing checked that the tips cover the "
                    "features you will actually analyse. That check is the reason this "
                    "stage exists")

    write_tips(os.path.join(outdir, "tree_tips.tsv"), tips, features)
    log.info("done: %s", rooted)


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        run(args)
    except TreeError as exc:
        if not log.handlers:
            setup_logging(None)
        log.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
