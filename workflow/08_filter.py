#!/usr/bin/env python3
"""Filter the table to the organisms and samples the study is about, and count the cost.

Three filters, in the order that makes each one interpretable, with the reads and
features surviving each step reported:

  1. **Taxonomy.** Keep the target domain, drop mitochondria and chloroplast. Host and
     plant organelle 16S amplifies happily and is not part of a bacterial community.
  2. **Samples.** Drop mocks and controls before any community statistic is computed.
     A mock left in the table shifts every between-sample distance.
  3. **Prevalence.** Drop features seen in too few samples, which is a choice about what
     counts as evidence rather than a technical step, so the number is recorded.

The representative sequences are filtered to match the table, because a tree built from
sequences the table no longer holds is a silent mismatch that surfaces much later as a
diversity metric that cannot be computed.

**The trap this stage was written around.** `qiime feature-table filter-features` has
`--p-filter-empty-samples` on by default, so filtering *features* can silently remove
*samples*. A prevalence filter can therefore change your sample count without saying so.
Every step here reports samples before and after, and losing a sample to a feature filter
is called out rather than absorbed.

**The prevalence threshold is a choice, not a standard.** For the Baxter reproduction the
lab kept OTUs present in at least 5% of samples, which is why `--min-samples-fraction`
exists and why its value lands in `filter_summary.tsv` next to the result.

Standard library only. Needs Python 3.8 or later, Linux or WSL, and conda with the
QIIME 2 environment (see workflow/setup_envs.sh).

Example:
    python workflow/08_filter.py -b ~/research/baxter2016/q2/collapsed/table_by_sample.qza \\
        -r ~/research/baxter2016/q2/dada2/rep_seqs.qza \\
        -t ~/research/baxter2016/q2/taxonomy_final/taxonomy_gg2_2024.09_v4.qza \\
        -m ~/research/baxter2016/q2/collapsed/sample_metadata.tsv \\
        --drop-where "in_study_metadata='no'" --min-samples-fraction 0.05 \\
        -o ~/research/baxter2016/q2/filtered
"""

from __future__ import annotations

import argparse
import csv
import logging
import math
import os
import platform
import shlex
import signal
import subprocess
import sys
from datetime import datetime, timezone

__version__ = "0.1.0"

INCLUDE = "Bacteria"                 # the domain the study is about
EXCLUDE = "mitochondria,chloroplast"  # organelle 16S, not community members

log = logging.getLogger("filter")


class FilterError(RuntimeError):
    """A check failed. The message says which one and why."""


def run_cmd(cmd: list[str], timeout: int) -> tuple[int, str, str]:
    """Run a command with stdin closed and a hard time limit."""
    log.info("run: %s", shlex.join(cmd))
    try:
        proc = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True, start_new_session=True)
    except FileNotFoundError as exc:
        raise FilterError(f"cannot start {cmd[0]}: {exc}") from exc
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.communicate()
        raise FilterError(f"timed out after {timeout} s and was killed: {shlex.join(cmd)}")
    return proc.returncode, out, err


def _tail(text: str, n: int = 15) -> str:
    return "\n".join(text.strip().splitlines()[-n:])


def export_counts(qza: str, outdir: str, tag: str, env: str,
                  timeout: int) -> tuple[dict[str, float], int]:
    """Reads per sample and the feature count, read straight off the table."""
    dest = os.path.join(outdir, f"export_{tag}")
    rc, _, err = run_cmd(["conda", "run", "-n", env, "qiime", "tools", "export",
                          "--input-path", qza, "--output-path", dest], timeout)
    if rc != 0:
        raise FilterError(f"export failed on {qza} with code {rc}:\n{_tail(err)}")
    biom = os.path.join(dest, "feature-table.biom")
    if not os.path.isfile(biom):
        raise FilterError(f"export finished but wrote no {biom}")
    tsv = os.path.join(dest, "feature-table.tsv")
    rc, _, err = run_cmd(["conda", "run", "-n", env, "biom", "convert", "-i", biom,
                          "-o", tsv, "--to-tsv"], timeout)
    if rc != 0:
        raise FilterError(f"biom convert exited with code {rc}:\n{_tail(err)}")
    per_sample: dict[str, float] = {}
    features = 0
    with open(tsv) as fh:
        header = None
        for line in fh:
            if line.startswith("# Const") or not line.strip():
                continue
            cells = line.rstrip("\n").split("\t")
            if header is None:
                header = cells[1:]
                per_sample = {name: 0.0 for name in header}
                continue
            features += 1
            if len(cells) != len(header) + 1:
                raise FilterError(f"{tsv}: row {cells[0]!r} has {len(cells) - 1} values "
                                  f"for {len(header)} columns")
            for name, value in zip(header, cells[1:]):
                try:
                    per_sample[name] += float(value)
                except ValueError as exc:
                    raise FilterError(f"{tsv}: {cells[0]} in {name} is not a "
                                      f"number") from exc
    if header is None:
        raise FilterError(f"{tsv}: no header")
    return per_sample, features


class Step:
    """One filter, with what it cost recorded so the report writes itself."""

    def __init__(self, name: str, detail: str, samples: int, features: int,
                 reads: float):
        self.name, self.detail = name, detail
        self.samples, self.features, self.reads = samples, features, reads


def write_summary(path: str, steps: list[Step], settings: dict[str, object]) -> None:
    with open(path, "w", newline="") as fh:
        fh.write("# Every filter applied, in order, and what it cost.\n")
        fh.write("# Thresholds here are choices recorded with the result, not standards.\n")
        for key, value in settings.items():
            fh.write(f"# {key}: {value}\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["step", "detail", "samples", "features", "reads",
                    "samples_lost", "features_lost", "reads_lost"])
        previous = None
        for step in steps:
            if previous is None:
                w.writerow([step.name, step.detail, step.samples, step.features,
                            int(step.reads), 0, 0, 0])
            else:
                w.writerow([step.name, step.detail, step.samples, step.features,
                            int(step.reads), previous.samples - step.samples,
                            previous.features - step.features,
                            int(previous.reads - step.reads)])
            previous = step


def write_per_sample(path: str, before: dict[str, float],
                     after: dict[str, float]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample-id", "reads_before", "reads_after", "fraction_kept",
                    "dropped"])
        for sample in sorted(before):
            b = before[sample]
            a = after.get(sample, 0.0)
            w.writerow([sample, int(b), int(a),
                        round(a / b, 6) if b else 0.0,
                        "yes" if sample not in after else "no"])


def setup_logging(outdir: str | None) -> None:
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(fmt)
    log.addHandler(console)
    if outdir:
        fh = logging.FileHandler(os.path.join(outdir, "filter_log.txt"), mode="a")
        fh.setFormatter(fmt)
        log.addHandler(fh)


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Filter the table by taxonomy, samples and prevalence, with "
                    "accounting at every step.")
    p.add_argument("-b", "--table", required=True, help="sample-level feature table")
    p.add_argument("-r", "--rep-seqs", required=True, help="ASV sequences")
    p.add_argument("-t", "--taxonomy", required=True, help="taxonomy artifact")
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("-m", "--metadata", help="sample metadata, needed for --drop-where")
    p.add_argument("--include", default=INCLUDE,
                   help=f"taxa to keep, comma separated (default {INCLUDE!r})")
    p.add_argument("--exclude", default=EXCLUDE,
                   help=f"taxa to drop, comma separated (default {EXCLUDE!r})")
    p.add_argument("--drop-where",
                   help="SQLite WHERE clause selecting samples to DROP, e.g. "
                        "\"in_study_metadata='no'\"")
    p.add_argument("--min-samples-fraction", type=float, default=0.0,
                   help="drop features seen in fewer than this fraction of samples. "
                        "A choice about evidence, recorded with the result (default 0, off)")
    p.add_argument("--min-sample-reads", type=int, default=0,
                   help="drop samples below this many reads after filtering "
                        "(default 0, off)")
    p.add_argument("--env", default="amplipub-qiime2-2025.7", help="conda env with QIIME 2")
    p.add_argument("--timeout", type=int, default=86400,
                   help="seconds before a step is killed (default 86400)")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p.parse_args(argv)


def run(args) -> None:
    paths = {}
    for key in ("table", "rep_seqs", "taxonomy"):
        value = os.path.abspath(os.path.expanduser(getattr(args, key)))
        if not os.path.isfile(value) or os.path.getsize(value) == 0:
            raise FilterError(f"{key.replace('_', '-')} not found or empty: {value}")
        paths[key] = value
    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    if not 0 <= args.min_samples_fraction < 1:
        raise FilterError("--min-samples-fraction must be at least 0 and below 1, got "
                          f"{args.min_samples_fraction}")
    if args.min_sample_reads < 0:
        raise FilterError(f"--min-sample-reads cannot be negative, got "
                          f"{args.min_sample_reads}")
    if args.drop_where and not args.metadata:
        raise FilterError("--drop-where needs --metadata")
    metadata = None
    if args.metadata:
        metadata = os.path.abspath(os.path.expanduser(args.metadata))
        if not os.path.isfile(metadata) or os.path.getsize(metadata) == 0:
            raise FilterError(f"metadata not found or empty: {metadata}")
    os.makedirs(outdir, exist_ok=True)
    setup_logging(outdir)
    log.info("=" * 70)
    log.info("filter %s | Python %s | %s", __version__, platform.python_version(),
             platform.platform())
    log.info("command: %s", shlex.join([sys.executable] + sys.argv))
    log.info("started: %s", datetime.now(timezone.utc).isoformat(timespec="seconds"))

    per_sample, features = export_counts(paths["table"], outdir, "00_input", args.env,
                                         args.timeout)
    start = dict(per_sample)
    steps = [Step("input", "as given", len(per_sample), features, sum(per_sample.values()))]
    log.info("input: %d samples, %d features, %d reads", steps[0].samples,
             steps[0].features, int(steps[0].reads))
    current = paths["table"]

    # 1. taxonomy
    taxa_out = os.path.join(outdir, "table_taxa.qza")
    cmd = ["conda", "run", "-n", args.env, "qiime", "taxa", "filter-table",
           "--i-table", current, "--i-taxonomy", paths["taxonomy"],
           "--p-mode", "contains", "--p-query-delimiter", ",",
           "--o-filtered-table", taxa_out]
    if args.include:
        cmd += ["--p-include", args.include]
    if args.exclude:
        cmd += ["--p-exclude", args.exclude]
    rc, _, err = run_cmd(cmd, args.timeout)
    if rc != 0:
        raise FilterError(f"taxa filter-table exited with code {rc}:\n{_tail(err)}")
    per_sample, features = export_counts(taxa_out, outdir, "01_taxa", args.env,
                                         args.timeout)
    steps.append(Step("taxonomy", f"include {args.include!r}, exclude {args.exclude!r}",
                      len(per_sample), features, sum(per_sample.values())))
    current = taxa_out

    # 2. samples
    if args.drop_where:
        samples_out = os.path.join(outdir, "table_samples.qza")
        rc, _, err = run_cmd(["conda", "run", "-n", args.env, "qiime", "feature-table",
                              "filter-samples", "--i-table", current,
                              "--m-metadata-file", metadata,
                              "--p-where", args.drop_where, "--p-exclude-ids",
                              "--o-filtered-table", samples_out], args.timeout)
        if rc != 0:
            raise FilterError(f"filter-samples exited with code {rc}:\n{_tail(err)}")
        dropped_names = sorted(set(per_sample))
        per_sample, features = export_counts(samples_out, outdir, "02_samples",
                                             args.env, args.timeout)
        dropped_names = [s for s in dropped_names if s not in per_sample]
        steps.append(Step("samples", f"dropped where {args.drop_where}",
                          len(per_sample), features, sum(per_sample.values())))
        current = samples_out
        if not dropped_names:
            log.warning("--drop-where %s matched no samples. Nothing was removed, which "
                        "is usually a wrong column or value rather than a clean table",
                        args.drop_where)
        else:
            log.info("dropped %d sample(s): %s%s", len(dropped_names),
                     ", ".join(dropped_names[:8]),
                     " ..." if len(dropped_names) > 8 else "")

    # 3. prevalence
    if args.min_samples_fraction > 0:
        min_samples = max(1, math.ceil(args.min_samples_fraction * len(per_sample)))
        before_samples = set(per_sample)
        prev_out = os.path.join(outdir, "table_prevalence.qza")
        rc, _, err = run_cmd(["conda", "run", "-n", args.env, "qiime", "feature-table",
                              "filter-features", "--i-table", current,
                              "--p-min-samples", str(min_samples),
                              "--o-filtered-table", prev_out], args.timeout)
        if rc != 0:
            raise FilterError(f"filter-features exited with code {rc}:\n{_tail(err)}")
        per_sample, features = export_counts(prev_out, outdir, "03_prevalence",
                                             args.env, args.timeout)
        steps.append(Step("prevalence",
                          f"kept features in >= {min_samples} samples "
                          f"({args.min_samples_fraction:.1%})",
                          len(per_sample), features, sum(per_sample.values())))
        current = prev_out
        # filter-features has --p-filter-empty-samples on by default, so a feature
        # filter can quietly remove samples. Say so rather than absorb it.
        lost = sorted(before_samples - set(per_sample))
        if lost:
            log.warning("%d sample(s) were removed by the *feature* filter because "
                        "nothing was left in them: %s%s. qiime feature-table "
                        "filter-features has --p-filter-empty-samples on by default",
                        len(lost), ", ".join(lost[:8]), " ..." if len(lost) > 8 else "")

    # 4. depth
    if args.min_sample_reads > 0:
        depth_out = os.path.join(outdir, "table_depth.qza")
        rc, _, err = run_cmd(["conda", "run", "-n", args.env, "qiime", "feature-table",
                              "filter-samples", "--i-table", current,
                              "--p-min-frequency", str(args.min_sample_reads),
                              "--o-filtered-table", depth_out], args.timeout)
        if rc != 0:
            raise FilterError(f"filter-samples exited with code {rc}:\n{_tail(err)}")
        before_samples = set(per_sample)
        per_sample, features = export_counts(depth_out, outdir, "04_depth", args.env,
                                             args.timeout)
        steps.append(Step("depth", f"kept samples with >= {args.min_sample_reads} reads",
                          len(per_sample), features, sum(per_sample.values())))
        current = depth_out
        lost = sorted(before_samples - set(per_sample))
        if lost:
            log.info("dropped %d sample(s) below %d reads: %s%s", len(lost),
                     args.min_sample_reads, ", ".join(lost[:8]),
                     " ..." if len(lost) > 8 else "")

    if not per_sample or not features:
        raise FilterError(
            f"filtering emptied the table: {len(per_sample)} samples, {features} "
            "features left. Check --include against the taxonomy you actually have; a "
            "domain label that does not appear removes everything")

    final_table = os.path.join(outdir, "table_filtered.qza")
    if os.path.abspath(current) != os.path.abspath(final_table):
        rc, _, err = run_cmd(["cp", current, final_table], args.timeout)
        if rc != 0:
            raise FilterError(f"could not write {final_table}:\n{_tail(err)}")

    # The sequences must match the table, or the tree will not match either.
    final_seqs = os.path.join(outdir, "rep_seqs_filtered.qza")
    rc, _, err = run_cmd(["conda", "run", "-n", args.env, "qiime", "feature-table",
                          "filter-seqs", "--i-data", paths["rep_seqs"],
                          "--i-table", final_table,
                          "--o-filtered-data", final_seqs], args.timeout)
    if rc != 0:
        raise FilterError(f"filter-seqs exited with code {rc}:\n{_tail(err)}")
    if not os.path.isfile(final_seqs):
        raise FilterError(f"filter-seqs finished but wrote no {final_seqs}")

    write_summary(os.path.join(outdir, "filter_summary.tsv"), steps, {
        "include": args.include, "exclude": args.exclude,
        "drop_where": args.drop_where or "(none)",
        "min_samples_fraction": args.min_samples_fraction,
        "min_sample_reads": args.min_sample_reads,
    })
    write_per_sample(os.path.join(outdir, "filter_per_sample.tsv"), start, per_sample)

    first, last = steps[0], steps[-1]
    for step in steps[1:]:
        log.info("%-10s %-44s -> %d samples, %d features, %d reads", step.name,
                 step.detail, step.samples, step.features, int(step.reads))
    log.info("kept %d of %d samples, %d of %d features, %d of %d reads (%.1f%%)",
             last.samples, first.samples, last.features, first.features,
             int(last.reads), int(first.reads),
             100 * last.reads / first.reads if first.reads else 0.0)
    log.info("thresholds in filter_summary.tsv are choices recorded with the result. "
             "The prevalence filter in particular decides what counts as evidence, so "
             "state it in the methods rather than leaving it to be inferred")
    log.info("done: %s  (%s)", final_table, final_seqs)


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        run(args)
    except FilterError as exc:
        if not log.handlers:
            setup_logging(None)
        log.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
