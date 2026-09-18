#!/usr/bin/env python3
"""Pool the sequencing runs of each sample, and prove no reads were lost doing it.

A feature table from 03_dada2.py has one column per **run**, not per sample. A
resequenced sample has several runs, and until they are pooled every per-sample number
downstream is wrong: the sample is counted twice, its depth is halved, and its diversity
is measured on a fraction of its reads. Nothing about that raises an error, which is why
this stage exists and why it checks its own arithmetic.

Give it the table, the run-to-sample map written by 00_fetch_sra.py, and optionally the
study metadata. The script

  1. checks the grouping ids are usable before spending compute on them,
  2. checks every run in the table is in the map and the other way round,
  3. pools with qiime feature-table group --p-mode sum,
  4. **proves the total read count is unchanged**, exactly, not approximately,
  5. joins the study metadata and reports what matched.

**Sample ids are not rewritten silently.** An id with a space in it is not a usable
QIIME 2 sample id, but quietly renaming somebody's samples is worse than stopping. The
stage fails, names the offenders, and `--sanitize-ids` is an explicit opt-in that writes
the before and after map so the change is on the record.

Standard library only. Needs Python 3.8 or later, Linux or WSL, and conda with the
QIIME 2 environment (see workflow/setup_envs.sh).

Example:
    python workflow/07_collapse.py -b ~/research/baxter2016/q2/dada2/table.qza \\
        -r ~/research/baxter2016/raw/run_to_sample.tsv \\
        --metadata ~/research/baxter2016/ref/metadata.tsv \\
        --expect-samples 495 -o ~/research/baxter2016/q2/collapsed
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
from datetime import datetime, timezone

__version__ = "0.1.0"

# QIIME 2 metadata ids: no leading #, no whitespace, not empty, not a reserved word.
RESERVED_IDS = {"id", "sampleid", "sample id", "sample-id", "featureid", "feature id",
                "feature-id", "#sampleid", "#otuid", "sample_name", "row id"}

log = logging.getLogger("collapse")


class CollapseError(RuntimeError):
    """A check failed. The message says which one and why."""


def run_cmd(cmd: list[str], timeout: int) -> tuple[int, str, str]:
    """Run a command with stdin closed and a hard time limit."""
    log.info("run: %s", shlex.join(cmd))
    try:
        proc = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True, start_new_session=True)
    except FileNotFoundError as exc:
        raise CollapseError(f"cannot start {cmd[0]}: {exc}") from exc
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.communicate()
        raise CollapseError(f"timed out after {timeout} s and was killed: "
                            f"{shlex.join(cmd)}")
    return proc.returncode, out, err


def _tail(text: str, n: int = 15) -> str:
    return "\n".join(text.strip().splitlines()[-n:])


def read_tsv(path: str) -> list[dict[str, str]]:
    """Read a TSV, tolerating the CR line endings old study metadata often carries."""
    raw = open(path, "rb").read().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    rows = list(csv.DictReader(raw.decode("utf-8").splitlines(), delimiter="\t"))
    if not rows:
        raise CollapseError(f"{path}: no rows")
    return rows


def read_run_map(path: str, column: str) -> dict[str, str]:
    """Read run id to group id. The first column is the run, as QIIME 2 expects."""
    rows = read_tsv(path)
    id_col = list(rows[0])[0]
    if column not in rows[0]:
        raise CollapseError(f"{path} has no column {column!r}. Columns: "
                            f"{', '.join(list(rows[0])[:12])}")
    out: dict[str, str] = {}
    for row in rows:
        run = (row[id_col] or "").strip()
        if not run or run.startswith("#"):
            continue
        group = (row[column] or "").strip()
        if not group:
            raise CollapseError(f"{path}: run {run} has no value in {column!r}, so it "
                                "cannot be assigned to a sample")
        if run in out:
            raise CollapseError(f"{path}: run {run} appears more than once")
        out[run] = group
    if not out:
        raise CollapseError(f"{path}: no runs")
    return out


def unsafe_ids(ids) -> list[str]:
    """Ids QIIME 2 will not accept as sample ids. Returns them, does not fix them."""
    bad = []
    for value in sorted(set(ids)):
        if (not value or value.startswith("#") or value != value.strip()
                or any(c.isspace() for c in value)
                or value.lower() in RESERVED_IDS):
            bad.append(value)
    return bad


def sanitize(value: str) -> str:
    """Deterministic, reversible-by-inspection: whitespace runs become one underscore."""
    return "_".join(value.split())


def export_table_tsv(qza: str, outdir: str, env: str, timeout: int) -> str:
    dest = os.path.join(outdir, os.path.basename(qza).replace(".qza", "_export"))
    rc, _, err = run_cmd(["conda", "run", "-n", env, "qiime", "tools", "export",
                          "--input-path", qza, "--output-path", dest], timeout)
    if rc != 0:
        raise CollapseError(f"qiime tools export failed on {qza} with code {rc}:"
                            f"\n{_tail(err)}")
    biom = os.path.join(dest, "feature-table.biom")
    if not os.path.isfile(biom):
        raise CollapseError(f"export finished but wrote no {biom}")
    tsv = os.path.join(dest, "feature-table.tsv")
    rc, _, err = run_cmd(["conda", "run", "-n", env, "biom", "convert", "-i", biom,
                          "-o", tsv, "--to-tsv"], timeout)
    if rc != 0:
        raise CollapseError(f"biom convert exited with code {rc}:\n{_tail(err)}")
    if not os.path.isfile(tsv):
        raise CollapseError(f"biom convert finished but wrote no {tsv}")
    return tsv


def read_table_totals(tsv: str) -> dict[str, float]:
    """Total reads per column, which is what the accounting check compares."""
    totals: dict[str, float] = {}
    with open(tsv) as fh:
        header = None
        for line in fh:
            if line.startswith("# Const") or not line.strip():
                continue
            cells = line.rstrip("\n").split("\t")
            if header is None:
                if cells[0].lstrip("#").strip().lower() not in ("otu id", "otu_id",
                                                                "featureid",
                                                                "feature id"):
                    raise CollapseError(f"{tsv}: first column is {cells[0]!r}, expected "
                                        "the feature id column")
                header = cells[1:]
                totals = {name: 0.0 for name in header}
                continue
            if len(cells) != len(header) + 1:
                raise CollapseError(f"{tsv}: row {cells[0]!r} has {len(cells) - 1} "
                                    f"values for {len(header)} columns")
            for name, value in zip(header, cells[1:]):
                try:
                    totals[name] += float(value)
                except ValueError as exc:
                    raise CollapseError(f"{tsv}: {cells[0]} in {name} is not a "
                                        f"number: {value!r}") from exc
    if not totals:
        raise CollapseError(f"{tsv}: no columns")
    return totals


def write_group_metadata(path: str, run_map: dict[str, str], column: str) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample-id", column])
        for run in sorted(run_map):
            w.writerow([run, run_map[run]])


def write_id_map(path: str, changes: dict[str, str]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["original_id", "sanitized_id"])
        for before in sorted(changes):
            w.writerow([before, changes[before]])


def write_runs_per_sample(path: str, run_map: dict[str, str],
                          totals: dict[str, float]) -> None:
    by_sample: dict[str, list[str]] = {}
    for run, sample in run_map.items():
        by_sample.setdefault(sample, []).append(run)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample-id", "n_runs", "runs", "reads"])
        for sample in sorted(by_sample, key=lambda s: (-len(by_sample[s]), s)):
            runs = sorted(by_sample[sample])
            w.writerow([sample, len(runs), ";".join(runs),
                        int(sum(totals.get(r, 0.0) for r in runs))])


def join_metadata(path: str, samples: list[str], id_column: str | None,
                  outdir: str) -> tuple[int, list[str], list[str]]:
    """Attach the study metadata to the pooled samples and report what matched."""
    rows = read_tsv(path)
    columns = list(rows[0])
    key = id_column or columns[0]
    if key not in columns:
        raise CollapseError(f"{path} has no column {key!r}. Columns: "
                            f"{', '.join(columns[:12])}")
    by_id: dict[str, dict[str, str]] = {}
    for row in rows:
        value = (row[key] or "").strip()
        if not value or value.startswith("#"):
            continue
        if value in by_id:
            raise CollapseError(f"{path}: {value} appears more than once in {key!r}")
        by_id[value] = row
    matched = [s for s in samples if s in by_id]
    if not matched:
        raise CollapseError(
            f"not one pooled sample matches {path} on {key!r}. The ids do not line up "
            f"(pooled ids look like {samples[0]!r}, metadata like "
            f"{sorted(by_id)[0]!r}). Grouping on a different column may be what you want")
    unmatched_samples = [s for s in samples if s not in by_id]
    unused_metadata = [m for m in sorted(by_id) if m not in set(samples)]
    other = [c for c in columns if c != key]
    out = os.path.join(outdir, "sample_metadata.tsv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample-id", "in_study_metadata"] + other)
        for sample in samples:
            row = by_id.get(sample)
            w.writerow([sample, "yes" if row else "no"]
                       + [(row[c] if row else "") for c in other])
    return len(matched), unmatched_samples, unused_metadata


def setup_logging(outdir: str | None) -> None:
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(fmt)
    log.addHandler(console)
    if outdir:
        fh = logging.FileHandler(os.path.join(outdir, "collapse_log.txt"), mode="a")
        fh.setFormatter(fmt)
        log.addHandler(fh)


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Pool the sequencing runs of each sample, with read accounting.")
    p.add_argument("-b", "--table", required=True, help="run-level feature table")
    p.add_argument("-r", "--run-map", required=True,
                   help="run to sample map, e.g. run_to_sample.tsv from 00_fetch_sra.py")
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("--group-column", default="sample_title",
                   help="run-map column holding the sample id (default sample_title)")
    p.add_argument("--metadata", help="study metadata to attach to the pooled samples")
    p.add_argument("--metadata-id-column",
                   help="metadata column holding the sample id (default its first)")
    p.add_argument("--expect-samples", type=int,
                   help="fail unless pooling gives exactly this many samples")
    p.add_argument("--sanitize-ids", action="store_true",
                   help="replace whitespace in sample ids with underscores. Writes the "
                        "before and after map, because renaming samples is on the record")
    p.add_argument("--env", default="amplipub-qiime2-2025.7", help="conda env with QIIME 2")
    p.add_argument("--timeout", type=int, default=86400,
                   help="seconds before a step is killed (default 86400)")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p.parse_args(argv)


def run(args) -> None:
    table = os.path.abspath(os.path.expanduser(args.table))
    run_map_path = os.path.abspath(os.path.expanduser(args.run_map))
    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    for path in (table, run_map_path):
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            raise CollapseError(f"not found or empty: {path}")
    if args.expect_samples is not None and args.expect_samples < 1:
        raise CollapseError(f"--expect-samples must be at least 1, got "
                            f"{args.expect_samples}")
    os.makedirs(outdir, exist_ok=True)
    setup_logging(outdir)
    log.info("=" * 70)
    log.info("collapse %s | Python %s | %s", __version__, platform.python_version(),
             platform.platform())
    log.info("command: %s", shlex.join([sys.executable] + sys.argv))
    log.info("started: %s", datetime.now(timezone.utc).isoformat(timespec="seconds"))

    run_map = read_run_map(run_map_path, args.group_column)
    groups = sorted(set(run_map.values()))
    log.info("%d runs map to %d samples on %r", len(run_map), len(groups),
             args.group_column)

    bad = unsafe_ids(groups)
    if bad and not args.sanitize_ids:
        raise CollapseError(
            f"{len(bad)} sample id(s) are not usable as QIIME 2 sample ids, for example "
            f"{bad[:3]}. Whitespace, a leading #, and reserved words are all rejected. "
            "Fix them in the run map, group on a different column, or pass "
            "--sanitize-ids to replace whitespace with underscores, which writes the "
            "before and after map")
    if bad:
        changes = {value: sanitize(value) for value in bad}
        clashes = [b for b, a in changes.items() if a in set(groups) - {b}]
        if clashes:
            raise CollapseError(f"sanitizing would make {clashes[0]!r} collide with an "
                                "existing sample id. Fix the ids by hand")
        run_map = {r: changes.get(s, s) for r, s in run_map.items()}
        groups = sorted(set(run_map.values()))
        write_id_map(os.path.join(outdir, "sanitized_ids.tsv"), changes)
        log.warning("%d sample id(s) were sanitized, see sanitized_ids.tsv: %s",
                    len(changes), ", ".join(f"{b!r} -> {a!r}"
                                            for b, a in list(changes.items())[:5]))

    before_tsv = export_table_tsv(table, outdir, args.env, args.timeout)
    before = read_table_totals(before_tsv)
    missing = sorted(set(before) - set(run_map))
    if missing:
        raise CollapseError(
            f"{len(missing)} run(s) in the table are not in the map, for example "
            f"{missing[0]}. feature-table group needs every id present, and a run with "
            "no sample would be dropped without a word")
    absent = sorted(set(run_map) - set(before))
    if absent:
        log.warning("%d run(s) in the map are not in the table and are ignored: %s%s",
                    len(absent), ", ".join(absent[:5]),
                    " ..." if len(absent) > 5 else "")
        run_map = {r: s for r, s in run_map.items() if r in before}
        groups = sorted(set(run_map.values()))

    group_md = os.path.join(outdir, "grouping.tsv")
    write_group_metadata(group_md, run_map, args.group_column)
    grouped = os.path.join(outdir, "table_by_sample.qza")
    rc, _, err = run_cmd(["conda", "run", "-n", args.env, "qiime", "feature-table",
                          "group", "--i-table", table, "--p-axis", "sample",
                          "--m-metadata-file", group_md,
                          "--m-metadata-column", args.group_column,
                          "--p-mode", "sum", "--o-grouped-table", grouped],
                         args.timeout)
    if rc != 0:
        raise CollapseError(f"feature-table group exited with code {rc}:\n{_tail(err)}")
    if not os.path.isfile(grouped):
        raise CollapseError(f"feature-table group finished but wrote no {grouped}")

    after = read_table_totals(export_table_tsv(grouped, outdir, args.env, args.timeout))
    # The whole point of the stage: pooling moves reads between columns, it must not
    # create or destroy any. Compared as integers, because these are counts.
    total_before = round(sum(before.values()))
    total_after = round(sum(after.values()))
    if total_before != total_after:
        raise CollapseError(
            f"pooling changed the read total: {total_before} before, {total_after} "
            f"after, a difference of {total_after - total_before}. --p-mode sum must "
            "conserve reads, so this is a bug or the wrong mode, not a rounding issue")
    if len(after) != len(groups):
        raise CollapseError(f"pooling gave {len(after)} samples but the map has "
                            f"{len(groups)} distinct ids")
    if args.expect_samples is not None and len(after) != args.expect_samples:
        raise CollapseError(f"expected {args.expect_samples} samples, pooling gave "
                            f"{len(after)}")
    for sample, value in after.items():
        expected = sum(before[r] for r, s in run_map.items() if s == sample)
        if round(value) != round(expected):
            raise CollapseError(f"{sample}: pooled to {round(value)} reads but its runs "
                                f"hold {round(expected)}")

    write_runs_per_sample(os.path.join(outdir, "runs_per_sample.tsv"), run_map, before)
    pooled = sum(1 for s in groups
                 if sum(1 for r in run_map if run_map[r] == s) > 1)
    log.info("%d runs pooled into %d samples; %d sample(s) had more than one run; "
             "%d reads in and %d out, unchanged", len(run_map), len(after), pooled,
             total_before, total_after)
    values = sorted(round(v) for v in after.values())
    log.info("per-sample reads after pooling: min %d, median %d, max %d",
             values[0], values[len(values) // 2], values[-1])

    if args.metadata:
        md = os.path.abspath(os.path.expanduser(args.metadata))
        if not os.path.isfile(md) or os.path.getsize(md) == 0:
            raise CollapseError(f"metadata not found or empty: {md}")
        n, unmatched, unused = join_metadata(md, sorted(after), args.metadata_id_column,
                                            outdir)
        log.info("study metadata: %d of %d pooled samples matched", n, len(after))
        if unmatched:
            log.warning("%d pooled sample(s) have no metadata row and are marked "
                        "in_study_metadata=no: %s%s", len(unmatched),
                        ", ".join(unmatched[:8]), " ..." if len(unmatched) > 8 else "")
        if unused:
            log.warning("%d metadata row(s) match no sequenced sample: %s%s",
                        len(unused), ", ".join(unused[:8]),
                        " ..." if len(unused) > 8 else "")
        log.info("controls and mocks normally land in the unmatched list; that is "
                 "expected, and the next stage is where they get removed")
    log.info("done: %s  (%s)", grouped, os.path.join(outdir, "runs_per_sample.tsv"))


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        run(args)
    except CollapseError as exc:
        if not log.handlers:
            setup_logging(None)
        log.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
