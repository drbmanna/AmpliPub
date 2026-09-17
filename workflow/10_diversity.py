#!/usr/bin/env python3
"""Alpha and beta diversity, with the rarefaction depth chosen from evidence.

The depth decision is the one people agonise over and then settle with a round number.
There is no published rule. What there is, verified 2026-09-12, is Schloss 2024 (mSphere,
PMC10900887), which found rarefaction the only approach that controls uneven sequencing
effort across common alpha and beta metrics, over datasets spanning 100-fold variation,
and which says of choosing the threshold:

    "My personal process for selecting a rarefaction threshold involves looking for a
    natural break in the distribution of the number of sequences."

So this stage computes that break, reports how many samples each candidate depth costs,
and reports **Good's coverage** so adequacy is measured rather than assumed. With
`--auto-depth` it picks the break and says which one it picked. Without a depth it refuses
to guess.

**Two things the standard tool does that are worth knowing.**

`core-metrics-phylogenetic` subsamples **once**. Schloss separates *rarefying*, a single
subsample, from *rarefaction*, repeating it 100 to 1,000 times and averaging, and argues
the conflation "was lost on many subsequent researchers". The single subsample is what
QIIME 2 gives you, so this stage also runs `alpha-rarefaction`, which does repeat
(`--p-iterations`), and the log says plainly which number came from which. Do not present
a single-subsample metric as a rarefied one.

**Good's coverage must be computed before prevalence filtering.** It is built on
singletons, and a prevalence filter removes exactly those, which inflates coverage toward
1.0 and makes the statistic meaningless. Pass the pre-filter table as
`--coverage-table`; if you do not, the stage says the number is not trustworthy rather
than printing it as if it were.

Standard library only. Needs Python 3.8 or later, Linux or WSL, and conda with the
QIIME 2 environment (see workflow/setup_envs.sh).

Example:
    python workflow/10_diversity.py -b ~/research/baxter2016/q2/filtered/table_filtered.qza \\
        -p ~/research/baxter2016/q2/tree/rooted_tree.qza \\
        -m ~/research/baxter2016/q2/collapsed/sample_metadata.tsv \\
        --coverage-table ~/research/baxter2016/q2/collapsed/table_by_sample.qza \\
        --depth 10000 --group-column dx -o ~/research/baxter2016/q2/diversity
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

ITERATIONS = 100      # Schloss 2024 says 100 to 1,000; the low end of his range
MAX_SAMPLE_LOSS = 0.1  # our choice: refuse to silently discard a tenth of the study
BREAK_WINDOW = 0.25   # look for the break in the lowest quarter of the distribution

log = logging.getLogger("diversity")


class DiversityError(RuntimeError):
    """A check failed. The message says which one and why."""


def run_cmd(cmd: list[str], timeout: int) -> tuple[int, str, str]:
    """Run a command with stdin closed and a hard time limit."""
    log.info("run: %s", shlex.join(cmd))
    try:
        proc = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True, start_new_session=True)
    except FileNotFoundError as exc:
        raise DiversityError(f"cannot start {cmd[0]}: {exc}") from exc
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.communicate()
        raise DiversityError(f"timed out after {timeout} s and was killed: "
                             f"{shlex.join(cmd)}")
    return proc.returncode, out, err


def _excerpt(text: str, head: int = 5, tail: int = 15) -> str:
    """Both ends of a failed command's output, not just the tail.

    QIIME 2 reports an argument error as a numbered list whose first entry is the reason
    and whose remaining entries are consequences of it, so a tail-only excerpt drops the
    one line that explains the failure.
    """
    lines = text.strip().splitlines()
    if len(lines) <= head + tail:
        return "\n".join(lines)
    omitted = len(lines) - head - tail
    return "\n".join(lines[:head] + [f"... {omitted} line(s) omitted ..."] + lines[-tail:])


def fresh_output_dir(path: str) -> None:
    """Clear the way for `qiime ... --output-dir`, which refuses an existing directory.

    A workflow engine creates the parent directories of a rule's declared outputs before
    the rule runs, so this stage can be handed an empty `core_metrics/` that it did not
    make. That empty directory is removed. A non-empty one holds real results and stops
    the stage rather than being deleted.
    """
    if not os.path.isdir(path):
        return
    contents = os.listdir(path)
    if contents:
        raise DiversityError(
            f"{path} already exists and is not empty ({len(contents)} entries). "
            "core-metrics-phylogenetic will not write into it. Move or delete it, or "
            "run into a fresh output directory")
    os.rmdir(path)


def read_from_artifact(qza: str, suffix: str) -> str:
    try:
        with zipfile.ZipFile(qza) as zf:
            hits = [n for n in zf.namelist() if n.endswith(suffix)]
            if len(hits) != 1:
                raise DiversityError(f"{qza}: expected one {suffix}, found {len(hits)}")
            return zf.read(hits[0]).decode("utf-8")
    except zipfile.BadZipFile as exc:
        raise DiversityError(f"{qza} is not a readable artifact: {exc}") from exc


def export_table(qza: str, outdir: str, tag: str, env: str, timeout: int) -> str:
    dest = os.path.join(outdir, f"export_{tag}")
    rc, _, err = run_cmd(["conda", "run", "-n", env, "qiime", "tools", "export",
                          "--input-path", qza, "--output-path", dest], timeout)
    if rc != 0:
        raise DiversityError(f"export failed on {qza} with code {rc}:\n{_excerpt(err)}")
    biom = os.path.join(dest, "feature-table.biom")
    if not os.path.isfile(biom):
        raise DiversityError(f"export finished but wrote no {biom}")
    tsv = os.path.join(dest, "feature-table.tsv")
    rc, _, err = run_cmd(["conda", "run", "-n", env, "biom", "convert", "-i", biom,
                          "-o", tsv, "--to-tsv"], timeout)
    if rc != 0:
        raise DiversityError(f"biom convert exited with code {rc}:\n{_excerpt(err)}")
    return tsv


def read_columns(tsv: str) -> dict[str, list[float]]:
    """Per-sample feature counts, which is what depth and coverage both need."""
    columns: dict[str, list[float]] = {}
    with open(tsv) as fh:
        header = None
        for line in fh:
            if line.startswith("# Const") or not line.strip():
                continue
            cells = line.rstrip("\n").split("\t")
            if header is None:
                header = cells[1:]
                columns = {name: [] for name in header}
                continue
            if len(cells) != len(header) + 1:
                raise DiversityError(f"{tsv}: row {cells[0]!r} has {len(cells) - 1} "
                                     f"values for {len(header)} columns")
            for name, value in zip(header, cells[1:]):
                try:
                    columns[name].append(float(value))
                except ValueError as exc:
                    raise DiversityError(f"{tsv}: {cells[0]} in {name} is not a "
                                         f"number") from exc
    if not columns:
        raise DiversityError(f"{tsv}: no samples")
    return columns


def depths(columns: dict[str, list[float]]) -> dict[str, int]:
    return {name: int(round(sum(values))) for name, values in columns.items()}


def count_singletons(columns: dict[str, list[float]]) -> int:
    return sum(1 for values in columns.values() for v in values if round(v) == 1)


def goods_coverage(columns: dict[str, list[float]]) -> dict[str, float]:
    """1 - singletons/reads. Built on singletons, so degenerate wherever they are absent.

    Two ways it becomes meaningless, both of which look like a perfect score:
    a prevalence filter removes singletons, and so does DADA2 by design. The DADA2
    maintainer, github.com/benjjneb/dada2 issue 1491, fetched 2026-09-12: "it does not
    infer *amplicon sequence variants* that are only supported by a single read -
    singletons are assumed too difficult to differentiate from errors. Hence no
    singletons in the output table of amplicon sequence variants." So on any ASV table
    from DADA2 this returns 1.0 for every sample, which is a fact about the denoiser and
    not evidence that anything was sequenced deeply enough. The caller checks for that.
    """
    out = {}
    for name, values in columns.items():
        total = sum(values)
        singletons = sum(1 for v in values if round(v) == 1)
        out[name] = (1 - singletons / total) if total else 0.0
    return out


def natural_break(values: list[int], window: float = BREAK_WINDOW) -> dict | None:
    """The largest relative gap in the low tail of the depth distribution.

    Schloss 2024's stated method, made countable: sort the depths and find where the
    jump from one sample to the next is largest. Only the low tail is searched, because
    a gap between two deep samples says nothing about where to cut.
    """
    ordered = sorted(values)
    if len(ordered) < 4:
        return None
    limit = max(3, int(len(ordered) * window))
    best = None
    for i in range(1, limit):
        low, high = ordered[i - 1], ordered[i]
        if low <= 0:
            continue
        ratio = high / low
        if best is None or ratio > best["ratio"]:
            best = {"ratio": ratio, "below": i, "low": low, "high": high}
    if best is None or best["ratio"] <= 1.0:
        return None
    return best


def write_depths(path: str, depth_by_sample: dict[str, int],
                 coverage: dict[str, float] | None, chosen: int) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample-id", "reads", "kept_at_depth", "goods_coverage"])
        for sample in sorted(depth_by_sample, key=lambda s: depth_by_sample[s]):
            w.writerow([sample, depth_by_sample[sample],
                        "yes" if depth_by_sample[sample] >= chosen else "no",
                        "" if coverage is None
                        else round(coverage.get(sample, 0.0), 6)])


def write_candidates(path: str, values: list[int], candidates: list[int]) -> None:
    total = len(values)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["depth", "samples_kept", "samples_lost", "fraction_kept",
                    "reads_used"])
        for depth in candidates:
            kept = sum(1 for v in values if v >= depth)
            w.writerow([depth, kept, total - kept,
                        round(kept / total, 6) if total else 0.0, depth * kept])


def setup_logging(outdir: str | None) -> None:
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(fmt)
    log.addHandler(console)
    if outdir:
        fh = logging.FileHandler(os.path.join(outdir, "diversity_log.txt"), mode="a")
        fh.setFormatter(fmt)
        log.addHandler(fh)


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Alpha and beta diversity, with the rarefaction depth chosen from "
                    "evidence rather than habit.")
    p.add_argument("-b", "--table", required=True, help="filtered feature table")
    p.add_argument("-p", "--phylogeny", required=True, help="rooted tree")
    p.add_argument("-m", "--metadata", required=True, help="sample metadata")
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("--depth", type=int,
                   help="rarefaction depth. Without it, nothing is guessed: either "
                        "pass one or pass --auto-depth")
    p.add_argument("--auto-depth", action="store_true",
                   help="use the natural break in the depth distribution, Schloss "
                        "2024's stated method, and log which depth that gave")
    p.add_argument("--coverage-table",
                   help="table from BEFORE prevalence filtering, for Good's coverage. "
                        "Coverage is built on singletons, which a prevalence filter "
                        "removes")
    p.add_argument("--group-column", action="append", default=[],
                   help="metadata column to test groups on. Repeat for more than one")
    p.add_argument("--iterations", type=int, default=ITERATIONS,
                   help=f"rarefaction curve iterations (default {ITERATIONS}; Schloss "
                        "2024 says 100 to 1,000)")
    p.add_argument("--max-sample-loss", type=float, default=MAX_SAMPLE_LOSS,
                   help="refuse a depth that drops more than this fraction of samples "
                        f"(default {MAX_SAMPLE_LOSS}, our choice)")
    p.add_argument("--threads", type=int, default=1,
                   help="threads for core-metrics (default 1, as QIIME 2 does)")
    p.add_argument("--env", default="qiime2-amplicon-2025.7", help="conda env with QIIME 2")
    p.add_argument("--timeout", type=int, default=86400,
                   help="seconds before a step is killed (default 86400)")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p.parse_args(argv)


def run(args) -> None:
    table = os.path.abspath(os.path.expanduser(args.table))
    tree = os.path.abspath(os.path.expanduser(args.phylogeny))
    metadata = os.path.abspath(os.path.expanduser(args.metadata))
    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    for label, path in (("table", table), ("phylogeny", tree), ("metadata", metadata)):
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            raise DiversityError(f"{label} not found or empty: {path}")
    if args.depth is None and not args.auto_depth:
        raise DiversityError(
            "no depth given. There is no standard rarefaction depth and this stage will "
            "not invent one: pass --depth, or --auto-depth to use the natural break in "
            "your own distribution (Schloss 2024's stated method)")
    if args.depth is not None and args.depth < 1:
        raise DiversityError(f"--depth must be at least 1, got {args.depth}")
    if args.iterations < 1:
        raise DiversityError(f"--iterations must be at least 1, got {args.iterations}")
    if not 0 <= args.max_sample_loss < 1:
        raise DiversityError("--max-sample-loss must be at least 0 and below 1, got "
                             f"{args.max_sample_loss}")
    os.makedirs(outdir, exist_ok=True)
    setup_logging(outdir)
    log.info("=" * 70)
    log.info("diversity %s | Python %s | %s", __version__, platform.python_version(),
             platform.platform())
    log.info("command: %s", shlex.join([sys.executable] + sys.argv))
    log.info("started: %s", datetime.now(timezone.utc).isoformat(timespec="seconds"))

    columns = read_columns(export_table(table, outdir, "table", args.env, args.timeout))
    depth_by_sample = depths(columns)
    values = sorted(depth_by_sample.values())
    fold = values[-1] // values[0] if values[0] else 0
    log.info("%d samples; depth min %d, median %d, max %d, %d-fold range",
             len(values), values[0], values[len(values) // 2], values[-1], fold)

    break_at = natural_break(values)
    if break_at:
        log.info("natural break in the low tail: %d -> %d (x%.3f), %d sample(s) below it",
                 break_at["low"], break_at["high"], break_at["ratio"],
                 break_at["below"])
    else:
        log.info("no natural break stands out in the low tail. That is informative: the "
                 "depth choice is arbitrary here, so report how sensitive your "
                 "conclusions are to it")

    candidates = sorted({1000, 2000, 5000, 10000, 20000}
                        | ({break_at["high"]} if break_at else set())
                        | ({args.depth} if args.depth else set()))
    write_candidates(os.path.join(outdir, "depth_candidates.tsv"), values, candidates)

    if args.depth is not None:
        depth = args.depth
        source = "given on the command line"
    elif break_at:
        depth = break_at["high"]
        source = (f"the natural break ({break_at['low']} -> {break_at['high']}, "
                  f"x{break_at['ratio']:.3f})")
    else:
        raise DiversityError(
            "--auto-depth was asked for but no natural break stands out in this "
            "distribution. Choose a depth yourself from depth_candidates.tsv and say "
            "in the methods why, since the data does not point anywhere")
    kept = sum(1 for v in values if v >= depth)
    lost = len(values) - kept
    log.info("depth %d, %s: keeps %d of %d samples, drops %d (%.1f%%)", depth, source,
             kept, len(values), lost, 100 * lost / len(values))
    if depth > values[-1]:
        raise DiversityError(f"depth {depth} is above every sample's read count "
                             f"(max {values[-1]}); nothing would be left")
    if lost / len(values) > args.max_sample_loss:
        raise DiversityError(
            f"depth {depth} drops {lost} of {len(values)} samples "
            f"({100 * lost / len(values):.1f}%), past the "
            f"{100 * args.max_sample_loss:.0f}% you allowed. Lower the depth, or raise "
            "--max-sample-loss deliberately and record why")

    coverage = None
    if args.coverage_table:
        cov_table = os.path.abspath(os.path.expanduser(args.coverage_table))
        if not os.path.isfile(cov_table) or os.path.getsize(cov_table) == 0:
            raise DiversityError(f"coverage table not found or empty: {cov_table}")
        cov_columns = read_columns(export_table(cov_table, outdir, "coverage", args.env,
                                                args.timeout))
        coverage = goods_coverage(cov_columns)
        shared = [s for s in depth_by_sample if s in coverage]
        if not shared:
            raise DiversityError("no sample in the coverage table matches the analysis "
                                 "table; are they from the same run?")
        singletons = count_singletons(cov_columns)
        got = sorted(coverage[s] for s in shared)
        log.info("Good's coverage over %d shared sample(s): min %.4f, median %.4f, "
                 "max %.4f", len(shared), got[0], got[len(got) // 2], got[-1])
        if singletons == 0:
            coverage_note = "degenerate: no singletons in the table"
            log.warning(
                "there is not one singleton in this table, so Good's coverage is 1.0 "
                "everywhere by construction and is NOT evidence that anything was "
                "sequenced deeply enough. DADA2 does not emit singleton ASVs by design "
                "(benjjneb/dada2 issue 1491: singletons are \"too difficult to "
                "differentiate from errors\"), so this statistic cannot work on an ASV "
                "table from it. Good's coverage belongs to OTU pipelines that keep "
                "singletons. Do not report it here")
            log.warning("the same cause makes rarefaction curves plateau early on DADA2 "
                        "output (benjjneb/dada2 issue 317), so read the curves with that "
                        "in mind rather than as evidence of saturation")
        else:
            coverage_note = f"{singletons} singleton(s) in the table"
            log.info("computed on the pre-filter table, which is the only place it can "
                     "mean anything: a prevalence filter removes the singletons it is "
                     "built on")
    else:
        coverage_note = "(not computed)"
        log.warning("no --coverage-table, so Good's coverage is not reported. Computing "
                    "it on a prevalence-filtered table would push it toward 1.0 and say "
                    "nothing, so it is left out rather than printed as if it were real")

    write_depths(os.path.join(outdir, "sample_depths.tsv"), depth_by_sample, coverage,
                 depth)

    # Repeated subsampling, which is what Schloss calls rarefaction.
    curves = os.path.join(outdir, "alpha_rarefaction.qzv")
    cmd = ["conda", "run", "-n", args.env, "qiime", "diversity", "alpha-rarefaction",
           "--i-table", table, "--i-phylogeny", tree,
           "--p-max-depth", str(depth), "--p-min-depth", "1",
           "--p-steps", "10", "--p-iterations", str(args.iterations),
           "--m-metadata-file", metadata, "--o-visualization", curves]
    rc, _, err = run_cmd(cmd, args.timeout)
    if rc != 0:
        raise DiversityError(f"alpha-rarefaction exited with code {rc}:\n{_excerpt(err)}")
    log.info("rarefaction curves: %d iterations per step, averaged. This is rarefaction "
             "in Schloss 2024's sense", args.iterations)

    core = os.path.join(outdir, "core_metrics")
    fresh_output_dir(core)
    rc, _, err = run_cmd(["conda", "run", "-n", args.env, "qiime", "diversity",
                          "core-metrics-phylogenetic", "--i-table", table,
                          "--i-phylogeny", tree, "--p-sampling-depth", str(depth),
                          "--m-metadata-file", metadata,
                          "--p-n-jobs-or-threads", str(args.threads),
                          "--output-dir", core], args.timeout)
    if rc != 0:
        raise DiversityError(f"core-metrics-phylogenetic exited with code {rc}:"
                             f"\n{_excerpt(err)}")
    if not os.path.isdir(core):
        raise DiversityError(f"core-metrics finished but wrote no {core}")
    log.warning("core-metrics-phylogenetic subsamples ONCE at depth %d. Schloss 2024 "
                "separates that (rarefying) from rarefaction, which repeats and "
                "averages. Report these as single-subsample metrics, and use the "
                "rarefaction curves where the distinction matters", depth)

    tests = []
    for column in args.group_column:
        alpha_out = os.path.join(outdir, f"alpha_{column}.qzv")
        for metric in ("faith_pd_vector", "shannon_vector", "observed_features_vector"):
            vector = os.path.join(core, f"{metric}.qza")
            if not os.path.isfile(vector):
                log.warning("%s not found in core-metrics output, skipping", metric)
                continue
            out = os.path.join(outdir, f"alpha_{metric.replace('_vector', '')}_{column}.qzv")
            rc, _, err = run_cmd(["conda", "run", "-n", args.env, "qiime", "diversity",
                                  "alpha-group-significance",
                                  "--i-alpha-diversity", vector,
                                  "--m-metadata-file", metadata,
                                  "--o-visualization", out], args.timeout)
            if rc != 0:
                raise DiversityError(f"alpha-group-significance failed on {metric} "
                                     f"with code {rc}:\n{_excerpt(err)}")
            tests.append(("alpha", metric, column, out))
        for matrix in ("unweighted_unifrac_distance_matrix",
                       "weighted_unifrac_distance_matrix",
                       "bray_curtis_distance_matrix"):
            dm = os.path.join(core, f"{matrix}.qza")
            if not os.path.isfile(dm):
                log.warning("%s not found in core-metrics output, skipping", matrix)
                continue
            out = os.path.join(outdir,
                               f"beta_{matrix.replace('_distance_matrix', '')}_{column}.qzv")
            rc, _, err = run_cmd(["conda", "run", "-n", args.env, "qiime", "diversity",
                                  "beta-group-significance",
                                  "--i-distance-matrix", dm,
                                  "--m-metadata-file", metadata,
                                  "--m-metadata-column", column,
                                  "--p-method", "permanova",
                                  "--p-permutations", "999", "--p-pairwise",
                                  "--o-visualization", out], args.timeout)
            if rc != 0:
                raise DiversityError(f"beta-group-significance failed on {matrix} "
                                     f"with code {rc}:\n{_excerpt(err)}")
            tests.append(("beta", matrix, column, out))
    if not args.group_column:
        log.warning("no --group-column, so no group test was run. The metrics are "
                    "computed but nothing was compared")
    else:
        log.info("%d group test(s) written across %d column(s), PERMANOVA with 999 "
                 "permutations", len(tests), len(args.group_column))

    with open(os.path.join(outdir, "diversity_settings.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["setting", "value", "source"])
        w.writerow(["sampling_depth", depth, source])
        w.writerow(["samples_kept", kept, "at this depth"])
        w.writerow(["samples_dropped", lost, "at this depth"])
        w.writerow(["curve_iterations", args.iterations,
                    "Schloss 2024 says 100 to 1,000"])
        w.writerow(["core_metrics_subsampling", "single",
                    "core-metrics-phylogenetic subsamples once, not repeatedly"])
        w.writerow(["permanova_permutations", 999, "QIIME 2 default"])
        w.writerow(["goods_coverage_source",
                    args.coverage_table or "(not computed)",
                    "must be a pre-prevalence-filter table"])
        w.writerow(["goods_coverage_status", coverage_note,
                    "1.0 everywhere means the table has no singletons, which DADA2 "
                    "never emits; not evidence of depth"])
    log.info("done: %s  (%s)", core, os.path.join(outdir, "diversity_settings.tsv"))


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        run(args)
    except DiversityError as exc:
        if not log.handlers:
            setup_logging(None)
        log.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
