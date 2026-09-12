#!/usr/bin/env python3
"""Compare the ASVs recovered from mock community samples with the known reference.

Give it the table and sequences from 03_dada2.py plus the mock reference FASTA. The
script

  1. finds the region between the primers in every reference record and collapses it to
     the set of distinct targets that could actually be recovered,
  2. exports the feature table and the ASV sequences from the QIIME 2 artifacts,
  3. for each mock sample, sorts its ASVs into exact matches and everything else,
  4. writes the numbers, with published figures alongside for comparison.

**Nothing here passes or fails.** No published source sets an acceptance threshold for a
mock community (six were checked on 2026-09-12; see workflow/README.md). Inventing one
and calling it standard would be worse than reporting the numbers and letting the reader
judge, so this stage reports. Published figures are printed next to ours so the
comparison is explicit rather than implied:

  - Kozich et al. 2013 (AEM 79:5112), same lab, same V4 region, same platform: V4 error
    rate 0.25-1.08% raw, 0.05-0.06% after filtering, 0.01% after preclustering. Mock OTUs
    against 20 expected: 22.8-23.5 with perfect chimera removal, 37.2-43.4 with UCHIME.
  - Callahan et al. 2016 (Nat Methods 13:581), HMP mock, merged reads: 40 exact reference
    matches, 2 spurious, 21 of 21 strains.

Those are OTUs in one case and ASVs in the other. They are not the same unit as ours and
the report says so rather than inviting a false comparison.

Three limits stated up front rather than buried.

First, an ASV tens of mismatches from every reference is a different organism, not a
miscalled base. Averaging its distance into an "error rate" measures nothing. This was
the first version's bug, caught on the real mocks on 2026-09-12 when it reported a 15%
"error rate" that was really the distance to unrelated organisms. So the rate is computed
only over ASVs within `--attributable-within` mismatches of a reference, at several
cutoffs so no single one carries the claim, with exact matches in the denominator because
they are reads with zero errors. Everything beyond the widest cutoff is reported as not
attributable to the reference, which is a statement about the sample, not the pipeline.

Second, the comparison is against the nearest reference of the **same length**, counting
differing positions. It is not mothur's `seq.error`, which aligns first, so do not report
it under that name; an ASV with an insertion or deletion has no same-length reference and
is counted as not attributable.

Third, a reference sequence whose primer sites carry mismatches may amplify poorly or not
at all, so a missing target is not automatically an error in the pipeline.

Standard library only. Needs Python 3.8 or later, Linux or WSL, and conda with the
QIIME 2 environment (see workflow/setup_envs.sh).

Example:
    python workflow/04_mock.py -b ~/research/baxter2016/q2/dada2/table.qza \\
        -r ~/research/baxter2016/q2/dada2/rep_seqs.qza \\
        -m ~/research/baxter2016/ref/HMP_MOCK.v35.fasta \\
        -s mock1,mock2,mock5,mock6,mock7 \\
        -o ~/research/baxter2016/q2/mock
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
from datetime import datetime, timezone

__version__ = "0.1.0"

PRIMER_F = "GTGCCAGCMGCCGCGGTAA"      # 515F
PRIMER_R = "GGACTACHVGGGTWTCTAAT"     # 806R, found as its reverse complement

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T",
         "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A",
              "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
              "B": "V", "V": "B", "D": "H", "H": "D", "N": "N"}

log = logging.getLogger("mock")


class MockError(RuntimeError):
    """A check failed. The message says which one and why."""


def run_cmd(cmd: list[str], timeout: int) -> tuple[int, str, str]:
    """Run a command with stdin closed and a hard time limit."""
    log.info("run: %s", shlex.join(cmd))
    try:
        proc = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True, start_new_session=True)
    except FileNotFoundError as exc:
        raise MockError(f"cannot start {cmd[0]}: {exc}") from exc
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.communicate()
        raise MockError(f"timed out after {timeout} s and was killed: {shlex.join(cmd)}")
    return proc.returncode, out, err


def _tail(text: str, n: int = 15) -> str:
    return "\n".join(text.strip().splitlines()[-n:])


# ---- sequence helpers ----------------------------------------------------

def revcomp(seq: str) -> str:
    try:
        return "".join(COMPLEMENT[b] for b in reversed(seq.upper()))
    except KeyError as exc:
        raise MockError(f"cannot complement base {exc.args[0]!r} in {seq!r}") from exc


def read_fasta(path: str) -> dict[str, str]:
    """Read a FASTA file. Duplicate names are an error, not a silent overwrite."""
    records: dict[str, str] = {}
    name = None
    parts: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(parts).upper()
                name = line[1:].split()[0]
                if name in records:
                    raise MockError(f"{path}: {name} appears more than once")
                parts = []
            else:
                if name is None:
                    raise MockError(f"{path}: sequence data before the first header")
                parts.append(line)
    if name is not None:
        records[name] = "".join(parts).upper()
    if not records:
        raise MockError(f"{path}: no sequences")
    return records


def mismatches(primer: str, window: str) -> int:
    """Count positions where the window does not satisfy the primer's IUPAC code."""
    n = 0
    for p, b in zip(primer, window):
        allowed = IUPAC.get(p)
        if allowed is None:
            raise MockError(f"{p!r} is not a IUPAC code in primer {primer!r}")
        if b not in allowed:
            n += 1
    return n


def find_primer(seq: str, primer: str, max_mismatch: int) -> tuple[int, int] | None:
    """Leftmost window matching the primer within max_mismatch. Returns (start, end)."""
    width = len(primer)
    best = None
    for i in range(len(seq) - width + 1):
        m = mismatches(primer, seq[i:i + width])
        if m == 0:
            return i, i + width
        if m <= max_mismatch and best is None:
            best = (i, i + width)
    return best


def extract_region(seq: str, fwd: str, rev: str, max_mismatch: int) -> str | None:
    """The sequence between the forward primer and the reverse primer's complement."""
    f = find_primer(seq, fwd, max_mismatch)
    if f is None:
        return None
    r = find_primer(seq[f[1]:], revcomp(rev), max_mismatch)
    if r is None:
        return None
    return seq[f[1]:f[1] + r[0]]


def build_targets(records: dict[str, str], fwd: str, rev: str,
                  max_mismatch: int) -> tuple[dict[str, list[str]], list[str]]:
    """Map each distinct target sequence to the reference names that produce it.

    Also returns the names where no primer pair was found at this mismatch budget.
    """
    targets: dict[str, list[str]] = {}
    missing = []
    for name, seq in records.items():
        region = extract_region(seq, fwd, rev, max_mismatch)
        if region is None or not region:
            missing.append(name)
            continue
        targets.setdefault(region, []).append(name)
    if not targets:
        raise MockError("no reference record yielded a region between the primers. "
                        "Check the primers and the reference orientation")
    return targets, missing


def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        raise MockError(f"hamming needs equal lengths, got {len(a)} and {len(b)}")
    return sum(1 for x, y in zip(a, b) if x != y)


def nearest_same_length(seq: str, targets: list[str]) -> tuple[int, int] | None:
    """Fewest differing positions against any target of the same length, and its length."""
    same = [t for t in targets if len(t) == len(seq)]
    if not same:
        return None
    best = min(hamming(seq, t) for t in same)
    return best, len(seq)


# ---- QIIME 2 artifacts ---------------------------------------------------

def export_artifacts(table: str, rep_seqs: str, workdir: str, env: str,
                     timeout: int) -> tuple[str, str]:
    """Export the table as TSV and the sequences as FASTA using the QIIME 2 tools."""
    tdir = os.path.join(workdir, "table")
    sdir = os.path.join(workdir, "seqs")
    for src, dest in ((table, tdir), (rep_seqs, sdir)):
        rc, _, err = run_cmd(["conda", "run", "-n", env, "qiime", "tools", "export",
                              "--input-path", src, "--output-path", dest], timeout)
        if rc != 0:
            raise MockError(f"qiime tools export failed on {src} with code {rc}:\n{_tail(err)}")
    biom = os.path.join(tdir, "feature-table.biom")
    fasta = os.path.join(sdir, "dna-sequences.fasta")
    for path in (biom, fasta):
        if not os.path.isfile(path):
            raise MockError(f"export finished but wrote no {path}")
    tsv = os.path.join(tdir, "feature-table.tsv")
    rc, _, err = run_cmd(["conda", "run", "-n", env, "biom", "convert", "-i", biom,
                          "-o", tsv, "--to-tsv"], timeout)
    if rc != 0:
        raise MockError(f"biom convert exited with code {rc}:\n{_tail(err)}")
    if not os.path.isfile(tsv):
        raise MockError(f"biom convert finished but wrote no {tsv}")
    return tsv, fasta


def parse_biom_tsv(text: str, name: str) -> dict[str, dict[str, float]]:
    """Read `biom convert --to-tsv` output into {sample: {feature: count}}."""
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("# Const")]
    if not lines:
        raise MockError(f"{name}: no rows")
    header = lines[0].split("\t")
    if header[0] not in ("#OTU ID", "#OTU_ID", "OTU ID"):
        raise MockError(f"{name}: first column is {header[0]!r}, expected '#OTU ID'")
    samples = header[1:]
    if not samples:
        raise MockError(f"{name}: the table has no samples")
    out: dict[str, dict[str, float]] = {s: {} for s in samples}
    for ln in lines[1:]:
        cells = ln.split("\t")
        if len(cells) != len(header):
            raise MockError(f"{name}: row {cells[0]!r} has {len(cells)} cells for "
                            f"{len(header)} columns")
        for sample, value in zip(samples, cells[1:]):
            try:
                count = float(value)
            except ValueError as exc:
                raise MockError(f"{name}: {cells[0]} in {sample} is not a number: "
                                f"{value!r}") from exc
            if count < 0:
                raise MockError(f"{name}: {cells[0]} in {sample} is negative")
            if count:
                out[sample][cells[0]] = count
    return out


# ---- the report ----------------------------------------------------------

def score_sample(counts: dict[str, float], seqs: dict[str, str],
                 targets: dict[str, list[str]], cutoffs: list[int]) -> dict:
    """Sort one sample's ASVs into exact matches, error variants, and foreign sequences.

    An ASV that differs from every reference by tens of positions is a different
    organism, not a miscalled base. Averaging its distance into an "error rate" measures
    nothing, so the mismatch rate is computed only over ASVs within `cutoffs` mismatches
    of a reference, reported at each cutoff so the number's sensitivity is visible, and
    reads beyond the widest cutoff are reported separately as not attributable to the
    reference.
    """
    target_list = list(targets)
    total = sum(counts.values())
    exact_reads = 0.0
    exact_bases = 0.0
    exact_hits: set[str] = set()
    others = []
    for feature, count in counts.items():
        seq = seqs.get(feature)
        if seq is None:
            raise MockError(f"feature {feature} is in the table but not in the sequences")
        if seq in targets:
            exact_reads += count
            exact_bases += len(seq) * count
            exact_hits.add(seq)
            continue
        near = nearest_same_length(seq, target_list)
        others.append({"feature": feature, "reads": count, "length": len(seq),
                       "mismatches": None if near is None else near[0]})
    other_reads = sum(o["reads"] for o in others)
    comparable = [o for o in others if o["mismatches"] is not None]
    widest = max(cutoffs)
    within = {}
    for k in cutoffs:
        near = [o for o in comparable if o["mismatches"] <= k]
        bases = sum(o["length"] * o["reads"] for o in near)
        errors = sum(o["mismatches"] * o["reads"] for o in near)
        within[k] = {
            "asvs": len(near),
            "reads": sum(o["reads"] for o in near),
            # exact matches are reads with zero mismatches, so they belong in the
            # denominator: the rate is errors per base over everything attributable
            # to the reference, which is what an error rate means
            "mismatch_rate": (errors / (bases + exact_bases)
                              if bases + exact_bases else 0.0),
            "mismatch_rate_variants_only": errors / bases if bases else 0.0,
        }
    foreign = [o for o in comparable if o["mismatches"] > widest]
    unattributable_reads = (sum(o["reads"] for o in foreign)
                            + sum(o["reads"] for o in others if o["mismatches"] is None))
    return {
        "reads": total,
        "asvs": len(counts),
        "targets_recovered": len(exact_hits),
        "targets_total": len(targets),
        "exact_reads": exact_reads,
        "exact_read_fraction": exact_reads / total if total else 0.0,
        "other_asvs": len(others),
        "other_reads": other_reads,
        "other_read_fraction": other_reads / total if total else 0.0,
        "asvs_without_same_length_reference": len(others) - len(comparable),
        "within": within,
        "cutoffs": list(cutoffs),
        "unattributable_asvs": len(foreign) + (len(others) - len(comparable)),
        "unattributable_reads": unattributable_reads,
        "unattributable_read_fraction": unattributable_reads / total if total else 0.0,
        "recovered": exact_hits,
        "others": others,
    }


def write_summary(path: str, rows: list[tuple[str, dict]], cutoffs: list[int]) -> None:
    base = ["reads", "asvs", "targets_recovered", "targets_total", "exact_read_fraction",
            "other_asvs", "other_read_fraction"]
    per_cutoff = []
    for k in cutoffs:
        per_cutoff += [f"variant_asvs_within_{k}", f"variant_read_fraction_within_{k}",
                       f"mismatch_rate_within_{k}"]
    tail = ["unattributable_asvs", "unattributable_read_fraction",
            "asvs_without_same_length_reference"]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample-id"] + base + per_cutoff + tail)
        for sample, r in rows:
            cells = [round(r[c], 6) if isinstance(r[c], float) else r[c] for c in base]
            for k in cutoffs:
                v = r["within"][k]
                cells += [v["asvs"], round(v["reads"] / r["reads"] if r["reads"] else 0.0, 6),
                          round(v["mismatch_rate"], 8)]
            cells += [r[c] if not isinstance(r[c], float) else round(r[c], 6) for c in tail]
            w.writerow([sample] + cells)


def write_missing(path: str, rows: list[tuple[str, dict]],
                  targets: dict[str, list[str]]) -> None:
    """Which reference targets were not recovered exactly, by sample."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample-id", "reference_names", "length"])
        for sample, r in rows:
            for seq, names in targets.items():
                if seq not in r["recovered"]:
                    w.writerow([sample, ";".join(sorted(names)), len(seq)])


def write_others(path: str, rows: list[tuple[str, dict]]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample-id", "feature", "reads", "length", "mismatches_to_nearest"])
        for sample, r in rows:
            for o in sorted(r["others"], key=lambda x: -x["reads"]):
                w.writerow([sample, o["feature"], o["reads"], o["length"],
                            "" if o["mismatches"] is None else o["mismatches"]])


def write_targets(path: str, targets: dict[str, list[str]], missing: list[str]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["reference_names", "length", "sequence"])
        for seq, names in sorted(targets.items(), key=lambda kv: sorted(kv[1])[0]):
            w.writerow([";".join(sorted(names)), len(seq), seq])
        for name in missing:
            w.writerow([name, "", "NO REGION FOUND BETWEEN THE PRIMERS"])


def setup_logging(outdir: str | None) -> None:
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(fmt)
    log.addHandler(console)
    if outdir:
        fh = logging.FileHandler(os.path.join(outdir, "mock_log.txt"), mode="a")
        fh.setFormatter(fmt)
        log.addHandler(fh)


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Compare mock community ASVs with a reference of known composition. "
                    "Reports; does not pass or fail.")
    p.add_argument("-b", "--table", required=True, help="feature table, e.g. table.qza")
    p.add_argument("-r", "--rep-seqs", required=True, help="ASV sequences, e.g. rep_seqs.qza")
    p.add_argument("-m", "--mock-reference", required=True, help="reference FASTA")
    p.add_argument("-s", "--mock-samples", required=True,
                   help="comma-separated sample ids, or a file with one per line")
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("--primer-f", default=PRIMER_F, help=f"forward primer (default {PRIMER_F}, 515F)")
    p.add_argument("--primer-r", default=PRIMER_R, help=f"reverse primer (default {PRIMER_R}, 806R)")
    p.add_argument("--max-primer-mismatch", type=int, default=1,
                   help="mismatches allowed per primer site in the reference (default 1). "
                        "Exact sites are always preferred; this only rescues records that "
                        "would otherwise be dropped")
    p.add_argument("--attributable-within", default="1,3,10",
                   help="mismatch cutoffs at which an ASV still counts as an error "
                        "variant of a reference rather than a different organism "
                        "(default 1,3,10). The mismatch rate is reported at each, so no "
                        "single cutoff carries the claim; reads beyond the widest are "
                        "reported as not attributable to the reference")
    p.add_argument("--low-depth-note", type=int, default=0,
                   help="log a note for any mock sample below this many reads (default 0, off)")
    p.add_argument("--env", default="qiime2-amplicon-2025.7", help="conda env with QIIME 2")
    p.add_argument("--timeout", type=int, default=3600,
                   help="seconds before an export step is killed (default 3600)")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p.parse_args(argv)


def read_sample_list(value: str) -> list[str]:
    if os.path.isfile(value):
        with open(value) as fh:
            names = [ln.split("#")[0].strip() for ln in fh]
    else:
        names = [v.strip() for v in value.split(",")]
    names = [n for n in names if n]
    if not names:
        raise MockError(f"no mock sample ids in {value!r}")
    seen = set()
    out = []
    for n in names:
        if n not in seen:
            seen.add(n)
            out.append(n)
    return out


def read_cutoffs(value: str) -> list[int]:
    """Parse --attributable-within into a sorted list of distinct non-negative cutoffs."""
    try:
        values = sorted({int(v) for v in value.split(",") if v.strip()})
    except ValueError as exc:
        raise MockError(f"--attributable-within must be whole numbers: {value!r}") from exc
    if not values:
        raise MockError("--attributable-within is empty")
    if values[0] < 1:
        raise MockError(f"--attributable-within must be 1 or more, got {values[0]}")
    return values


def run(args) -> None:
    table = os.path.abspath(os.path.expanduser(args.table))
    rep_seqs = os.path.abspath(os.path.expanduser(args.rep_seqs))
    reference = os.path.abspath(os.path.expanduser(args.mock_reference))
    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    for path in (table, rep_seqs, reference):
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            raise MockError(f"not found or empty: {path}")
    if args.max_primer_mismatch < 0:
        raise MockError(f"--max-primer-mismatch cannot be negative, got {args.max_primer_mismatch}")
    cutoffs = read_cutoffs(args.attributable_within)
    wanted = read_sample_list(args.mock_samples)
    os.makedirs(outdir, exist_ok=True)
    setup_logging(outdir)
    log.info("=" * 70)
    log.info("mock %s | Python %s | %s", __version__, platform.python_version(),
             platform.platform())
    log.info("command: %s", shlex.join([sys.executable] + sys.argv))
    log.info("started: %s", datetime.now(timezone.utc).isoformat(timespec="seconds"))
    log.info("this stage reports and does not pass or fail: no published source sets a "
             "mock community acceptance threshold")

    records = read_fasta(reference)
    exact_targets, exact_missing = build_targets(records, args.primer_f, args.primer_r, 0)
    targets, missing = build_targets(records, args.primer_f, args.primer_r,
                                     args.max_primer_mismatch)
    log.info("reference: %d records, %d distinct targets with exact primer sites, "
             "%d with up to %d mismatch(es) allowed",
             len(records), len(exact_targets), len(targets), args.max_primer_mismatch)
    rescued = sorted(set(exact_missing) - set(missing))
    if rescued:
        log.info("only found with a relaxed primer site: %s. A mismatched primer site can "
                 "mean the template amplifies poorly, so treat a missing target with care",
                 ", ".join(rescued))
    if missing:
        log.warning("no region between the primers in: %s", ", ".join(sorted(missing)))
    lengths = sorted({len(t) for t in targets})
    log.info("target lengths: %s", ", ".join(str(x) for x in lengths))
    write_targets(os.path.join(outdir, "targets.tsv"), targets, missing)

    tsv, fasta = export_artifacts(table, rep_seqs, outdir, args.env, args.timeout)
    with open(tsv) as fh:
        counts = parse_biom_tsv(fh.read(), tsv)
    seqs = read_fasta(fasta)
    absent = [s for s in wanted if s not in counts]
    if absent:
        raise MockError(f"mock sample(s) not in the table: {', '.join(absent)}. "
                        f"The table has {len(counts)} samples, for example "
                        f"{', '.join(sorted(counts)[:3])}")

    rows = [(s, score_sample(counts[s], seqs, targets, cutoffs)) for s in wanted]
    write_summary(os.path.join(outdir, "mock_summary.tsv"), rows, cutoffs)
    write_missing(os.path.join(outdir, "mock_missing_targets.tsv"), rows, targets)
    write_others(os.path.join(outdir, "mock_other_asvs.tsv"), rows)

    for sample, r in rows:
        log.info("%s: %d reads, %d ASVs, %d of %d targets recovered exactly "
                 "(%.2f%% of reads)", sample, r["reads"], r["asvs"],
                 r["targets_recovered"], r["targets_total"],
                 100 * r["exact_read_fraction"])
        for k in cutoffs:
            v = r["within"][k]
            log.info("    within %d mismatch(es): %d ASV(s), %.2f%% of reads, "
                     "error rate %.4f%% over the reads attributable to the reference",
                     k, v["asvs"], 100 * v["reads"] / r["reads"] if r["reads"] else 0.0,
                     100 * v["mismatch_rate"])
        log.info("    not attributable to the reference: %d ASV(s), %.2f%% of reads. "
                 "Beyond %d mismatches an ASV is a different organism, not a miscalled "
                 "base, so it is excluded from the error rate rather than averaged into it",
                 r["unattributable_asvs"], 100 * r["unattributable_read_fraction"],
                 max(cutoffs))
        if r["asvs_without_same_length_reference"]:
            log.info("    %d ASV(s) have no reference of the same length (an indel), "
                     "counted as not attributable",
                     r["asvs_without_same_length_reference"])
        if args.low_depth_note and r["reads"] < args.low_depth_note:
            log.warning("  %s has %d reads, below the %d you asked to be told about. "
                        "Read these numbers as a description of this library, not as "
                        "evidence about the pipeline", sample, r["reads"],
                        args.low_depth_note)
    log.info("for comparison, not as a threshold: Kozich et al. 2013 report a V4 error "
             "rate of 0.01%% after preclustering and 37.2-43.4 mock OTUs against 20 "
             "expected with UCHIME; Callahan et al. 2016 report 40 exact matches and 2 "
             "spurious on the HMP mock. Those are OTUs and ASVs from other datasets, not "
             "the same unit as these numbers")
    log.info("done: %s", os.path.join(outdir, "mock_summary.tsv"))


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        run(args)
    except MockError as exc:
        if not log.handlers:
            setup_logging(None)
        log.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
