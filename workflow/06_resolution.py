#!/usr/bin/env python3
"""What can this amplicon region actually tell apart, and what does it not.

Give it a reference database and one or more primer pairs. For each region the script

  1. finds what lies between the primers in every reference record,
  2. groups records whose region sequence is identical,
  3. reports, per taxonomic rank, how many taxa the region can separate,
  4. writes out every ambiguity set: the taxa that collapse into one sequence.

**Why this exists.** Every 16S study picks a region and then reports species names.
Almost none checks whether the chosen region can distinguish the species it is naming.
It usually cannot. Two organisms with an identical sequence over the amplified region are
not hard to classify, they are *impossible* to classify apart, and no classifier, model or
database fixes that. This stage measures it, before sequencing rather than after.

The output is useful two ways round:

  - Choosing a region. Run several primer pairs against the same reference and compare.
    The answer is a count of taxa each region resolves, not an assertion that longer is
    better.
  - Reading a result. The ambiguity sets say which species-level calls from a finished
    run are safe to make, and which are one arbitrary pick out of several equal
    candidates.

**An ambiguity set is a property of the reference and the region, not of your data.** If
a set holds three species and your sample contains one of them, the call is still
ambiguous, because nothing in the amplicon says which. Report the set.

A record that yields no region is reported separately and never silently dropped: a
primer site with a real mismatch can mean the template amplifies poorly, which is a
different fact from being indistinguishable.

Standard library only. Needs Python 3.8 or later. No QIIME 2, no conda, no network.

Example:
    python workflow/06_resolution.py -m ~/ref/HMP_MOCK.v35.fasta \\
        -p v4=GTGCCAGCMGCCGCGGTAA,GGACTACHVGGGTWTCTAAT \\
        -p v3v4=CCTACGGGNGGCWGCAG,GACTACHVGGGTATCTAATCC \\
        --full-length -o ~/research/resolution
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import platform
import shlex
import sys
from datetime import datetime, timezone

__version__ = "0.1.0"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from amplicon_regions import (  # noqa: E402
    RegionError, find_primer, group_by_region, read_fasta, revcomp,
)

# Rank prefixes as SILVA, Greengenes2 and GTDB write them in QIIME 2 taxonomy files.
RANKS = [("d", "domain"), ("p", "phylum"), ("c", "class"), ("o", "order"),
         ("f", "family"), ("g", "genus"), ("s", "species")]
TOP_ALIASES = {"k": "d"}
EMPTY = {"", "unassigned", "unclassified", "uncultured", "metagenome"}

log = logging.getLogger("resolution")

ResolutionError = RegionError


def split_taxon(taxon: str) -> dict[str, str]:
    """Split a taxon string into ranks. Prefixed labels win; otherwise use position."""
    parts = [p.strip() for p in (taxon or "").split(";") if p.strip()]
    by_rank: dict[str, str] = {}
    for i, part in enumerate(parts):
        if len(part) >= 3 and part[1:3] == "__":
            code = TOP_ALIASES.get(part[0], part[0])
            label = part[3:].strip()
        else:
            code = RANKS[i][0] if i < len(RANKS) else None
            label = part
        if code is None:
            continue
        name = dict(RANKS).get(code)
        if name and label.lower() not in EMPTY:
            by_rank[name] = label
    return by_rank


def read_taxonomy(path: str) -> dict[str, dict[str, str]]:
    """Read a two-column id/taxon table, with or without a header."""
    out: dict[str, dict[str, str]] = {}
    with open(path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 2 or not row[0].strip():
                continue
            key = row[0].strip()
            if key in ("Feature ID", "#OTU ID", "id") or key.startswith("#q2:types"):
                continue
            if key in out:
                raise ResolutionError(f"{path}: {key} appears more than once")
            out[key] = split_taxon(row[1])
    if not out:
        raise ResolutionError(f"{path}: no taxonomy rows")
    return out


def parse_primers(values: list[str]) -> dict[str, tuple[str, str]]:
    """Read NAME=FORWARD,REVERSE specifications."""
    out: dict[str, tuple[str, str]] = {}
    for value in values:
        name, sep, rest = value.partition("=")
        if not sep or not name.strip():
            raise ResolutionError(f"--primers wants NAME=FORWARD,REVERSE, got {value!r}")
        fwd, comma, rev = rest.partition(",")
        if not comma or not fwd.strip() or not rev.strip():
            raise ResolutionError(f"--primers wants NAME=FORWARD,REVERSE, got {value!r}")
        name = name.strip()
        if name in out:
            raise ResolutionError(f"region name {name!r} given more than once")
        out[name] = (fwd.strip().upper(), rev.strip().upper())
    return out


def primer_presence(records: dict[str, str], primer: str,
                    max_mismatch: int) -> tuple[int, int]:
    """How many records carry this primer site, as written and as reverse complement."""
    as_written = sum(1 for seq in records.values()
                     if find_primer(seq, primer, max_mismatch))
    rc = revcomp(primer)
    as_rc = sum(1 for seq in records.values() if find_primer(seq, rc, max_mismatch))
    return as_written, as_rc


def explain_no_region(records: dict[str, str], fwd: str, rev: str,
                      max_mismatch: int) -> str:
    """Say which primer site is missing. "No region found" on its own helps nobody.

    A reference trimmed to start inside the amplicon has lost its forward primer site,
    which says nothing about whether the region itself resolves anything.
    """
    n = len(records)
    f_fwd, f_rc = primer_presence(records, fwd, max_mismatch)
    r_fwd, r_rc = primer_presence(records, rev, max_mismatch)
    notes = [f"forward primer present in {f_fwd}/{n} records as written and {f_rc}/{n} "
             f"reverse complemented",
             f"reverse primer present in {r_fwd}/{n} as written and {r_rc}/{n} "
             f"reverse complemented"]
    if f_fwd == 0 and f_rc == 0:
        notes.append("the forward site is absent everywhere, so this reference cannot "
                     "answer for this region. A reference trimmed to start inside the "
                     "amplicon looks exactly like this")
    elif r_fwd == 0 and r_rc == 0:
        notes.append("the reverse site is absent everywhere, same conclusion")
    elif f_rc > f_fwd or r_fwd > r_rc:
        notes.append("the orientation looks reversed: try swapping the primers")
    return "; ".join(notes)


def resolve(groups: dict[str, list[str]], taxonomy: dict[str, dict[str, str]] | None,
            rank: str) -> dict:
    """How many labels at this rank the region separates, and which ones it does not.

    A label is resolved when every group containing it contains no other label. A group
    holding two labels is an ambiguity set: both produce the same sequence here.
    """
    if taxonomy is None:
        return {}
    labelled = 0
    ambiguous_sets = []
    label_in_mixed: set[str] = set()
    all_labels: set[str] = set()
    for seq, names in groups.items():
        labels = {taxonomy[n][rank] for n in names
                  if n in taxonomy and rank in taxonomy[n]}
        all_labels |= labels
        if labels:
            labelled += 1
        if len(labels) > 1:
            label_in_mixed |= labels
            ambiguous_sets.append({"length": len(seq), "records": sorted(names),
                                   "labels": sorted(labels)})
    return {
        "labels": len(all_labels),
        "resolved": len(all_labels - label_in_mixed),
        "unresolved": len(label_in_mixed),
        "resolved_fraction": (len(all_labels - label_in_mixed) / len(all_labels)
                              if all_labels else 0.0),
        "groups_with_a_label": labelled,
        "ambiguity_sets": ambiguous_sets,
    }


def analyse(records: dict[str, str], fwd: str | None, rev: str | None,
            max_mismatch: int, taxonomy: dict | None) -> dict:
    """One region: group the records, then resolve every rank."""
    if fwd is None:
        # Full length: the whole record is the region, so nothing is trimmed away.
        groups: dict[str, list[str]] = {}
        for name, seq in records.items():
            groups.setdefault(seq, []).append(name)
        missing: list[str] = []
    else:
        groups, missing = group_by_region(records, fwd, rev, max_mismatch)
    collapsed = sum(len(v) for v in groups.values() if len(v) > 1)
    lengths = sorted(len(s) for s in groups)
    return {
        "records": len(records),
        "with_region": len(records) - len(missing),
        "no_region": len(missing),
        "no_region_names": sorted(missing),
        "distinct_sequences": len(groups),
        "collapse_groups": sum(1 for v in groups.values() if len(v) > 1),
        "records_collapsed": collapsed,
        "min_length": lengths[0] if lengths else 0,
        "max_length": lengths[-1] if lengths else 0,
        "median_length": lengths[len(lengths) // 2] if lengths else 0,
        "ranks": {rank: resolve(groups, taxonomy, rank) for _, rank in RANKS},
        "groups": groups,
    }


def empty_result(records: dict[str, str]) -> dict:
    """A region that yielded nothing still belongs in the report, as zeroes."""
    return {"records": len(records), "with_region": 0, "no_region": len(records),
            "no_region_names": sorted(records), "distinct_sequences": 0,
            "collapse_groups": 0, "records_collapsed": 0, "min_length": 0,
            "max_length": 0, "median_length": 0,
            "ranks": {rank: {} for _, rank in RANKS}, "groups": {}}


def write_summary(path: str, results: dict[str, dict], with_taxonomy: bool) -> None:
    cols = ["region", "records", "with_region", "no_region", "distinct_sequences",
            "collapse_groups", "records_collapsed", "min_length", "median_length",
            "max_length"]
    if with_taxonomy:
        for _, rank in RANKS:
            cols += [f"{rank}_labels", f"{rank}_resolved", f"{rank}_resolved_fraction"]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for name, r in results.items():
            row = [name] + [r[c] for c in cols[1:10]]
            if with_taxonomy:
                for _, rank in RANKS:
                    d = r["ranks"][rank]
                    row += [d.get("labels", 0), d.get("resolved", 0),
                            round(d.get("resolved_fraction", 0.0), 6)]
            w.writerow(row)


def write_ambiguity(path: str, results: dict[str, dict]) -> None:
    """Every set of taxa that collapse into one sequence. The heart of the report."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["region", "rank", "n_labels", "labels", "n_records", "records",
                    "length"])
        for name, r in results.items():
            for _, rank in RANKS:
                for s in r["ranks"].get(rank, {}).get("ambiguity_sets", []):
                    w.writerow([name, rank, len(s["labels"]), "; ".join(s["labels"]),
                                len(s["records"]), "; ".join(s["records"]), s["length"]])


def write_groups(path: str, results: dict[str, dict]) -> None:
    """Sequence-level collapse groups, with no taxonomy needed to read them."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["region", "n_records", "records", "length"])
        for name, r in results.items():
            for seq, names in sorted(r["groups"].items(), key=lambda kv: -len(kv[1])):
                if len(names) > 1:
                    w.writerow([name, len(names), "; ".join(sorted(names)), len(seq)])


def write_no_region(path: str, results: dict[str, dict]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["region", "record"])
        for name, r in results.items():
            for record in r["no_region_names"]:
                w.writerow([name, record])


def setup_logging(outdir: str | None) -> None:
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(fmt)
    log.addHandler(console)
    if outdir:
        fh = logging.FileHandler(os.path.join(outdir, "resolution_log.txt"), mode="a")
        fh.setFormatter(fmt)
        log.addHandler(fh)


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Measure what an amplicon region can and cannot tell apart.")
    p.add_argument("-m", "--reference", required=True, help="reference FASTA")
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("-p", "--primers", action="append", default=[],
                   metavar="NAME=FWD,REV", help="a region and its primer pair. Repeat "
                                                "to compare several")
    p.add_argument("--full-length", action="store_true",
                   help="also report the whole record as a region, for comparison")
    p.add_argument("-t", "--taxonomy",
                   help="two-column id and taxon string table, to resolve by rank")
    p.add_argument("--max-primer-mismatch", type=int, default=1,
                   help="mismatches allowed per primer site (default 1). Exact sites "
                        "always win")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p.parse_args(argv)


def run(args) -> None:
    reference = os.path.abspath(os.path.expanduser(args.reference))
    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    if not os.path.isfile(reference) or os.path.getsize(reference) == 0:
        raise ResolutionError(f"reference not found or empty: {reference}")
    if args.max_primer_mismatch < 0:
        raise ResolutionError("--max-primer-mismatch cannot be negative, got "
                              f"{args.max_primer_mismatch}")
    regions = parse_primers(args.primers)
    if not regions and not args.full_length:
        raise ResolutionError("give at least one --primers NAME=FWD,REV, or --full-length")
    os.makedirs(outdir, exist_ok=True)
    setup_logging(outdir)
    log.info("=" * 70)
    log.info("resolution %s | Python %s | %s", __version__, platform.python_version(),
             platform.platform())
    log.info("command: %s", shlex.join([sys.executable] + sys.argv))
    log.info("started: %s", datetime.now(timezone.utc).isoformat(timespec="seconds"))

    records = read_fasta(reference)
    log.info("reference: %d records from %s", len(records), reference)
    taxonomy = None
    if args.taxonomy:
        taxonomy = read_taxonomy(os.path.abspath(os.path.expanduser(args.taxonomy)))
        shared = set(records) & set(taxonomy)
        if not shared:
            raise ResolutionError(
                "not one reference record has a taxonomy row. The ids do not match "
                f"(reference has for example {sorted(records)[0]!r}, taxonomy has "
                f"{sorted(taxonomy)[0]!r})")
        if len(shared) < len(records):
            log.warning("%d of %d records have no taxonomy row and are grouped but not "
                        "resolved", len(records) - len(shared), len(records))

    todo = dict(regions)
    if args.full_length:
        todo["full_length"] = (None, None)

    results = {}
    failed = []
    for name, (fwd, rev) in todo.items():
        try:
            results[name] = analyse(records, fwd, rev, args.max_primer_mismatch, taxonomy)
        except RegionError as exc:
            # One region that this reference cannot answer for must not stop the
            # comparison; the whole point is to compare regions side by side.
            log.warning("%s: %s", name, exc)
            log.warning("%s: %s", name, explain_no_region(records, fwd, rev,
                                                          args.max_primer_mismatch))
            results[name] = empty_result(records)
            failed.append(name)
            continue
        r = results[name]
        log.info("%s: %d of %d records yield a region, %d distinct sequences, "
                 "%d group(s) hold more than one record (%d records collapsed), "
                 "length %d to %d",
                 name, r["with_region"], r["records"], r["distinct_sequences"],
                 r["collapse_groups"], r["records_collapsed"], r["min_length"],
                 r["max_length"])
        if r["no_region"]:
            log.info("    %d record(s) yield no region here. A primer site with a real "
                     "mismatch can mean poor amplification, which is a different fact "
                     "from being indistinguishable", r["no_region"])
        if taxonomy:
            for _, rank in RANKS:
                d = r["ranks"][rank]
                if d.get("labels"):
                    log.info("    %-7s %d label(s), %d resolved (%.1f%%), %d lost to "
                             "%d ambiguity set(s)", rank, d["labels"], d["resolved"],
                             100 * d["resolved_fraction"], d["unresolved"],
                             len(d["ambiguity_sets"]))

    if failed and len(failed) == len(todo):
        raise ResolutionError(
            f"no region could be measured: {', '.join(failed)} all failed. See the "
            "per-region messages above for which primer site is missing")
    if failed:
        log.warning("%d of %d region(s) could not be measured against this reference "
                    "(%s) and are reported as zeroes, not as regions that resolve "
                    "nothing", len(failed), len(todo), ", ".join(failed))

    write_summary(os.path.join(outdir, "resolution_summary.tsv"), results,
                  taxonomy is not None)
    write_groups(os.path.join(outdir, "resolution_groups.tsv"), results)
    write_no_region(os.path.join(outdir, "resolution_no_region.tsv"), results)
    if taxonomy:
        write_ambiguity(os.path.join(outdir, "resolution_ambiguity_sets.tsv"), results)

    measured = [n for n in todo if n not in failed]
    if len(measured) > 1 and taxonomy:
        best = max(measured,
                   key=lambda n: results[n]["ranks"]["species"].get("resolved", 0))
        log.info("most species-level labels resolved: %s. That is a count from this "
                 "reference, not a general claim about the region", best)
    log.info("an ambiguity set is a property of the reference and the region, not of "
             "your samples: if a set holds three species and your sample contains one, "
             "the call is still ambiguous, because nothing in the amplicon says which")
    log.info("done: %s", os.path.join(outdir, "resolution_summary.tsv"))


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        run(args)
    except ResolutionError as exc:
        if not log.handlers:
            setup_logging(None)
        log.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
